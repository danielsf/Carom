/*
spock:

maybe you should just set burn_in

use one interval of that many  steps to find the eigen vectors

then use another interval of that (sub_divided) to set _factor
--yourself
*/

#include "mcmc/mcmc.h"

mcmc::~mcmc(){
   if(_dice!=NULL){
       delete _dice;
   }
}

void mcmc::set_min(int dex, double dd){
    _guess_min.set(dex,dd);
}

void mcmc::set_max(int dex, double dd){
    _guess_max.set(dex,dd);
}

void mcmc::initialize(){
    _bases.set_name("mcmc_bases");
    _sigma.set_name("mcmc_sigma");
    _guess_min.set_name("mcmc_guess_min");
    _guess_max.set_name("mcmc_guess_max");
    _chisq = NULL;
    _burn_in=0;
    _dice=NULL;
    _factor=2.38;
    _time_started=double(time(NULL));

    sprintf(_name_root,"chain");

}

mcmc::mcmc(int nchains, int seed, chisquared *fn){
    initialize();

    if(seed<0){
        _dice=new Ran(int(time(NULL)));
    }
    else{
        _dice=new Ran(seed);
    }

    _chisq=fn;
    _chains.initialize(nchains,_chisq->get_dim(),_dice);

    _bases.set_cols(_chisq->get_dim());

}

void mcmc::write_timing(char *msg){
   char name[2*letters];
    sprintf(name,"%s_timing.txt",_name_root);

    FILE *output;

    output=fopen(name,"a");
    fprintf(output,"at %d %s\n",_chisq->get_called(),msg);
    fclose(output);

}

void mcmc::write_timing(int overwrite){
    char name[2*letters];
    sprintf(name,"%s_timing.txt",_name_root);

    FILE *output;

    if(overwrite==1){
        output=fopen(name,"w");
        fprintf(output,"#calls time timeper timeperRaw overhead acceptance factor thinby totalpts\n");
    }
    else{
        output=fopen(name,"a");
    }

    double overhead,timeSpent,timePer,timePerRaw;

    timeSpent=double(time(NULL))-_time_started;
    timePer=timeSpent/double(_chisq->get_called());

    timePerRaw=_chisq->get_time_spent()/double(_chisq->get_called());
    overhead = timePer-timePerRaw;

    fprintf(output,"%d %e %e %e %e %e %e %d %d\n",
    _chisq->get_called(),timeSpent,timePer,timePerRaw,overhead,acceptance_rate(),
    _factor,_chains.get_thinby(0.1,0.0),_chains.get_points());

    fclose(output);

}

void mcmc::set_burnin(int bb){
    //_check_every will be however often we assess bases
    _burn_in=bb;
}

void mcmc::set_name_root(char *word){
    int i;
    for(i=0;i<letters-10 && word[i]!=0;i++){
        _name_root[i]=word[i];
    }
    _name_root[i]=0;

}

double mcmc::acceptance_rate(){
    int ichain,ipt,ct,rows;
    ct=0;
    rows=0;
    for(ichain=0;ichain<_chains.get_n_chains();ichain++){
        for(ipt=0;ipt<_chains(ichain)->get_rows();ipt++){
            rows++;
            ct+=_chains(ichain)->get_degeneracy(ipt);
        }
    }

    if(ct==0){
        return 0.0;
    }

    return double(rows)/double(ct);

}

double mcmc::update_bases(){
    array_2d<double> covar;
    covar.set_name("mcmc_update_bases_covar");

    array_2d<double> old_basis,evecs;
    array_1d<double> evals;

    old_basis.set_name("mcmc_update_bases_old_basis");
    evecs.set_name("mcmc_update_bases_evecs");
    evals.set_name("mcmc_update_bases_evals");

    try{
        _chains.get_covariance_matrix(covar);

        eval_symm(covar,evecs,evals,0.1);
    }
    catch(int iex){
        return -1.0;
    }

    old_basis.set_cols(_chisq->get_dim());

    int ix,iy;
    for(ix=0;ix<_chisq->get_dim();ix++){
        for(iy=0;iy<_chisq->get_dim();iy++){
            old_basis.set(ix,iy,_bases.get_data(ix,iy));
            _bases.set(ix,iy,evecs.get_data(ix,iy));
        }
        _bases(ix)->normalize();
        _sigma.set(ix,sqrt(fabs(evals.get_data(ix))));
    }

    double dotMax,dot,dotMin;
    int i;
    dotMax=-1.0;
    for(ix=0;ix<_chisq->get_dim();ix++){
        dotMin=2.0*exception_value;
        for(iy=0;iy<_chisq->get_dim();iy++){
            dot=0.0;
            for(i=0;i<_chisq->get_dim();i++){
                dot+=old_basis.get_data(ix,i)*_bases.get_data(iy,i);
            }

            if(fabs(1.0-fabs(dot))<dotMin){
                dotMin=fabs(1.0-fabs(dot));
            }
        }

        if(dotMin>dotMax){
            dotMax=dotMin;
        }
    }

    validate_bases();
    return dotMax;
}

void mcmc::sample(int nSamples){

    write_timing(1);

    int i,j;
    if(_bases.get_rows()==0){
        for(i=0;i<_chisq->get_dim();i++){
            _sigma.set(i,0.1*(_guess_max.get_data(i)-_guess_min.get_data(i)));
            for(j=0;j<_chisq->get_dim();j++){
                if(i==j){
                    _bases.set(i,j,1.0);
                }
                else{
                    _bases.set(i,j,0.0);
                }
            }
        }
    }

    validate_bases();

    array_1d<double> trial,dir;
    trial.set_name("mcmc_sample_trial");
    dir.set_name("mcmc_sample_dir");

    int iChain,ix,iy,thinby,total_points,update_ct;
    int final_ct,burn_ct,total_ct,last_wrote,last_assessed;
    int last_updated,last_updated_factor,something_changed;
    int last_dumped;

    double sqrtD=sqrt(_chisq->get_dim());
    double mu,acceptance,norm;

    char message[3*letters];

    for(iChain=0;iChain<_chains.get_n_chains();iChain++){
        _chains(iChain)->set_output_name_root(_name_root);
        _chains(iChain)->set_chain_label(iChain);
    }

    final_ct=0;
    burn_ct=0;
    total_ct=0;
    last_updated=0;
    last_updated_factor=0;
    last_wrote=0;
    last_dumped=0;
    update_ct=0;

    while(final_ct<nSamples){
        for(iChain=0;iChain<_chains.get_n_chains();iChain++){
            mu=2.0*exception_value;
            while(mu>exception_value){
                if(_chains(iChain)->is_current_point_valid()==1){
                    for(ix=0;ix<_chisq->get_dim();ix++){
                        trial.set(ix,_chains(iChain)->get_current_point(ix));
                        dir.set(ix,normal_deviate(_dice,0.0,1.0));
                    }
                    norm=_factor/sqrtD;

                    for(ix=0;ix<_chisq->get_dim();ix++){
                        for(iy=0;iy<_chisq->get_dim();iy++){
                            trial.add_val(iy,norm*dir.get_data(ix)*_sigma.get_data(ix)*_bases.get_data(ix,iy));
                        }
                    }
                }
                else{
                    for(ix=0;ix<_chisq->get_dim();ix++){
                        trial.set(ix,_guess_min.get_data(ix)+_dice->doub()*(_guess_max.get_data(ix)-_guess_min.get_data(ix)));
                    }
                }

                mu=_chisq[0](trial);
            }
            _chains(iChain)->add_point(trial,mu);

        }
        total_ct++;

        something_changed=0;
        final_ct++;

        if(final_ct<_burn_in){
            if(final_ct>=last_updated+500){

                mu=update_bases();
                sprintf(message,"updating bases -- dotMax %e;",mu);
                write_timing(message);
                last_updated=final_ct;
                last_updated_factor=final_ct;
                last_dumped=final_ct;
                last_wrote=final_ct;
                write_timing(0);
                for(iChain=0;iChain<_chains.get_n_chains();iChain++){
                    _chains(iChain)->write_chain(1);
                }
            }
            else if(final_ct>=last_updated_factor+100){
                acceptance=acceptance_rate();
                if(fabs(1.0/acceptance-3.0)>1.0){
                    if(acceptance<0.3333333){
                         _factor*=0.9;
                    }
                    else{
                        _factor*=1.2;
                    }
                }
                last_updated_factor=final_ct;

            }
        }
        else if(final_ct>last_wrote+1000){
            if(final_ct>last_dumped+5000){
                ix=1;
                last_dumped=final_ct;
            }
            else{
                ix=0;
            }

            write_timing(0);
            for(iChain=0;iChain<_chains.get_n_chains();iChain++){
                _chains(iChain)->write_chain(ix);
            }
            last_wrote=final_ct;
        }

        /*if(something_changed==1){
            for(iChain=0;iChain<_chains.get_n_chains();iChain++){
                _chains(iChain)->write_chain(1);
                _chains(iChain)->increment_iteration();
            }
            last_wrote=0;
            last_assessed=0;
            final_ct=0;
        }*/
   }

    printf("done %d %d -- %d\n",final_ct,nSamples,_chisq->get_called());

   for(iChain=0;iChain<_chains.get_n_chains();iChain++){
       _chains(iChain)->write_chain(1);
   }

}


void mcmc::find_fisher_matrix(array_2d<double> &covar, array_1d<double> &centerOut, double *minOut){

    array_1d<double> trial,norm;
    trial.set_name("node_findCovar_trial");
    norm.set_name("node_findCovar_norm");

    array_1d<double> center;
    array_2d<double> seed;
    double fcenter;
    center.set_name("mcmc_find_fisher_center");
    seed.set_name("mcmc_find_fisher_seed");

    int i,j;
    array_1d<double> min,max;
    seed.set_cols(_chisq->get_dim());

    min.set_name("mcmc_find_fisher_min");
    max.set_name("mcmc_find_fisher_max");

    for(i=0;i<_chisq->get_dim();i++){
        if(_chisq->get_min(i)<exception_value){
            min.set(i,_chisq->get_min(i));
        }
        else{
            min.set(i,_guess_min.get_data(i));
        }

        if(_chisq->get_max(i)>-1.0*exception_value){
            max.set(i,_chisq->get_max(i));
        }
        else{
            max.set(i,_guess_max.get_data(i));
        }

        for(j=0;j<_chisq->get_dim()+1;j++){
            seed.set(j,i,min.get_data(i)+_dice->doub()*(max.get_data(i)-min.get_data(i)));
        }
    }

    simplex_minimizer f_min;
    f_min.set_chisquared(_chisq);
    f_min.set_minmax(min,max);
    f_min.set_dice(_dice);
    f_min.use_gradient();
    f_min.find_minimum(seed, center);
    fcenter=f_min.get_minimum();
    minOut[0]=fcenter;
    for(i=0;i<_chisq->get_dim();i++){
        centerOut.set(i,center.get_data(i));
    }


    double sgn;
    array_1d<double> temp_dir,temp_pt;
    temp_dir.set_name("mcmc_find_fisher_temp_dir");
    temp_pt.set_name("mcmc_find_fisher_temp_pt");

    for(i=0;i<_chisq->get_dim();i++){
        temp_dir.set(i,0.0);
    }

    for(i=0;i<_chisq->get_dim();i++){
        for(sgn=-1.0;sgn<1.1;sgn+=2.0){
            temp_dir.set(i,sgn);
            mcmc_bisection(-1.0,center,fcenter,temp_dir,temp_pt);

            if(sgn<0.0){
                min.set(i,temp_pt.get_data(i));
            }
            else{
                max.set(i,temp_pt.get_data(i));
            }

        }
        temp_dir.set(i,0.0);
    }


    for(i=0;i<_chisq->get_dim();i++){
        norm.set(i,max.get_data(i)-min.get_data(i));
    }

    array_2d<double> fpp,fpm,fmp,fmm;
    fpp.set_name("mcmc_find_fisher_fpp");
    fpm.set_name("mcmc_find_fisher_fpm");
    fmp.set_name("mcmc_find_fisher_fmp");
    fmm.set_name("mcmc_find_fisher_fmm");

    fpp.set_cols(_chisq->get_dim());
    fpm.set_cols(_chisq->get_dim());
    fmp.set_cols(_chisq->get_dim());
    fmm.set_cols(_chisq->get_dim());

    array_1d<double> f2p,f2m;
    f2p.set_name("mcmc_find_fisher_f2p");
    f2m.set_name("mcmc_find_fisher_f2m");

    double mu;
    array_1d<double> dx;
    dx.set_name("mcmc_find_fisher_dx");
    for(i=0;i<_chisq->get_dim();i++){
        dx.set(i,1.0e-2);
    }

    int ix,iy,keepGoing,ctAbort,ctAbortMax,calledMax;
    int ibefore=_chisq->get_called();

    double tol=-1.0;
    double dCenter;

    int abortValue,abortProx;

    abortValue=0;
    abortProx=0;

    ctAbort=0;
    ctAbortMax=100;
    calledMax = 10*_chisq->get_dim()*_chisq->get_dim();

    for(i=0;i<_chisq->get_dim();i++){
        trial.set(i,center.get_data(i));
    }

    for(ix=0;ix<_chisq->get_dim();ix++){
        for(i=0;i<_chisq->get_dim();i++){
            trial.set(i,center.get_data(i));
        }

        if(_chisq->get_called()-ibefore>calledMax ||
           ctAbort>=ctAbortMax){
                printf("Could not find CoVar; aborting\n");
                printf("ctAbort %d -- %d %d\n",ctAbort,abortValue,abortProx);
                printf("called %d\n",_chisq->get_called()-ibefore);
                throw -1;
        }

        if(fcenter>exception_value){
            printf("Center is invalid; aborting\n");
            throw -1;
        }

        keepGoing=1;
        trial.set(ix,center.get_data(ix)+2.0*dx.get_data(ix)*norm.get_data(ix));
        mu=_chisq[0](trial);
        dCenter=fabs(mu-fcenter);

        if(mu<exception_value && dCenter>tol){
            f2p.set(ix,mu);
        }
        else if(mu<exception_value){
            dx.multiply_val(ix,1.5);
            keepGoing=0;
            ix--;
            ctAbort++;
            abortValue++;
        }
        else{
            center.subtract_val(ix,2.5*dx.get_data(ix)*norm.get_data(ix));
            fcenter=_chisq[0](center);
            keepGoing=0;
            ix--;
            ctAbort++;
            abortProx++;
        }

        if(keepGoing==1){
            trial.set(ix,center.get_data(ix)-2.0*dx.get_data(ix)*norm.get_data(ix));
            mu=_chisq[0](trial);
            dCenter=fabs(mu-fcenter);

            if(mu<exception_value && dCenter>tol){
                f2m.set(ix,mu);
            }
            else if(mu<exception_value){
                dx.multiply_val(ix,1.5);
                keepGoing=0;
                ix--;
                ctAbort++;
                abortValue++;
            }
            else{
                center.add_val(ix,2.5*dx.get_data(ix)*norm.get_data(ix));
                fcenter=_chisq[0](center);
                keepGoing=0;
                ix--;
                ctAbort++;
                abortProx++;
            }
        }

        for(iy=ix-1;iy>=0 && keepGoing==1;iy--){
            for(i=0;i<_chisq->get_dim();i++){
                trial.set(i,center.get_data(i));
            }

            if(_chisq->get_called()-ibefore>calledMax ||
               ctAbort>=ctAbortMax){
                printf("Could not find CoVar; aborting\n");
                printf("ctAbort %d -- %d %d\n",ctAbort,abortValue,abortProx);
                printf("called %d\n",_chisq->get_called()-ibefore);
                throw -1;
            }

            if(fcenter>exception_value){
                printf("center is invalid; aborting\n");
                throw -1;
            }

            trial.set(ix,center.get_data(ix)+dx.get_data(ix)*norm.get_data(ix));
            trial.set(iy,center.get_data(iy)+dx.get_data(iy)*norm.get_data(iy));
            mu=_chisq[0](trial);
            dCenter=fabs(mu-fcenter);
            if(mu<exception_value && dCenter>tol){
                fpp.set(ix,iy,mu);
            }
            else if(mu<exception_value){
                dx.multiply_val(ix,1.5);
                dx.multiply_val(iy,1.5);
                keepGoing=0;
                ix--;
                ctAbort++;
                abortValue++;
            }
            else{
                center.subtract_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                center.subtract_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                fcenter=_chisq[0](center);
                keepGoing=0;
                ix--;
                ctAbort++;
                abortProx++;
            }

            if(keepGoing==1){
               trial.set(iy,center.get_data(iy)-dx.get_data(iy)*norm.get_data(iy));
               mu=_chisq[0](trial);
               dCenter=fabs(mu-fcenter);
               if(mu<exception_value && dCenter>tol){
                   fpm.set(ix,iy,mu);
               }
               else if(mu<exception_value){
                   dx.multiply_val(ix,1.5);
                   dx.multiply_val(iy,1.5);
                   keepGoing=0;
                   ix--;
                   ctAbort++;
                   abortValue++;
               }
               else{
                   center.subtract_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                   center.add_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                   fcenter=_chisq[0](center);
                   keepGoing=0;
                   ix--;
                   ctAbort++;
                   abortProx++;
               }
            }

            if(keepGoing==1){
                trial.set(ix,center.get_data(ix)-dx.get_data(ix)*norm.get_data(ix));
                mu=_chisq[0](center);
                dCenter=fabs(mu-fcenter);
                if(mu<exception_value && dCenter>tol){
                    fmm.set(ix,iy,mu);
                }
                else if(mu<exception_value){
                    dx.multiply_val(ix,1.5);
                    dx.multiply_val(iy,1.5);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                    abortValue++;
                }
                else{
                    center.add_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                    center.add_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                    fcenter=_chisq[0](center);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                    abortProx++;
                }
            }

            if(keepGoing==1){
                trial.set(iy,center.get_data(iy)+dx.get_data(iy)*norm.get_data(iy));
                mu=_chisq[0](trial);
                dCenter=fabs(mu-fcenter);
                if(mu<exception_value && dCenter>tol){
                    fmp.set(ix,iy,mu);
                }
                else if(mu<exception_value){
                    dx.multiply_val(ix,1.5);
                    dx.multiply_val(iy,1.5);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                    abortValue++;
                }
                else{
                    center.add_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                    center.subtract_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                    fcenter=_chisq[0](center);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                    abortProx++;
                }
            }

        }
    }

    covar.set_cols(_chisq->get_dim());
    for(ix=0;ix<_chisq->get_dim();ix++){
        covar.set(ix,ix,0.25*(f2p.get_data(ix)+f2m.get_data(ix)-2.0*fcenter)/power(dx.get_data(ix)*norm.get_data(ix),2));
    }

    double num,denom;
    for(ix=0;ix<_chisq->get_dim();ix++){
        for(iy=ix-1;iy>=0;iy--){
            num=0.25*(fpp.get_data(ix,iy)+fmm.get_data(ix,iy)-fmp.get_data(ix,iy)-fpm.get_data(ix,iy));
            denom=dx.get_data(ix)*norm.get_data(ix)*dx.get_data(iy)*norm.get_data(iy);
            covar.set(ix,iy,num/denom);
            covar.set(iy,ix,num/denom);
        }
    }

    printf("successfully got covar -- %e %d\n",fcenter,_chisq->get_called());

}


void mcmc::find_fisher_eigen(array_2d<double> &bases, array_1d<double> &centerOut, double *minVal){

    array_2d<double> covar;
    covar.set_name("mcmc_guess_bases_covar");
    covar.set_cols(_chisq->get_dim());
    int ix,iy;
    double covarmax=-1.0;

    find_fisher_matrix(covar,centerOut,minVal);

    for(ix=0;ix<_chisq->get_dim();ix++){
        for(iy=ix;iy<_chisq->get_dim();iy++){
            if(fabs(covar.get_data(ix,iy))>covarmax){
                covarmax=fabs(covar.get_data(ix,iy));
            }
        }
    }

    for(ix=0;ix<_chisq->get_dim();ix++){
        for(iy=0;iy<_chisq->get_dim();iy++){
            covar.divide_val(ix,iy,covarmax);
        }
    }

    printf("assembled covariance matrix\n");

    array_2d<double> evecs;
    evecs.set_name("mcmc_guess_bases_evecs");
    array_1d<double> evals;
    evals.set_name("mcmc_guess_bases_evals");

    evecs.set_cols(_chisq->get_dim());

    try{
        eval_symm(covar,evecs,evals,0.001);
    }
    catch(int iex){
        printf("Guess failed on eigen vectors\n");
        throw -1;
    }

    for(ix=0;ix<_chisq->get_dim();ix++){
        for(iy=0;iy<_chisq->get_dim();iy++){
            bases.set(ix,iy,evecs.get_data(ix,iy));
        }
        bases(ix)->normalize();
    }

    printf("got first batch of guessed bases\n");

}

void mcmc::mcmc_bisection(double deltachi, array_1d<double> &lowball_in, double flow, array_1d<double> &dir, array_1d<double> &found){

    double target=flow+deltachi;

    if(target<flow){
        target=flow+double(_chisq->get_dim());
    }

    array_1d<double> highball,lowball;
    highball.set_name("mcmc_bisection_highball");
    lowball.set_name("mcmc_bisection_lowball");
    int i;
    for(i=0;i<_chisq->get_dim();i++){
        lowball.set(i,lowball_in.get_data(i));
        highball.set(i,lowball.get_data(i));
    }

    double fhigh;
    fhigh=-2.0*exception_value;
    while(fhigh<target){
        for(i=0;i<_chisq->get_dim();i++){
            highball.add_val(i,dir.get_data(i));
        }
        fhigh=_chisq[0](highball);
    }

    array_1d<double> trial;
    trial.set_name("mcmc_bisection_trial");
    double fbest=flow;
    int ct=0;

    double tol=0.1*double(_chisq->get_dim());
    double ftrial;

    for(i=0;i<_chisq->get_dim();i++){
        found.set(i,lowball.get_data(i));
    }

    while(ct<20 && fabs(fbest-target)>tol){
        ct++;
        for(i=0;i<_chisq->get_dim();i++){
            trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
        }
        ftrial=_chisq[0](trial);

        if(ftrial<target){
            for(i=0;i<_chisq->get_dim();i++){
                lowball.set(i,trial.get_data(i));
            }
        }
        else{
            for(i=0;i<_chisq->get_dim();i++){
                highball.set(i,trial.get_data(i));
            }
        }

        if(fabs(ftrial-target)<fabs(fbest-target)){
            fbest=ftrial;
            for(i=0;i<_chisq->get_dim();i++){
                found.set(i,trial.get_data(i));
            }
        }
    }

}

void mcmc::guess_bases(double deltaChi, int seedPoints){

    array_2d<double> temp_bases;
    temp_bases.set_name("mcmc_guess_bases_temp_bases");
    temp_bases.set_cols(_chisq->get_dim());

    int i,j;
    array_1d<double> center;
    double minVal;

    center.set_name("mcmc_guess_bases_center");


    try{
        find_fisher_eigen(temp_bases,center,&minVal);
        for(i=0;i<_chisq->get_dim();i++){
            for(j=0;j<_chisq->get_dim();j++){
                _bases.set(i,j,temp_bases.get_data(i,j));
            }
        }
    }
    catch(int iex){
        printf("guessing bases failed %d\n",_chisq->get_called());
        return;
    }

    array_1d<double> temp_dir,temp_pt;
    temp_dir.set_name("mcmc_guess_bases_temp_dir");
    temp_pt.set_name("mcmc_guess_bases_temp_pt");

    int ix;
    double sgn,dd,d;
    for(ix=0;ix<_chisq->get_dim();ix++){
        dd=0.0;
        for(sgn=-1.0;sgn<1.1;sgn+=2.0){
            for(i=0;i<_chisq->get_dim();i++){
                temp_dir.set(i,sgn*_bases.get_data(ix,i));
            }
            mcmc_bisection(deltaChi, center,minVal,temp_dir,temp_pt);

            d=0.0;
            for(i=0;i<_chisq->get_dim();i++){
                d+=power(center.get_data(i)-temp_pt.get_data(i),2);
            }
            dd+=sqrt(d);

        }
        _sigma.set(ix,dd);
    }

    printf("finished guessing bases; %d\n",_chisq->get_called());
    for(i=0;i<_chisq->get_dim();i++){
        printf("%.3e -- ",_sigma.get_data(i));
        for(j=0;j<_chisq->get_dim();j++)printf("%.3e ",_bases.get_data(i,j));
        printf("\n");
    }
    printf("\n");

    int k;
    double mu,norm;

    array_1d<double> dir;
    dir.set_name("mcmc_guess_bases_dir");

    if(seedPoints==1){
        temp_pt.reset();
        for(i=0;i<_chains.get_n_chains();i++){
            mu=2.0*exception_value;
            while(mu>exception_value){
                norm=fabs(normal_deviate(_dice,0.0,2.0))+2.0;
                for(j=0;j<_chisq->get_dim();j++){
                    temp_pt.set(j,center.get_data(j));
                    dir.set(j,normal_deviate(_dice,0.0,1.0));
                }
                dir.normalize();

                for(j=0;j<_chisq->get_dim();j++){
                    for(k=0;k<_chisq->get_dim();k++){
                        temp_pt.add_val(k,norm*dir.get_data(j)*_sigma.get_data(j)*_bases.get_data(j,k));
                    }
                }

                mu=_chisq[0](temp_pt);
            }
            _chains(i)->add_point(temp_pt,mu);

        }
    }

}

void mcmc::validate_bases(){
    array_1d<double> temp;
    temp.set_name("mcmc_validate_bases");

    int ix,iy,i;
    double err,maxerr,dot,tol=0.01;
    int isOrth,maxIsOrth;
    for(ix=0;ix<_chisq->get_dim();ix++){
        _bases(ix)->normalize();
    }

    maxerr=-1.0;
    for(ix=0;ix<_chisq->get_dim();ix++){
        for(iy=ix;iy<_chisq->get_dim();iy++){
            dot=0.0;
            for(i=0;i<_chisq->get_dim();i++){
                dot+=_bases.get_data(ix,i)*_bases.get_data(iy,i);
            }

            if(ix==iy){
                isOrth=0;
                err=fabs(1.0-sqrt(dot));
            }
            else{
                isOrth=1;
                err=sqrt(fabs(dot));
            }

            if(err>maxerr){
                maxIsOrth=isOrth;
                maxerr=err;
            }

        }
    }

    if(maxerr>tol){
        printf("WARNING bases max err %e\n",maxerr);
        printf("isOrth %d\n",maxIsOrth);
        exit(1);
    }
}
