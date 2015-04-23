#include "mcmc/mcmc.h"

mcmc::~mcmc(){
   if(dice!=NULL){
       delete dice;
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
    _centerpt.set_name("mcmc_centerpt");
    _chisq = NULL;
    _burn_in=1000;
    _last_set=0;
    _max_covar=-2.0*exception_value;
    _covar_lim=0.1;
    dice=NULL;
    
    sprintf(_name_root,"chain");
}

mcmc::mcmc(int nchains, chisquared *fn){
    initialize();
    _chisq=fn;
    _chains.initialize(nchains,_chisq->get_dim());
    
    _bases.set_cols(_chisq->get_dim());
    int i,j;
    for(i=0;i<_chisq->get_dim();i++){
        _sigma.set(i,1.0);
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

void mcmc::set_burn_in(int ii){
    //_burn_in will be however often we assess bases
    _burn_in=ii;
}

void mcmc::set_covar_lim(double dd){
    //_covar_lim will be the value of _max_covar at which we are satisfied
    //with the bases
    _covar_lim=dd;
}

void mcmc::set_name_root(char *word){
    int i;
    for(i=0;i<letters-10 && word[i]!=0;i++){
        _name_root[i]=word[i];
    }
    _name_root[i]=0;
}

void mcmc::find_fisher_matrix(array_2d<double> &covar, array_1d<double> &centerOut){

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
        
        for(j=0;j<_chisq->get_dim();j++){
            seed.set(j,i,min.get_data(i)+dice->doub()*(max.get_data(i)-min.get_data(i)));
        }
    }
    
    simplex_minimizer f_min;
    f_min.set_chisquared(_chisq);
    f_min.find_minimum(seed, center);
    fcenter=f_min.get_minimum();
    _min_val=fcenter;
    
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
    
    double tol=1.0e-10;
    double dCenter;
    
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
                printf("ctAbort %d\n",ctAbort);
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
        }
        else{
            center.subtract_val(ix,2.5*dx.get_data(ix)*norm.get_data(ix));
            fcenter=_chisq[0](center);
            keepGoing=0;
            ix--;
            ctAbort++;
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
            }
            else{
                center.add_val(ix,2.5*dx.get_data(ix)*norm.get_data(ix));
                fcenter=_chisq[0](center);
                keepGoing=0;
                ix--;
                ctAbort++;
            }
        }
        
        for(iy=ix-1;iy>=0 && keepGoing==1;iy--){
            for(i=0;i<_chisq->get_dim();i++){
                trial.set(i,center.get_data(i));
            }
            
            if(_chisq->get_called()-ibefore>calledMax ||
               ctAbort>=ctAbortMax){
                printf("Could not find CoVar; aborting\n");
                printf("ctAbort %d\n",ctAbort);
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
            }
            else{
                center.subtract_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                center.subtract_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                fcenter=_chisq[0](center);
                keepGoing=0;
                ix--;
                ctAbort++;
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
               }
               else{
                   center.subtract_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                   center.add_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                   fcenter=_chisq[0](center);
                   keepGoing=0;
                   ix--;
                   ctAbort++;
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
                }
                else{
                    center.add_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                    center.add_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                    fcenter=_chisq[0](center);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
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
                }
                else{
                    center.add_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                    center.subtract_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                    fcenter=_chisq[0](center);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
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
    for(i=0;i<_chisq->get_dim();i++){
        centerOut.set(i,center.get_data(i));
    }
}


void mcmc::find_fisher_eigen(array_2d<double> &bases){
    
    array_2d<double> covar;
    covar.set_name("mcmc_guess_bases_covar");
    covar.set_cols(_chisq->get_dim());
    int ix,iy;
    double covarmax=-1.0;
    
    find_fisher_matrix(covar,_centerpt);
    
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
    
    int i1=_chisq->get_dim()/2;
    try{
        eval_symm(covar,evecs,evals,i1,_chisq->get_dim(),1,0.001);
    }
    catch(int iex){
        printf("Guess failed on first batch of eigen vectors\n");
        throw -1;
    }
    
    for(ix=0;ix<i1;ix++){
        for(iy=0;iy<_chisq->get_dim();iy++){
            bases.set(ix,iy,evecs.get_data(iy,ix));
        }
        bases(ix)->normalize();
    }
    
    printf("got first batch of guessed bases\n");
    
    try{
        eval_symm(covar,evecs,evals,_chisq->get_dim()-i1,_chisq->get_dim(),-1,0.001);
    }
    catch(int iex){
        printf("Guess failed on second batch of eigen vectors\n");
        throw -1;
    }
    
    for(ix=i1;ix<_chisq->get_dim();ix++){
        for(iy=0;iy<_chisq->get_dim();iy++){
            bases.set(ix,iy,evecs.get_data(iy,ix-i1));
        }
        bases(ix)->normalize();
    }
    
    printf("got second batch of guessed bases\n");
    
}

void mcmc::bisection(array_1d<double> &lowball_in, double flow, array_1d<double> &dir, array_1d<double> &found){

    double target=flow+double(_chisq->get_dim());
    
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
        for(i=0;i<_chisq->get_dim();i++){
            trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
        }
        ftrial=_chisq[0](trial);
        
        if(ftrial<flow){
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

void mcmc::guess_bases(){

    array_2d<double> temp_bases;
    temp_bases.set_name("mcmc_guess_bases_temp_bases");
    
    int i,j;
    
    try{
        find_fisher_eigen(temp_bases);
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
            bisection(_centerpt,_min_val,temp_dir,temp_pt);
            
            d=0.0;
            for(i=0;i<_chisq->get_dim();i++){
                d+=power(_centerpt.get_data(i)-temp_pt.get_data(i),2);
            }
            dd+=sqrt(d);
            
        }
        _sigma.set(ix,dd/6.0);
    }
    
    printf("finished guessing bases; %d\n",_chisq->get_called());
} 
