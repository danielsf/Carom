#include "dalex.h"

void dalex::build(chisq_wrapper *cc){

    _chifn=cc;

    _explorers.set_chifn(_chifn);
    _min_explorers.set_chifn(_chifn);

    int i,j;
    _basis_vectors.set_cols(_chifn->get_dim());
    for(i=0;i<_chifn->get_dim();i++){
        _basis_norm.set(i,_chifn->get_characteristic_length(i));
        for(j=0;j<_chifn->get_dim();j++){
            if(i==j){
                _basis_vectors.set(i,j,1.0);
            }
            else{
                _basis_vectors.set(i,j,0.0);
            }
        }
    }

}

double dalex::get_norm(int dex){

    if(_good_points.get_dim()==0){
        return _chifn->get_characteristic_length(dex);
    }

    double min,max;
    int i,ip;
    for(i=0;i<_good_points.get_dim();i++){

        ip=_good_points.get_data(i);

        if(i==0 || _chifn->get_pt(ip,dex)<min){
            min=_chifn->get_pt(ip,dex);
        }

        if(i==0 || _chifn->get_pt(ip,dex)>max){
            max=_chifn->get_pt(ip,dex);
        }
    }

    if(max-min<1.0e-20){
        return _chifn->get_characteristic_length(dex);
    }

    return max-min;
}

void dalex::search(){
    safety_check("search");
    int i,j,i_found;
    double mu;
    int pts_0=_chifn->get_pts();
    array_1d<double> pt;
    array_1d<int> to_use,to_kick;
    assess_good_points();

    int has_explored=0;

    _chifn->set_search_type(_type_refine);
    iterate_on_minimum();

    _chifn->set_search_type(_type_explore);
    int is_outside;
    int need_kick=1;
    while(to_use.get_dim()==0 && (_limit<0 || _chifn->get_pts()<_limit)){
        explore(need_kick);
        need_kick=0;
        for(i=0;i<_explorers.get_n_particles();i++){
            _explorers.get_pt(i,pt);
            evaluate(pt,&mu,&i_found);
            if(mu<target()){
                to_kick.add(i);
                is_outside=1;
                for(j=0;j<_exclusion_zones.ct() && is_outside==1;j++){
                    if(_exclusion_zones(j)->contains(pt)==1){
                        is_outside=0;
                    }
                }
                if(is_outside==1){
                    to_use.add(i_found);
                }
            }
        }
    }

    array_1d<double> dd_min,dd_min_sorted;
    for(i=0;i<to_use.get_dim();i++){
        dd_min.set(i,0.0);
        for(j=0;j<_chifn->get_dim();j++){
            mu=_chifn->get_pt(to_use.get_data(i),j)-_chifn->get_pt(mindex(),j);
            dd_min.add_val(i,power(mu/_chifn->get_characteristic_length(j),2));
        }
    }

    sort(dd_min, dd_min_sorted, to_use);

    int i_end;

    _chifn->set_search_type(_type_tendril);
    for(i=0;i<to_use.get_dim();i++){
        is_outside=1;
        for(j=0;j<_exclusion_zones.ct() && is_outside==1;j++){
            if(_exclusion_zones(j)->contains(_chifn->get_pt(to_use.get_data(i)))==1){
                is_outside=0;
            }
        }
        if(is_outside==1){
            tendril_search(to_use.get_data(i));
        }
        if(_limit>0 && _chifn->get_pts()>_limit){
            break;
        }
    }
    _update_good_points(pts_0);

}


void dalex::simplex_search(){
    array_1d<int> empty;
    simplex_search(empty);
}

void dalex::simplex_search(int ii){
    array_1d<int> specified;
    specified.add(ii);
    simplex_search(specified);
}

void dalex::simplex_search(array_1d<int> &specified){

    printf("\n    doing dalex_simplex %e -- %e %e\n",chimin(),target(),_target_factor);
    safety_check("simplex_search");
    array_1d<double> min,max;
    min.set_name("dalex_simplex_min");
    max.set_name("dalex_simplex_max");
    array_2d<double> seed;
    seed.set_name("dalex_simplex_seed");
    array_1d<double> trial,grad;

    int i,j;
    int pt_0=_chifn->get_pts();

    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,0.0);
        max.set(i,_chifn->get_characteristic_length(i));
    }

    array_1d<int> chosen_seed;

    for(i=0;i<specified.get_dim();i++){
        seed.add_row(_chifn->get_pt(specified.get_data(i)));
    }

    while(seed.get_rows()<_chifn->get_dim()+1){
        i=_chifn->random_int()%_min_explorers.get_n_particles();
        if(chosen_seed.contains(i)==0){
            _min_explorers.get_pt(i,trial);
            seed.add_row(trial);
            chosen_seed.add(i);
        }
    }

    simplex_minimizer ffmin;
    ffmin.set_chisquared(_chifn);
    ffmin.set_minmax(min,max);
    ffmin.set_dice(_chifn->get_dice());
    ffmin.use_gradient();
    ffmin.find_minimum(seed,trial);

    printf("    after dalex_simplex chimin %e\n",chimin());
    _simplex_mindex=mindex();

    array_1d<double> min_pt;
    ffmin.get_minpt(min_pt);
    double mu;
    int i_found;
    evaluate(min_pt, &mu, &i_found);

    _update_good_points(pt_0);

}


int dalex::bisection(int ilow, const array_1d<double> &dir, double local_target, double tol){
    safety_check("bisection(int, arr)");
    array_1d<double> trial_high;
    trial_high.set_name("dalex_bisection_trial_high");
    double rr=1.0;
    int ii,i_found;
    double mu=-2.0*exception_value;
    for(ii=0;ii<_chifn->get_dim();ii++){
        trial_high.set(ii,_chifn->get_pt(ilow,ii));
    }

    while(mu<=local_target){
        for(ii=0;ii<_chifn->get_dim();ii++){
            trial_high.add_val(ii,rr*dir.get_data(ii));
        }
        evaluate(trial_high, &mu, &i_found);
        rr*=2.0;
    }

    return bisection(_chifn->get_pt(ilow), trial_high, local_target, tol);
}


int dalex::bisection(int ilow, int ihigh, double local_target, double tol){
    safety_check("bisection(int, int)");
    return bisection(_chifn->get_pt(ilow), _chifn->get_pt(ihigh),
                     local_target, tol);
}


int dalex::bisection(const array_1d<double>& lowball_in, const array_1d<double>& highball_in,
                     double local_target, double tol){

    safety_check("bisection(arr, arr)");
    array_1d<double> trial,lowball,highball;
    trial.set_name("dalex_bisection_trial");
    lowball.set_name("dalex_bisection_lowball");
    highball.set_name("dalex_bisection_highball");
    int i_found,ii,jj,i_trial;

    double mu,flow,fhigh;
    evaluate(lowball_in, &flow, &ii);
    evaluate(highball_in, &fhigh, &jj);

    if(flow>local_target){
        printf("WARNING flow is greater than target %e %e\n",flow,local_target);
        printf("min %e target %e\n",chimin(),target());
        exit(1);
    }

    i_found=ii;

    for(ii=0;ii<_chifn->get_dim();ii++){
        lowball.set(ii,lowball_in.get_data(ii));
        highball.set(ii,highball_in.get_data(ii));
    }

    int ct;
    double wgt_low;
    for(ct=0;ct<20 &&
       (ct<5 || fabs(_chifn->get_fn(i_found)-local_target)>tol); ct++){

        wgt_low=(fhigh-local_target)/(fhigh-flow);
        if(wgt_low<0.25 || wgt_low>0.75){
            wgt_low=0.5;
        }

        for(ii=0;ii<_chifn->get_dim();ii++){
            trial.set(ii, wgt_low*lowball.get_data(ii)+(1.0-wgt_low)*highball.get_data(ii));
        }

        evaluate(trial,&mu,&i_trial);

        if(mu<local_target){
            flow=mu;
            for(ii=0;ii<_chifn->get_dim();ii++){
                lowball.set(ii,trial.get_data(ii));
            }
        }
        else{
            fhigh=mu;
            for(ii=0;ii<_chifn->get_dim();ii++){
                highball.set(ii,trial.get_data(ii));
            }
        }

        if(fabs(mu-local_target)<fabs(_chifn->get_fn(i_found)-local_target)){
            i_found=i_trial;
        }
    }

    return i_found;

}


void dalex::find_trial_bases(int idim, array_1d<double> &dx, array_2d<double> &bases_out){

    int i,j;

    array_1d<double> trial_dir;
    trial_dir.set_name("dalex_find_trial_bases_trial_dir");

    if(_basis_vectors.get_rows()==0){
        _basis_vectors.reset();
        _basis_vectors.set_cols(_chifn->get_dim());
        for(i=0;i<_chifn->get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                _basis_vectors.set(i,j,normal_deviate(_chifn->get_dice(),0.0,1.0));
            }
        }
    }

    if(dx.get_dim()==0){
        for(i=0;i<_chifn->get_dim();i++){
            dx.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dx.normalize();
    }

    bases_out.set_cols(_basis_vectors.get_cols());
    for(i=0;i<_chifn->get_dim();i++){
        for(j=0;j<_chifn->get_dim();j++){
            bases_out.set(i,j,_basis_vectors.get_data(i,j));
        }
    }

    for(i=0;i<_chifn->get_dim();i++){
        bases_out.add_val(idim,i,dx.get_data(i));
    }
    bases_out.normalize_row(idim);

    int ix,jx;
    double mu;

    for(ix=idim+1;ix!=idim;){
        if(ix>=_chifn->get_dim()){
            ix=0;
        }

        for(jx=idim;jx!=ix;){
            if(ix>=_chifn->get_dim()){
                jx=0;
            }

            mu=0.0;
            for(i=0;i<_chifn->get_dim();i++){
                mu+=bases_out.get_data(ix,i)*bases_out.get_data(jx,i);
            }
            for(i=0;i<_chifn->get_dim();i++){
                bases_out.subtract_val(ix,i,mu*bases_out.get_data(jx,i));
            }

            if(jx<_chifn->get_dim()-1)jx++;
            else jx=0;
        }

        bases_out.normalize_row(ix);

        if(ix<_chifn->get_dim()-1)ix++;
        else ix=0;

    }

    validate_bases(bases_out,"dalex_find_trial_bases");


}

void dalex::validate_bases(array_2d<double> &bases, char *whereami){
    int ix,i,jx;
    double mu;
    /////////////////testing
    for(ix=0;ix<_chifn->get_dim();ix++){
        bases.normalize_row(ix);
        mu=0.0;
        for(i=0;i<_chifn->get_dim();i++){
            mu+=bases.get_data(ix,i)*bases.get_data(ix,i);
        }
        if(fabs(mu-1.0)>1.0e-6){
            printf("WARNING in %s, square norm %e\n",whereami,mu);
            exit(1);
        }

        for(jx=ix+1;jx<_chifn->get_dim();jx++){
            mu=0.0;
            for(i=0;i<_chifn->get_dim();i++){
                mu+=bases.get_data(ix,i)*bases.get_data(jx,i);
            }

            if(fabs(mu)>1.0e-6){
                printf("WARNING in %s, dot product %e\n",whereami,mu);
                exit(1);
            }
        }
    }

}

void dalex::guess_bases(array_2d<double> &bases){
    safety_check("guess_bases");

    array_2d<double> covar;
    covar.set_name("dalex_guess_bases_covar");
    covar.set_cols(_chifn->get_dim());
    bases.set_cols(_chifn->get_dim());
    int ix,iy;
    double covarmax=-1.0;

    find_covariance_matrix(mindex(),covar);

    for(ix=0;ix<_chifn->get_dim();ix++){
        for(iy=ix;iy<_chifn->get_dim();iy++){
            if(fabs(covar.get_data(ix,iy))>covarmax){
                covarmax=fabs(covar.get_data(ix,iy));
            }
        }
    }

    for(ix=0;ix<_chifn->get_dim();ix++){
        for(iy=0;iy<_chifn->get_dim();iy++){
            covar.divide_val(ix,iy,covarmax);
        }
    }

    printf("assembled covariance matrix\n");

    array_2d<double> evecs;
    evecs.set_name("dalex_guess_bases_evecs");
    array_1d<double> evals;
    evals.set_name("dalex_guess_bases_evals");

    evecs.set_cols(_chifn->get_dim());

    try{
        eval_symm(covar,evecs,evals,0.1);
    }
    catch(int iex){
        printf("Guess failed on of eigen vectors\n");
        throw -1;
    }

    for(ix=0;ix<_chifn->get_dim();ix++){
        for(iy=0;iy<_chifn->get_dim();iy++){
            bases.set(ix,iy,evecs.get_data(ix,iy));
        }
        bases.normalize_row(ix);
    }

    validate_bases(bases,"dalex_guess_bases");

    printf("validated guessed bases\n");
}



void dalex::find_bases(){

    safety_check("find_bases");

    _minimizers.reset_preserving_room();

    array_1d<double> dir;
    int i,j;
    for(i=0;i<_basis_associates.get_dim();i++){
        if(_chifn->get_fn(_basis_associates.get_data(i))>target()){
            _basis_associates.remove(i);
            i--;
        }
    }
    int i_pt;
    double mu;
    int n_0=_basis_associates.get_dim();
    int old_search_type = _chifn->get_search_type();
    _chifn->set_search_type(_type_find_bases);
    while(_basis_associates.get_dim()<4*_chifn->get_dim()+n_0){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dir.normalize();
        i_pt=bisection(mindex(),dir,0.5*(target()+chimin()),0.1);

        /*if(fabs(_chifn->get_fn(i_pt)-0.5*(target()+chimin()))>1.0 && chimin()<500.0){
            printf("WARNING failed to get associate within tol %e %e-- %e %e\n",
            0.5*(target()+chimin()),_chifn->get_fn(i_pt),chimin(),target());
            printf("first\n");
            exit(1);
        }*/

        if(i_pt!=mindex() && _basis_associates.contains(i_pt)==0){
            _basis_associates.add(i_pt);
            i_pt=bisection(mindex(),dir,0.25*chimin()+0.75*target(),0.1);
            if(i_pt!=mindex() && _basis_associates.contains(i_pt)==0){
                _basis_associates.add(i_pt);

                /*if(fabs(_chifn->get_fn(i_pt)-(0.75*target()+0.25*chimin()))>1.0 && chimin()<500.0){
                    printf("WARNING failed to get associate within tol %e %e -- %e %e\n",
                    0.75*target()+0.25*chimin(),_chifn->get_fn(i_pt),chimin(),target());
                    printf("second\n");
                    exit(1);
                }*/

            }
        }
    }
    _chifn->set_search_type(old_search_type);

    array_2d<double> trial_bases;
    array_1d<double> trial_model,dx;

    trial_bases.set_name("dalex_find_bases_trial_bases");
    trial_model.set_name("dalex_find_bases_trial_model");
    dx.set_name("dalex_find_bases_dx");

    if(_basis_vectors.get_rows()==0){
        find_trial_bases(0,dx,trial_bases);
        for(i=0;i<_chifn->get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                _basis_vectors.set(i,j,trial_bases.get_data(i,j));
            }
        }
    }

    int dimsq=_chifn->get_dim()*_chifn->get_dim();

    if(_basis_associates.get_dim()==0){
        printf("WARNING _basis associates is empty\n");
        exit(1);
    }


    int ct,idim,aborted,max_abort,total_aborted,changed_bases;
    double error0,error,errorBest,stdev,stdevlim,error1;

    stdev=0.1/sqrt(double(_chifn->get_dim()));
    stdevlim=1.0e-5/sqrt(double(_chifn->get_dim()));
    max_abort=_chifn->get_dim()*100;

    error=basis_error(_basis_vectors,_basis_model);
    error0=error;
    error1=error;
    errorBest=error;
    aborted=0;
    total_aborted=0;
    changed_bases=0;
    ct=0;

    printf("center %e min %e\n",_chifn->get_fn(mindex()),chimin());
    printf("error0 %e\n", error0);

    if(1==1){
        printf("guessing basis from second derivative\n");
        try{
            guess_bases(trial_bases);
            error=basis_error(trial_bases,trial_model);
            printf("guess got error %e\n",error);
            if(error<error0){
                changed_bases=1;
                for(i=0;i<_chifn->get_dim();i++){
                    _basis_model.set(i,trial_model.get_data(i));
                    for(j=0;j<_chifn->get_dim();j++){
                        _basis_vectors.set(i,j,trial_bases.get_data(i,j));
                    }
                }
                error1=error;
                errorBest=error;
            }
        }
        catch(int iex){
            printf("never mind; guess did not work\n");
        }
    }

    double penalty=0.0;
    for(i=0;i<_chifn->get_dim();i++){
        if(_basis_model.get_data(i)<0.0){
            penalty+=1.0;
        }
    }

    double pp0=penalty;

    while(stdev>stdevlim && aborted<max_abort){
        ct++;
        idim=-1;
        while(idim>=_chifn->get_dim() || idim<0){
            idim=_chifn->random_int()%_chifn->get_dim();
        }

        for(i=0;i<_chifn->get_dim();i++){
            dx.set(i,normal_deviate(_chifn->get_dice(),0.0,stdev));
        }

        find_trial_bases(idim,dx,trial_bases);
        error=basis_error(trial_bases,trial_model);

        if(error<errorBest){
            penalty=0.0;
            for(i=0;i<_chifn->get_dim();i++){
                if(trial_model.get_data(i)<0.0){
                    penalty+=1.0;
                }
            }
            if(error1-error>1.0e-5*error){
                aborted=0;
            }
            else{
                aborted++;
                total_aborted++;
            }

            changed_bases=1;
            for(i=0;i<_chifn->get_dim();i++){
                _basis_model.set(i,trial_model.get_data(i));
                for(j=0;j<_chifn->get_dim();j++){
                    _basis_vectors.set(i,j,trial_bases.get_data(i,j));
                }
            }
            errorBest=error;
        }
        else{
            aborted++;
            total_aborted++;
        }

        if(ct%(max_abort/2)==0){
            if(total_aborted<(3*ct)/4)stdev*=1.5;
            else if(total_aborted>(3*ct)/4)stdev*=0.5;
        }

        if(ct%1000==0){
            error1=errorBest;
            printf("    ct %d error %.2e pp %.1e from %.2e %.1e min %.3e pts %d\n",
            ct,errorBest,penalty,error0,pp0,chimin(),_basis_associates.get_dim());
        }
    }

    int k;
    array_1d<double> min,max,focused_min,focused_max;
    min.set_name("find_bases_min");
    max.set_name("find_bases_max");
    focused_min.set_name("find_bases_focused_min");
    focused_max.set_name("find_bases_focused_max");
    if(changed_bases==1){
        for(i=0;i<_good_points.get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                mu=0.0;
                for(k=0;k<_chifn->get_dim();k++){
                    mu+=_chifn->get_pt(_good_points.get_data(i),k)*_basis_vectors.get_data(j,k);
                }
                if(i==0 || mu<min.get_data(j)){
                    min.set(j,mu);
                }
                if(i==0 || mu>max.get_data(j)){
                    max.set(j,mu);
                }

                if(_basis_associates.contains(_good_points.get_data(i))==1){
                    if(j>=focused_min.get_dim() || mu<focused_min.get_data(j)){
                        focused_min.set(j,mu);
                    }
                    if(j>=focused_max.get_dim() || mu>focused_max.get_data(j)){
                        focused_max.set(j,mu);
                    }
                }
            }
        }
        for(i=0;i<_chifn->get_dim();i++){
            _basis_norm.set(i,max.get_data(i)-min.get_data(i));
            if(_basis_norm.get_data(i)<1.0e-20){
                _basis_norm.set(i,1.0);
            }

            _basis_associate_norm.set(i,focused_max.get_data(i)-focused_min.get_data(i));
            if(_basis_associate_norm.get_data(i)<1.0e-20){
                printf("WARNING basis associate norm %d %e\n",
                i,_basis_associate_norm.get_data(i));
            }
        }

    }

    _basis_chimin=chimin();
    printf("done finding bases\n");


}


double dalex::basis_error(array_2d<double> &trial_bases, array_1d<double> &trial_model){

    if(_basis_associates.get_dim()<=0){
        printf("WARNING cannot calculate basis error there are only %d associates\n",
        _basis_associates.get_dim());

        exit(1);
    }

    safety_check("basis_error");

    trial_model.zero();
    if(_basis_ddsq.get_rows()>_basis_associates.get_dim()){
        _basis_ddsq.reset();
    }

    if(_basis_ddsq.get_cols()!=_chifn->get_dim()){
        _basis_ddsq.set_cols(_chifn->get_dim());
    }

    if(_basis_mm.get_dim()!=_chifn->get_dim()*_chifn->get_dim()){
        _basis_mm.set_dim(_chifn->get_dim()*_chifn->get_dim());
    }

    if(_basis_bb.get_dim()!=_chifn->get_dim()){
        _basis_bb.set_dim(_chifn->get_dim());
    }

    if(_basis_vv.get_dim()!=_chifn->get_dim()){
        _basis_vv.set_dim(_chifn->get_dim());
    }

    _basis_mm.zero();
    _basis_bb.zero();
    _basis_vv.zero();
    _basis_ddsq.zero();

    int i,j,ix;
    double mu;
    for(ix=0;ix<_basis_associates.get_dim();ix++){
        for(i=0;i<_chifn->get_dim();i++){
            mu=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                mu+=(_chifn->get_pt(_basis_associates.get_data(ix),j)-_chifn->get_pt(mindex(),j))*trial_bases.get_data(i,j);
            }
            _basis_ddsq.set(ix,i,mu*mu);
        }
    }

    for(i=0;i<_chifn->get_dim();i++){
        for(j=0;j<_basis_associates.get_dim();j++){
            _basis_bb.add_val(i,_basis_ddsq.get_data(j,i)*(_chifn->get_fn(_basis_associates.get_data(j))-chimin()));
        }
    }

    int k;
    for(i=0;i<_chifn->get_dim();i++){
        for(j=i;j<_chifn->get_dim();j++){
            ix=i*_chifn->get_dim()+j;
            for(k=0;k<_basis_associates.get_dim();k++){
                _basis_mm.add_val(ix,_basis_ddsq.get_data(k,i)*_basis_ddsq.get_data(k,j));
            }
            if(j!=i){
                _basis_mm.set(j*_chifn->get_dim()+i,_basis_mm.get_data(ix));
            }
        }
    }

    try{
        naive_gaussian_solver(_basis_mm,_basis_bb,trial_model,_chifn->get_dim());
    }
    catch(int iex){
        printf("WARNING basis_error was no good\n");
        return 2.0*exception_value;
    }

    double error=0.0,chi_model;
    for(i=0;i<_basis_associates.get_dim();i++){
        chi_model=chimin();
        for(j=0;j<_chifn->get_dim();j++){
            chi_model+=trial_model.get_data(j)*_basis_ddsq.get_data(i,j);
        }
        error+=power(_chifn->get_fn(_basis_associates.get_data(i))-chi_model,2);
    }

    return error/double(_basis_associates.get_dim());

}


void dalex::find_covariance_matrix(int iCenter, array_2d<double> &covar){
    safety_check("find_covariance_matrix");


    array_1d<double> trial,center,norm;
    trial.set_name("dalex_findCovar_trial");
    norm.set_name("dalex_findCovar_norm");
    center.set_name("dalex_findCovar_center");

    double fcenter;
    int i;

    fcenter=_chifn->get_fn(iCenter);

    for(i=0;i<_chifn->get_dim();i++){
        center.set(i,_chifn->get_pt(iCenter,i));
        norm.set(i,1.0);
    }

    array_2d<double> fpp,fpm,fmp,fmm;
    fpp.set_name("dalex_findCovar_fpp");
    fpm.set_name("dalex_findCovar_fpm");
    fmp.set_name("dalex_findCovar_fmp");
    fmm.set_name("dalex_findCovar_fmm");

    fpp.set_cols(_chifn->get_dim());
    fpm.set_cols(_chifn->get_dim());
    fmp.set_cols(_chifn->get_dim());
    fmm.set_cols(_chifn->get_dim());

    array_1d<double> f2p,f2m;
    f2p.set_name("dalex_findCovar_f2p");
    f2m.set_name("dalex_findCovar_f2m");

    int ifpp,ifpm,ifmp,ifmm,if2p,if2m;

    double mu;
    array_1d<double> dx;
    dx.set_name("dalex_findCovar_dx");
    for(i=0;i<_chifn->get_dim();i++){
        dx.set(i,1.0e-2);
    }

    int ix,iy,keepGoing,ctAbort,ctAbortMax,calledMax;
    int ibefore=_chifn->get_called();

    ctAbort=0;
    ctAbortMax=100;
    calledMax = 10*_chifn->get_dim()*_chifn->get_dim();

    for(i=0;i<_chifn->get_dim();i++){
        trial.set(i,center.get_data(i));
    }

    for(ix=0;ix<_chifn->get_dim();ix++){
        for(i=0;i<_chifn->get_dim();i++){
            trial.set(i,center.get_data(i));
        }

        if(_chifn->get_called()-ibefore>calledMax ||
           ctAbort>=ctAbortMax){
                printf("Could not find CoVar; aborting\n");
                printf("ctAbort %d\n",ctAbort);
                printf("called %d\n",_chifn->get_called()-ibefore);
                throw -1;
        }

        if(iCenter<0){
            printf("Center is invalid; aborting\n");
            throw -1;
        }

        keepGoing=1;
        trial.set(ix,center.get_data(ix)+2.0*dx.get_data(ix)*norm.get_data(ix));
        evaluate(trial,&mu,&if2p);

        if(if2p>=0 && if2p!=iCenter){
            f2p.set(ix,mu);
        }
        else if(if2p==iCenter){
            dx.multiply_val(ix,1.5);
            keepGoing=0;
            ix--;
            ctAbort++;
        }
        else if(if2p<0){
            center.subtract_val(ix,2.5*dx.get_data(ix)*norm.get_data(ix));
            evaluate(center,&fcenter,&iCenter);
            keepGoing=0;
            ix--;
            ctAbort++;
        }

        if(keepGoing==1){
            trial.set(ix,center.get_data(ix)-2.0*dx.get_data(ix)*norm.get_data(ix));
            evaluate(trial,&mu,&if2m);

            if(if2m>=0 && if2m!=iCenter){
                f2m.set(ix,mu);
            }
            else if(if2m==iCenter){
                dx.multiply_val(ix,1.5);
                keepGoing=0;
                ix--;
                ctAbort++;
            }
            else if(if2m<0){
                center.add_val(ix,2.5*dx.get_data(ix)*norm.get_data(ix));
                evaluate(center,&fcenter,&iCenter);
                keepGoing=0;
                ix--;
                ctAbort++;
            }
        }

        for(iy=ix-1;iy>=0 && keepGoing==1;iy--){
            for(i=0;i<_chifn->get_dim();i++){
                trial.set(i,center.get_data(i));
            }

            if(_chifn->get_called()-ibefore>calledMax ||
               ctAbort>=ctAbortMax){
                printf("Could not find CoVar; aborting\n");
                printf("ctAbort %d\n",ctAbort);
                printf("called %d\n",_chifn->get_called()-ibefore);
                throw -1;
            }

            if(iCenter<0){
                printf("center is invalid; aborting\n");
                throw -1;
            }

            trial.set(ix,center.get_data(ix)+dx.get_data(ix)*norm.get_data(ix));
            trial.set(iy,center.get_data(iy)+dx.get_data(iy)*norm.get_data(iy));
            evaluate(trial,&mu,&ifpp);
            if(ifpp>=0 && ifpp!=iCenter){
                fpp.set(ix,iy,mu);
            }
            else if(ifpp==iCenter){
                dx.multiply_val(ix,1.5);
                dx.multiply_val(iy,1.5);
                keepGoing=0;
                ix--;
                ctAbort++;
            }
            else if(ifpp<0){
                center.subtract_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                center.subtract_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                evaluate(center,&fcenter,&iCenter);
                keepGoing=0;
                ix--;
                ctAbort++;
            }

            if(keepGoing==1){
               trial.set(iy,center.get_data(iy)-dx.get_data(iy)*norm.get_data(iy));
               evaluate(trial,&mu,&ifpm);
               if(ifpm>=0 && ifpm!=iCenter){
                   fpm.set(ix,iy,mu);
               }
               else if(ifpm==iCenter){
                   dx.multiply_val(ix,1.5);
                   dx.multiply_val(iy,1.5);
                   keepGoing=0;
                   ix--;
                   ctAbort++;
               }
               else if(ifpm<0){
                   center.subtract_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                   center.add_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                   evaluate(center,&fcenter,&iCenter);
                   keepGoing=0;
                   ix--;
                   ctAbort++;
               }
            }

            if(keepGoing==1){
                trial.set(ix,center.get_data(ix)-dx.get_data(ix)*norm.get_data(ix));
                evaluate(trial,&mu,&ifmm);
                if(ifmm>=0 && ifmm!=iCenter){
                    fmm.set(ix,iy,mu);
                }
                else if(ifmm==iCenter){
                    dx.multiply_val(ix,1.5);
                    dx.multiply_val(iy,1.5);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                }
                else if(ifmm<0){
                    center.add_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                    center.add_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                    evaluate(center,&fcenter,&iCenter);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                }
            }

            if(keepGoing==1){
                trial.set(iy,center.get_data(iy)+dx.get_data(iy)*norm.get_data(iy));
                evaluate(trial,&mu,&ifmp);
                if(ifmp>=0 && ifmp!=iCenter){
                    fmp.set(ix,iy,mu);
                }
                else if(ifmp==iCenter){
                    dx.multiply_val(ix,1.5);
                    dx.multiply_val(iy,1.5);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                }
                else if(ifmp<0){
                    center.add_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                    center.subtract_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                    evaluate(center,&fcenter,&iCenter);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                }
            }

        }
    }

    covar.set_cols(_chifn->get_dim());
    for(ix=0;ix<_chifn->get_dim();ix++){
        covar.set(ix,ix,0.25*(f2p.get_data(ix)+f2m.get_data(ix)-2.0*fcenter)/power(dx.get_data(ix)*norm.get_data(ix),2));
    }

    double num,denom;
    for(ix=0;ix<_chifn->get_dim();ix++){
        for(iy=ix-1;iy>=0;iy--){
            num=0.25*(fpp.get_data(ix,iy)+fmm.get_data(ix,iy)-fmp.get_data(ix,iy)-fpm.get_data(ix,iy));
            denom=dx.get_data(ix)*norm.get_data(ix)*dx.get_data(iy)*norm.get_data(iy);
            covar.set(ix,iy,num/denom);
            covar.set(iy,ix,num/denom);
        }
    }


}

void dalex::get_negative_gradient(int i_origin, cost_fn &cost, ellipse &dummy_ellipse, array_1d<double> &out_dir){
    double mu_pos, mu_neg;
    array_1d<double> trial;
    trial.set_name("get_neg_grad_trial");
    int i,j;
    out_dir.reset_preserving_room();
    for(i=0;i<_chifn->get_dim();i++){
        out_dir.set(i,0.0);
    }
    double rat=0.01;
    double delta;
    for(i=0;i<_chifn->get_dim();i++){
        delta=0.0;
        for(j=0;j<_chifn->get_dim();j++){
            delta+=(_chifn->get_pt(i_origin,j)-_chifn->get_pt(mindex(),j))*dummy_ellipse.bases(i,j);
        }
        delta*=rat;

        for(j=0;j<_chifn->get_dim();j++){
            trial.set(j,_chifn->get_pt(i_origin,j)+rat*dummy_ellipse.radii(i)*dummy_ellipse.bases(i,j));
        }
        mu_pos=cost(trial);
        for(j=0;j<_chifn->get_dim();j++){
            trial.set(j,_chifn->get_pt(i_origin,j)-rat*dummy_ellipse.radii(i)*dummy_ellipse.bases(i,j));
        }
        mu_neg=cost(trial);
        for(j=0;j<_chifn->get_dim();j++){
            out_dir.add_val(j,(mu_neg-mu_pos)*dummy_ellipse.bases(i,j)/(2.0*rat*dummy_ellipse.radii(i)));
        }
    }
}

int dalex::simplex_boundary_search(const int specified, const int i_origin,
                                   ellipse_list &exclusion_zones, int *i_next){

    safety_check("simplex_boundary_search");
    printf("\ndoing dalex.simplex_boundary_search() %d\n",_chifn->get_pts());
    int pt_start=_chifn->get_pts();
    assess_good_points();

    if(_chifn->get_dim()>9 && specified>=0){
        printf("    starting from %e -- %e %e\n",
        _chifn->get_fn(specified),_chifn->get_pt(specified,6), _chifn->get_pt(specified,9));
    }

    int i_node,i_pt;
    int i,j,k;
    double xmin,xmax,xx;

    int n_good_0=_good_points.get_dim();

    array_1d<int> associates;
    associates.set_name("dalex_simplex_boundary_associates");
    int i_start;
    int n_thin=-1;
    int ip,io;

    double mu;
    array_1d<double> trial;
    trial.set_name("simplex_boundary_trial");

    if(_tendril_path.get_rows()>20000){
        n_thin=_tendril_path.get_rows()/20000;
        printf("    thinning %d\n",n_thin);
    }

    if(n_thin==1){
       n_thin=2;
    }

    for(i=0;i<_tendril_path.get_rows();i++){
        ip=_tendril_path.get_data(i,0);
        io=_tendril_path.get_data(i,1);
        if(n_thin<0 || i%n_thin==0){
            if(specified<=0 || _has_struck==0){
                associates.add(ip);
            }
            else{
                for(j=0;j<_chifn->get_dim();j++){
                    trial.set(j,0.5*(_chifn->get_pt(specified,j)+_chifn->get_pt(io,j)));
                }
                evaluate(trial,&mu,&k);
                if(mu<_chifn->target()){
                    associates.add(ip);
                }
            }
        }
    }

    if(associates.get_dim()==0){
        for(i=0;i<_good_points.get_dim();i++){
            associates.add(_good_points.get_data(i));
        }
    }

    cost_fn dchifn(_chifn,associates);

    printf("    associates %d path %d\n", associates.get_dim(),_tendril_path.get_rows());

    simplex_minimizer ffmin;
    ffmin.set_chisquared(&dchifn);
    ffmin.set_dice(_chifn->get_dice());
    array_1d<double> min,max;
    min.set_name("dalex_simplex_search_min");
    max.set_name("dalex_simplex_search_min");

    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,0.0);
        max.set(i,_chifn->get_characteristic_length(i));
    }

    ffmin.set_minmax(min,max);
    ffmin.use_gradient();

    int old_type=_chifn->get_search_type();
    _chifn->set_search_type(_type_tendril_seed);

    array_2d<double> seed;
    seed.set_name("dalex_simplex_search_seed");
    array_1d<double> trial1,trial2;
    trial1.set_name("dalex_simplex_trial1");
    trial2.set_name("dalex_simplex_trial2");

    int d_step = _strikes+1;

    array_2d<double> ellipse_pts;
    ellipse_pts.set_name("dalex_simplex_boundary_ellipse_pts");
    ellipse dummy_ellipse;
    for(i=i_origin;i<specified;i++){
        if(_chifn->get_fn(i)<target()){
            ellipse_pts.add_row(_chifn->get_pt(i));
        }
    }

    for(i=i_origin-1;i>0 && ellipse_pts.get_rows()<2*_chifn->get_dim();i--){
        if(_chifn->get_fn(i)<target()){
            ellipse_pts.add_row(_chifn->get_pt(i));
        }
    }

    dummy_ellipse.build(ellipse_pts);

    int i_anchor;
    array_1d<double> base_dir;
    base_dir.set_name("dalex_simplex_boundary_base_dir");
    for(i=0;i<_chifn->get_dim();i++){
        base_dir.set(i,_chifn->get_pt(specified,i)-_chifn->get_pt(i_origin,i));
    }
    base_dir.normalize();

    i_anchor=specified;

    array_1d<double> anchor;
    anchor.set_name("simplex_boundary_anchor");
    for(i=0;i<_chifn->get_dim();i++){
        anchor.set(i,_chifn->get_pt(i_anchor,i));
    }
    array_1d<double> bisect_dir,epsilon;
    bisect_dir.set_name("simplex_boundary_bisect_dir");
    epsilon.set_name("simplex_boundary_bisect_dir");
    int i_bisect1,i_bisect2,i_chosen;
    double mu1,mu2;
    double component,rat;
    array_1d<int> kept_dex;
    kept_dex.set_name("simplex_boundary_kept_dex");
    array_1d<double> avg_pt;
    avg_pt.set_name("simplex_boundary_avg_pt");

    double local_target;

    if(specified>=0){
        i_anchor=specified;
        if(_chifn->get_fn(i_anchor)<target()){
            local_target=target();
        }
        else{
            local_target=_chifn->get_fn(i_anchor)+0.5*(target()-chimin());
        }
        while(seed.get_rows()!=_chifn->get_dim()+1){
            for(i=0;i<_chifn->get_dim();i++){
                for(j=0;j<_chifn->get_dim();j++){
                    bisect_dir.set(j,dummy_ellipse.bases(i,j)+base_dir.get_data(j));
                }
                i_bisect1=bisection(i_anchor,bisect_dir,local_target,0.001);
                for(j=0;j<_chifn->get_dim();j++){
                    bisect_dir.set(j,-1.0*dummy_ellipse.bases(i,j)+base_dir.get_data(j));
                }
                i_bisect2=bisection(i_anchor,bisect_dir,local_target,0.001);
                mu1=dchifn(_chifn->get_pt(i_bisect1));
                mu2=dchifn(_chifn->get_pt(i_bisect2));

                i_chosen=-1;
                if(i_bisect1==i_anchor && i_bisect2!=i_anchor){
                    i_chosen=i_bisect2;
                }
                else if(i_bisect1!=i_anchor && i_bisect2==i_anchor){
                    i_chosen=i_bisect1;
                }
                else if(i_bisect1!=i_anchor && i_bisect2!=i_anchor){
                    if(mu1<mu2){
                        i_chosen=i_bisect1;
                    }
                    else{
                        i_chosen=i_bisect2;
                    }
                }
                else{
                    i_anchor--;
                    while(_chifn->get_fn(i_anchor)>target() || i_anchor==i_origin){
                        i_anchor--;
                    }
                    seed.reset_preserving_room();
                    kept_dex.reset_preserving_room();
                    break;
                }

                if(i_chosen>=0){
                    for(j=0;j<_chifn->get_dim();j++){
                        bisect_dir.set(j,0.5*(_chifn->get_pt(i_anchor,j)+_chifn->get_pt(i_chosen,j)));
                    }
                    evaluate(bisect_dir,&mu1,&j);
                    seed.add_row(bisect_dir);
                    kept_dex.add(i_chosen);
                }
            }

            if(seed.get_rows()==_chifn->get_dim()){
                for(i=0;i<_chifn->get_dim();i++){
                    avg_pt.set(i,0.0);
                }
                for(i=0;i<kept_dex.get_dim();i++){
                   for(j=0;j<_chifn->get_dim();j++){
                       avg_pt.add_val(j,_chifn->get_pt(kept_dex.get_data(i),j));
                   }
                }
                for(i=0;i<_chifn->get_dim();i++){
                    avg_pt.divide_val(i,double(kept_dex.get_dim()));
                }
                evaluate(avg_pt,&mu1,&i);
                seed.add_row(avg_pt);
            }
        }
    }
    else{
        printf("calling _explorers.get_seed(); did not expect that\n");
        exit(1);
        _explorers.get_seed(seed);
    }

    printf("fn anchor %e; %d %d %d\n",
    _chifn->get_fn(i_anchor),specified,i_origin,i_anchor);

    _chifn->set_search_type(old_type);

    double mu_min;
    double start_min,start_max;
    int i_start_min;
    array_1d<double> start_vals,start_vals_sorted;
    array_1d<int> start_vals_dex;
    start_vals.set_name("simplex_boundary_start_vals");
    start_vals_sorted.set_name("simplex_boundary_start_vals_sorted");
    start_vals_dex.set_name("simplex_boundary_start_vals_dex");
    for(i=0;i<seed.get_rows();i++){
        mu=dchifn(seed(i));
        start_vals.add(mu);
        start_vals_dex.add(i);
        if(i==0 || mu<start_min){
            start_min=mu;
            _chifn->evaluate(seed(i),&mu,&i_start_min);
        }
        if(i==0 || mu>start_max){
            start_max=mu;
        }
    }
    sort(start_vals,start_vals_sorted,start_vals_dex);
    mu=start_vals_sorted.get_data(start_vals_dex.get_dim()/2);
    printf("    starting from %e; %e; %e\n",start_min,mu,start_max);
    printf("    starting from (delta) %e; %e; %e\n",start_min-target(),mu-target(),start_max-target());

    array_1d<double> minpt;
    minpt.set_name("dalex_simplex_search_minpt");

    for(i=pt_start;i<_chifn->get_pts();i++){
        mu=dchifn(_chifn->get_pt(i));
    }

    ffmin.find_minimum(seed,minpt);

    array_1d<int> path_row;
    path_row.set_name("path_row");

    int start_path;
    if(pt_start<i_start_min){
        start_path=i_start_min;
    }
    else{
        start_path=pt_start;
    }

    i_next[0]=-1;
    double fn,cost_min;
    int valid_cache;
    cost_min=2.0*exception_value;
    for(i=pt_start;i<_chifn->get_pts();i++){
        if(_chifn->get_fn(i)<target()){
            valid_cache=dchifn.get_cached_values(i,&fn);
            if(valid_cache==0){
                printf("bad cache %d %d\n",i,pt_start);
                exit(1);
            }
            if(valid_cache==1 && fn<cost_min){
                i_next[0]=i;
                cost_min=fn;
            }
        }
    }
    if(i_next[0]<0){
        printf("WARNING could not set i_next\n");
        exit(1);
    }

    int pre_fill=_chifn->get_pts();

    // *** do bisection along directions perpendicular to ***
    // *** i_next - specified ***
    base_dir.reset_preserving_room();
    array_1d<double> trial_dir;
    base_dir.set_name("dalex_boundary_base_dir");
    trial_dir.set_name("dalex_boundary_trial_dir");
    array_2d<double> perp_dir;
    perp_dir.set_name("dalex_boundary_perp_dir");
    array_1d<int> good_dexes;
    good_dexes.set_name("dalex_boundary_good_dexes");
    int i_midst;
    int pre_midpt=_chifn->get_pts();
    double sgn;
    if(i_next[0]!=specified){
        for(i=0;i<_chifn->get_dim();i++){
            base_dir.set(i,_chifn->get_pt(i_next[0],i)-_chifn->get_pt(specified,i));
        }
        base_dir.normalize();
        while(perp_dir.get_rows()<_chifn->get_dim()-1){
            for(i=0;i<_chifn->get_dim();i++){
                trial_dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
            }
            component=0.0;
            for(i=0;i<_chifn->get_dim();i++){
                component+=trial_dir.get_data(i)*base_dir.get_data(i);
            }
            for(i=0;i<_chifn->get_dim();i++){
                trial_dir.subtract_val(i,component*base_dir.get_data(i));
            }
            for(i=0;i<perp_dir.get_rows();i++){
                component=0.0;
                for(j=0;j<_chifn->get_dim();j++){
                    component+=trial_dir.get_data(j)*perp_dir.get_data(i,j);
                }
                for(j=0;j<_chifn->get_dim();j++){
                    trial_dir.subtract_val(j,component*perp_dir.get_data(i,j));
                }
            }
            component=trial_dir.normalize();
            if(component>1.0e-20){
                perp_dir.add_row(trial_dir);
            }
        }

        for(i=specified;i<i_next[0];i++){
            if(_chifn->get_fn(i)<target()){
                good_dexes.add(i);
            }
        }
        i_midst=good_dexes.get_data(good_dexes.get_dim()/2);
        for(i=0;i<perp_dir.get_rows();i++){
            for(sgn=-1.0;sgn<1.1;sgn+=2.0){
                for(j=0;j<_chifn->get_dim();j++){
                    trial_dir.set(j,sgn*perp_dir.get_data(i,j));
                }
                bisection(i_midst,trial_dir,target(),0.01);
            }
        }
        j=0;
        for(i=pre_midpt;i<_chifn->get_pts();i++){
            if(_chifn->get_fn(i)<target()){
                j++;
            }
        }
        printf("bisection sampled %d good %d\n",
               _chifn->get_pts()-pre_midpt,j);
    }

    // *** try to fill in the local ellipse ***
    ellipse_pts.reset_preserving_room();
    for(i=specified;i<_chifn->get_pts();i++){
        if(_chifn->get_fn(i)<target()){
            ellipse_pts.add_row(_chifn->get_pt(i));
        }
    }

    dummy_ellipse.build(ellipse_pts);
    if(_ellipse_sampler.is_initialized()==0){
        _ellipse_sampler.initialize(_chifn->get_dim(), _chifn->random_int()%1000000+1);
    }

    array_1d<double> ell_pt,center;
    ell_pt.set_name("dalex_simplex_boundary_ell_pt");
    center.set_name("dalex_simplex_boundary_center");
    trial.reset_preserving_room();
    int n_fill=(_chifn->get_pts()-specified)/10;
    int n_good=0;
    int n_bisect=0;
    array_1d<int> new_good;
    new_good.set_name("dalex_simplex_new_good");
    int i_fill_start=_chifn->get_pts();
    while(_chifn->get_pts()-i_fill_start<n_fill){
        _ellipse_sampler.get_pt(ell_pt);
        for(j=0;j<_chifn->get_dim();j++){
            trial.set(j,dummy_ellipse.center(j));
        }
        for(j=0;j<_chifn->get_dim();j++){
            for(k=0;k<_chifn->get_dim();k++){
                trial.add_val(k,ell_pt.get_data(j)*dummy_ellipse.radii(j)*dummy_ellipse.bases(j,k));
            }
        }
        evaluate(trial,&mu,&j);
        if(mu<target()){
            n_good++;
        }
        else{
            for(i=0;i<_chifn->get_dim();i++){
                center.set(i,dummy_ellipse.center(i));
            }
            i=bisection(center,trial,target(),0.001);
            if(new_good.contains(i)==0){
                n_bisect++;
                new_good.add(i);
            }
        }
    }
    printf("    n_fill %d n_good %d n_bisect %d -- %d\n",n_fill,n_good,n_bisect,new_good.get_dim());

    int i_good_start;

    _update_good_points(pt_start);

    for(i=pre_fill;i<_chifn->get_pts();i++){
        if(_chifn->get_fn(i)<target()){
            fn=dchifn(_chifn->get_pt(i));
            if(fn<cost_min){
                i_next[0]=i;
                cost_min=fn;
            }
        }
    }
    if(i_next[0]<0){
        printf("WARNING could not set i_next\n");
        exit(1);
    }

    if(_chifn->get_dim()>9){
        printf("    actually found %e -- %.3e %.3e; %.3e %.3e\n",
        _chifn->get_fn(i_next[0]),_chifn->get_pt(i_next[0],6), _chifn->get_pt(i_next[0],9),
        _chifn->get_pt(i_next[0],0), _chifn->get_pt(i_next[0], 1));
    }

    printf("    adjusted %e from %e\n",
    dchifn(_chifn->get_pt(i_next[0])),_chifn->get_fn(i_next[0]));

    printf("    min is %e target %e\n",chimin(),target());
    if(_chifn->get_dim()>9){
       printf("    minpt at %e %e\n",
       _chifn->get_pt(mindex(),6),
       _chifn->get_pt(mindex(),9));
    }


    for(i=start_path;i<_chifn->get_pts();i++){
        if(_tendril_path.get_rows()==0){
            _tendril_path.set_cols(2);
        }
        if(_chifn->get_fn(i)<_chifn->target() &&
           _chifn->get_search_type_log(i)==_type_tendril){
            path_row.set(0,i);
            if(distance(i,i_next[0])<distance(i,i_start_min)){
                path_row.set(1,i_next[0]);
            }
            else{
                path_row.set(1,i_start_min);
            }
            _tendril_path.add_row(path_row);
        }

    }

    int is_a_strike=0;

    for(i=0;i<exclusion_zones.ct() && is_a_strike==0;i++){
        if(exclusion_zones(i)->contains(_chifn->get_pt(i_next[0]))==1){
            is_a_strike=1;
        }
    }

    if(is_a_strike==1){
        return 1;
    }
    else{
        return 0;
    }

}

void dalex::explore(){
    explore(0);
}

void dalex::explore(int with_kick){
    printf("\nexploring\n");
    int pt_0=_chifn->get_pts();

    if(_explorers.get_n_particles()<2*_chifn->get_dim()){
        _explorers.set_n_particles(2*_chifn->get_dim());
    }
    array_1d<int> associates;
    associates.set_name("dalex_explore_associates");
    int i;
    for(i=0;i<_good_points.get_dim();i++){
        if(_chifn->get_search_type_log(_good_points.get_data(i))<_type_tendril){
            associates.add(_good_points.get_data(i));
        }
    }

    _explorers.set_associates(associates);
    _explorers.set_envelope(4.0);
    _explorers.sample(4*_chifn->get_dim(),with_kick);

    _update_good_points(pt_0);
}

void dalex::min_explore(int n_particles, int n_steps){
    printf("\nmin exploring\n");
    int pt_0=_chifn->get_pts();

    _min_explorers.set_n_particles(n_particles);
    array_1d<int> associates;
    associates.set_name("dalex_explore_associates");
    int skip;
    if(_good_points.get_dim()<5000){
        skip=-1;
    }
    else{
        skip=5;
    }
    int i;
    for(i=0;i<_good_points.get_dim();i++){
        if(skip<0 || i%skip==0){
            associates.add(_good_points.get_data(i));
        }
    }

    _min_explorers.set_associates(associates);
    _min_explorers.sample(n_steps);

    _update_good_points(pt_0);
}

void dalex::tendril_search(int specified){

    _has_struck=0;

    int i,j,k;
    int pt_0=_chifn->get_pts();
    assess_good_points();
    _update_good_points();
    int n_good_0=_good_points.get_dim();

    double mu;
    int i_found;

    int i_exclude;
    int i_particle;
    array_1d<int> end_pts;
    end_pts.set_name("dalex_tendril_end_pts");
    array_2d<double> exclusion_points;
    array_1d<int> exclusion_dex;
    ellipse local_ellipse;

    int i_origin=0;

    simplex_boundary_search(specified, i_origin, _exclusion_zones, &i_particle);
    i_origin=specified;
    end_pts.add(i_particle);
    for(i=pt_0;i<_chifn->get_pts();i++){
        if(_chifn->get_fn(i)<target()){
            exclusion_points.add_row(_chifn->get_pt(i));
            exclusion_dex.add(i);
        }
    }
    local_ellipse.build(exclusion_points);
    i_exclude=_chifn->get_pts();
    _update_good_points();

    assess_good_points();

    array_1d<double> dir1,dir2,trial_center;
    dir1.set_name("dalex_simplex_boundary_dir1");
    dir2.set_name("dalex_simplex_boundary_dir2");
    trial_center.set_name("dalex_simplex_boundary_trial_center");

    int fall_back;
    int fall_back_origin;
    int ct_last;

    _strikes=0;
    int iteration=0;
    int is_a_strike;
    int i_next;
    double volume,volume_0;

    volume=1.0;
    for(i=0;i<_chifn->get_dim();i++){
        volume*=local_ellipse.radii(i);
    }
    volume_0=volume;

    int in_old_ones;
    double old_volume;

    fall_back=i_particle;
    fall_back_origin=specified;

    while(_strikes<3 && (_limit<0 || _chifn->get_pts()<_limit)){

        iteration++;

        ct_last=_chifn->get_pts();
        in_old_ones=simplex_boundary_search(i_particle, i_origin, _exclusion_zones, &i_next);
        end_pts.add(i_next);

        is_a_strike=in_old_ones;

        for(i=i_exclude;i<_chifn->get_pts();i++){
            if(_chifn->get_fn(i)<target()){
                exclusion_points.add_row(_chifn->get_pt(i));
                exclusion_dex.add(i);
            }
        }
        i_exclude=_chifn->get_pts();

        if(local_ellipse.contains(_chifn->get_pt(i_next), 1)==1){
            is_a_strike=1;
        }

        local_ellipse.build(exclusion_points);

        printf("    exclusion zones %d\n",_exclusion_zones.ct());

        i_origin=i_particle;
        i_particle=i_next;

        old_volume=volume_0;
        volume=1.0;
        for(i=0;i<_chifn->get_dim();i++){
            volume*=local_ellipse.radii(i);
        }
        if(in_old_ones==0 && volume>volume_0){
            is_a_strike=0;
            volume_0=volume;
        }

        if(is_a_strike==0 && in_old_ones==1){
            printf("WARNING is_a_strike %d; in_old_ones %d; should not happen\n",
                   is_a_strike, in_old_ones);
            exit(1);
        }

        if(is_a_strike==1){
            _strikes++;
            _has_struck=1;
            if(_strikes<3){
                i_particle=fall_back;
                i_origin=fall_back_origin;
            }
        }
        else{
            _strikes=0;
            fall_back=i_particle;
            fall_back_origin=i_origin;
        }

        printf("    volume %e from %e-- %d; chifn(i_next) %e\n",
               volume,old_volume,_exclusion_zones.ct(),_chifn->get_fn(i_next));
        printf("    strikes %d has struck %d\n",_strikes,_has_struck);

    }

    for(i=i_exclude;i<_chifn->get_pts();i++){
            if(_chifn->get_fn(i)<target()){
                exclusion_points.add_row(_chifn->get_pt(i));
                exclusion_dex.add(i);
            }
    }

    local_ellipse.build(exclusion_points);

    _exclusion_zones.add(local_ellipse);
    printf("\n    strike out (%d strikes; %d pts)\n",
           _strikes,_chifn->get_pts());

    _chifn->write_pts();

}

void dalex::iterate_on_minimum(){
    printf("iterating with %e\n",chimin());
    double min_0=chimin();
    double min_00=min_0;

    double min_1=-2.0*exception_value;
    int i,j;

    if(_reset_chimin>exception_value){
        _reset_chimin=chimin();
    }

    if(_good_points.get_dim()==0){
        find_bases();
    }

    while(min_1<min_0){
        min_0=chimin();
        min_explore(2*_chifn->get_dim(), 4*_chifn->get_dim());
        simplex_search(mindex());
        min_1=chimin();
    }

    if(chimin()<_reset_chimin-_reset_threshold){
        _reset_chimin=chimin();
        _good_points.reset_preserving_room();
        _explorers.reset();
        _tendril_path.reset_preserving_room();
        _exclusion_zones.reset();
        find_bases();
    }
    else if(chimin()<min_00-0.01){
        find_bases();
    }

    printf("done iterating %e %d\n",chimin(),_chifn->get_called());

}

void dalex::refine_minimum(){

    printf("    refining minimum %e\n",chimin());

    int n_particles=2*_chifn->get_dim();
    int n_steps=5*_chifn->get_dim();
    int accepted=0;
    int rejected=0;

    array_1d<double> trial,coord;
    trial.set_name("refine_minimum_trial");
    coord.set_name("refine_minimum_coord");
    double roll,ratio,mu;
    int i_found;
    int i,j;

    _minimizers.reset_preserving_room();
    while(_minimizers.get_dim()!=n_particles){
        for(i=0;i<_chifn->get_dim();i++){
           coord.set(i,2.0*(_chifn->random_double()-0.5)*
                       _basis_associate_norm.get_data(i));
        }

        for(i=0;i<_chifn->get_dim();i++){
            trial.set(i,_chifn->get_pt(mindex(),i));
        }

        for(i=0;i<_chifn->get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                trial.add_val(j,coord.get_data(i)*_basis_vectors.get_data(i,j));
            }
        }

        evaluate(trial, &mu, &i_found);
        if(i_found>=0){
            _minimizers.add(i_found);
        }
    }

    int i_step,ip,i0;
    int i_dim,accept_it;
    double rr;
    for(i_step=0;i_step<n_steps;i_step++){
        for(ip=0;ip<_minimizers.get_dim();ip++){
            i0=_minimizers.get_data(ip);
            i_dim=_chifn->random_int()%_chifn->get_dim();
            rr=normal_deviate(_chifn->get_dice(),0.0,1.0)*0.5*_basis_associate_norm.get_data(i_dim);
            for(i=0;i<_chifn->get_dim();i++){
                trial.set(i,_chifn->get_pt(i0,i));
            }
            for(i=0;i<_chifn->get_dim();i++){
                trial.add_val(i,rr*_basis_vectors.get_data(i_dim,i));
            }
            evaluate(trial,&mu,&i_found);
            accept_it=0;
            if(i_found>=0){
                if(mu<_chifn->get_fn(i0)){
                    accept_it=1;
                }
                else{
                    roll=_chifn->random_double();
                    ratio=exp(-0.5*(mu-_chifn->get_fn(i0)));
                    if(roll<ratio){
                        accept_it=1;
                    }
                }
            }

            if(accept_it==1){
                accepted++;
                _minimizers.set(ip,i_found);
            }
            else{
                rejected++;
            }
        }
    }

    double f_min=exception_value;
    for(i=0;i<_minimizers.get_dim();i++){
        mu=_chifn->get_fn(_minimizers.get_data(i));
        if(mu<f_min){
            f_min=mu;
        }
    }

    printf("    refined to %e %d %d %e\n",chimin(),accepted,rejected,f_min);

}

void dalex::compass_search(ellipse &ee){
    array_1d<double> center;
    center.set_name("compass_search_center");
    array_1d<double> dir;
    dir.set_name("compass_search_dir");
    int i_found;
    double mu;
    int i;
    for(i=0;i<_chifn->get_dim();i++){
        center.set(i,ee.center(i));
    }
    evaluate(center,&mu,&i_found);
    if(mu>target()){
        printf("   cannot do compass: %e\n",mu);
        return;
    }

    int i_dim;
    double sgn;
    for(i_dim=0;i_dim<_chifn->get_dim();i_dim++){
        for(sgn=-1.0;sgn<1.1;sgn+=2.0){
            for(i=0;i<_chifn->get_dim();i++){
                dir.set(i,sgn*ee.bases(i_dim,i));
            }
            bisection(i_found, dir, target(), 0.01);
        }
    }
}
