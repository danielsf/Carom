#include "dalex.h"

void dalex::write_to_log(char *msg){
    if(_log_file_name[0]==0){
        printf("%s",msg);
        return;
    }

    FILE *log_file;
    log_file = fopen(_log_file_name, "a");
    fprintf(log_file,"%s",msg);
    fclose(log_file);
}

void dalex::build(chisq_wrapper *cc){

    _chifn=cc;

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
    array_1d<int> to_use,explorer_dex;
    assess_good_points();

    int has_explored=0;

    while(_chifn->get_pts()<_limit){
        _chifn->set_search_type(_type_refine);
        iterate_on_minimum();

        _chifn->set_search_type(_type_tendril);
        tendril_search();
        _update_good_points();
    }

}


void dalex::simplex_search(int i_specified){

    printf("\n    doing dalex_simplex %e -- %e %e\n",chimin(),target(),_target_factor);
    safety_check("simplex_search");
    array_1d<double> min,max;
    min.set_name("dalex_simplex_min");
    max.set_name("dalex_simplex_max");
    array_2d<double> seed;
    array_1d<int> seed_dex;
    seed_dex.set_name("dalex_simplex_seed_dex");
    seed.set_name("dalex_simplex_seed");
    array_1d<double> trial,grad;

    int i,j;
    int pt_0=_chifn->get_pts();

    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,0.0);
        max.set(i,_chifn->get_characteristic_length(i));
    }

    ellipse local_ellipse;
    array_2d<double> ellipse_pts;
    ellipse_pts.set_name("simplex_search_ellipse_pts");
    int thin_to=5000;
    int n_thin=-1;
    if(_good_points.get_dim()>2*thin_to){
        n_thin=_good_points.get_dim()/thin_to;
    }
    for(i=0;i<_good_points.get_dim();i++){
        if(n_thin<0 || i%n_thin==0){
            ellipse_pts.add_row(_chifn->get_pt(_good_points.get_data(i)));
        }
    }
    local_ellipse.build(ellipse_pts);

    array_1d<double> dir;
    dir.set_name("simplex_min_dir");
    double target_chi;
    double sgn;
    int pos_ct=0;
    int neg_ct=0;
    int k;
    double dot_product;
    for(i=0;seed.get_rows()<_chifn->get_dim();i++){
        if(i<_chifn->get_dim()){
            for(j=0;j<_chifn->get_dim();j++){
                dir.set(j,local_ellipse.bases(i,j));
            }
        }
        else{
            for(j=0;j<_chifn->get_dim();j++){
                 dir.set(j,_chifn->random_double());
            }
            dir.normalize();
        }
        target_chi=_chifn->chimin()+0.25*_chifn->get_deltachi();
        target_chi += 0.25*_chifn->random_double()*_chifn->get_deltachi();
        pos_ct=0;
        neg_ct=0;
        for(j=0;j<ellipse_pts.get_rows();j++){
            dot_product=0.0;
            for(k=0;k<_chifn->get_dim();k++){
                dot_product+=dir.get_data(k)*(ellipse_pts.get_data(j,k)-_chifn->get_pt(mindex(),k));
            }
            if(dot_product<0.0){
                neg_ct++;
            }
            else{
                pos_ct++;
            }
        }

        if(pos_ct>neg_ct){
            sgn=-1.0;
        }
        else{
            sgn=1.0;
        }

        for(j=0;j<_chifn->get_dim();j++){
            dir.set(j,dir.get_data(j)*sgn);
        }
        j=bisection(mindex(),dir,target_chi,0.05*_chifn->get_deltachi());
        if(j>=0 && seed_dex.contains(j)==0){
            seed_dex.add(j);
            seed.add_row(_chifn->get_pt(j));
        }
    }

    array_1d<double> avg_pt;
    avg_pt.set_name("simplex_search_avg_pt");
    int i_avg;
    double mu_avg;
    for(i=0;i<_chifn->get_dim();i++){
        avg_pt.set(i,0.0);
        for(j=0;j<seed.get_rows();j++){
            avg_pt.add_val(i,seed.get_data(j,i));
        }
        avg_pt.divide_val(i,float(seed.get_rows()));
    }
    evaluate(avg_pt, &mu_avg, &i_avg);
    seed.add_row(avg_pt);

    simplex_minimizer ffmin;
    ffmin.set_bases(_basis_vectors);
    ffmin.set_chisquared(_chifn);
    ffmin.set_minmax(min,max);
    ffmin.set_dice(_chifn->get_dice());
    ffmin.use_gradient();
    ffmin.find_minimum(seed,trial);

    printf("    after dalex_simplex chimin %e\n",chimin());
    _simplex_mindex=mindex();

    array_1d<double> min_pt;
    min_pt.set_name("dalex_simplex_min_pt");
    ffmin.get_minpt(min_pt);
    double mu;
    int i_found;
    evaluate(min_pt, &mu, &i_found);

    _update_good_points(pt_0);
    printf("done with this\n");

}


int dalex::bisection(int ilow, const array_1d<double> &dir, double local_target, double tol){
    safety_check("bisection(int, arr)");
    array_1d<double> trial_high;
    trial_high.set_name("dalex_bisection_trial_high");
    double rr=1.0;
    int ii,i_found;
    double mu=-2.0*exception_value;

    if(_chifn->get_fn(ilow)>local_target){
        printf("WARNING in bisection(i, dir), chifn(ilow) %e > target %e\n",
               _chifn->get_fn(ilow), local_target);
        exit(1);
    }

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
    int i_found,ii,i_low_0,jj,i_trial;
    double mu_found;

    double mu,flow,fhigh;
    evaluate(lowball_in, &flow, &ii);
    evaluate(highball_in, &fhigh, &jj);
    i_low_0=ii;

    if(flow>local_target){
        printf("WARNING flow is greater than target %e %e\n",flow,local_target);
        printf("min %e target %e\n",chimin(),target());
        exit(1);
    }

    i_found=-1;
    mu_found=2.0*exception_value;

    for(ii=0;ii<_chifn->get_dim();ii++){
        lowball.set(ii,lowball_in.get_data(ii));
        highball.set(ii,highball_in.get_data(ii));
    }

    int ct;
    double wgt_low;
    for(ct=0;ct<20 &&
       (ct<5 || fabs(mu_found-local_target)>tol); ct++){

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

        if(i_found<0 || fabs(mu-local_target)<fabs(mu_found-local_target)){
            i_found=i_trial;
            mu_found=mu;
        }
    }

    if(i_found<0){
        i_found=i_low_0;
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
    array_1d<double> focused_min,focused_max;
    focused_min.set_name("find_bases_focused_min");
    focused_max.set_name("find_bases_focused_max");

    array_1d<double> norm_dir;
    norm_dir.set_name("find_bases_norm_dir");
    int i_found_pos;
    int i_found_neg;

    if(changed_bases==1){
        for(i=0;i<_good_points.get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                mu=0.0;
                for(k=0;k<_chifn->get_dim();k++){
                    mu+=_chifn->get_pt(_good_points.get_data(i),k)*_basis_vectors.get_data(j,k);
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
            _basis_associate_norm.set(i,focused_max.get_data(i)-focused_min.get_data(i));
            if(_basis_associate_norm.get_data(i)<1.0e-20){
                printf("WARNING basis associate norm %d %e\n",
                i,_basis_associate_norm.get_data(i));
            }
        }

        for(i=0;i<_chifn->get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                norm_dir.set(j,_basis_vectors.get_data(i,j));
            }
            i_found_pos=bisection(mindex(),norm_dir,target(),0.01);
            for(j=0;j<_chifn->get_dim();j++){
                norm_dir.set(j,-1.0*_basis_vectors.get_data(i,j));
            }
            i_found_neg=bisection(mindex(),norm_dir,target(),0.01);
            if(i_found_pos==i_found_neg && i_found_pos!=mindex()){
                printf("Something is wrong ipos %d ineg %d mindex %d\n",
                i_found_pos,i_found_neg,mindex());
                exit(1);
            }
            mu=cardinal_distance(i_found_pos,i_found_neg);
            if(mu>1.0e-20){
                _basis_norm.set(i,0.5*mu);
            }
            else{
                _basis_norm.set(i,1.0);
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
                                   ellipse_list &exclusion_zones, int *i_next,
                                   int use_relative_norm){

    char log_message[letters];

    safety_check("simplex_boundary_search");
    sprintf(log_message,"\ndoing dalex.simplex_boundary_search() %d\n",_chifn->get_pts());
    write_to_log(log_message);
    int pt_start=_chifn->get_pts();
    assess_good_points();

    int i_node,i_pt;
    int i,j,k;
    double xmin,xmax,xx;

    int n_good_0=_good_points.get_dim();

    array_1d<int> associates;
    associates.set_name("dalex_simplex_boundary_associates");
    int i_start;
    int n_thin=-1;
    int thin_ct=0;
    int ip,io;

    double mu;
    array_1d<double> trial;
    trial.set_name("simplex_boundary_trial");

    array_1d<int> end_points_checked;
    array_1d<double> end_points_chi;
    end_points_checked.set_name("end_points_checked");
    end_points_chi.set_name("end_points_chi");

    int accept_it;
    int associate_iterations=0;
    int thin_to=20000;

    accept_it=0;
    while(accept_it==0){
        associate_iterations++;
        if(associate_iterations>2){
            printf("WARNING associate_iterations %d\n",associate_iterations);
            exit(1);
        }
        for(i=0;i<_tendril_path.get_rows();i++){
            ip=_tendril_path.get_data(i,0);
            io=_tendril_path.get_data(i,1);
            if(end_points_checked.contains(io)==0){
                for(j=0;j<_chifn->get_dim();j++){
                    trial.set(j,0.5*(_chifn->get_pt(specified,j)+_chifn->get_pt(io,j)));
                }

                evaluate(trial,&mu,&k);
                end_points_checked.add(io);
                end_points_chi.add(mu);
            }
            else{
                for(j=0;j<end_points_checked.get_dim();j++){
                    if(end_points_checked.get_data(j)==io){
                        mu=end_points_chi.get_data(j);
                        break;
                    }
                }
            }

            if(mu<target() || specified<=0){
                if(associates.contains(io)==0){
                    associates.add(io);
                }
            }

            if(associates.contains(ip)==0){
                if(specified<=0 || _has_struck==0 || mu<target()){
                    thin_ct++;
                    if(n_thin<0 || thin_ct%n_thin==0){
                        associates.add(ip);
                    }
                }
            }
        }

        if(associates.get_dim()-end_points_checked.get_dim()<2*thin_to){
            accept_it=1;
        }
        else{
            n_thin=associates.get_dim()/thin_to;
            if(n_thin==1){
                n_thin=2;
            }
            thin_ct=0;
            associates.reset();
        }
    }

    int associates_from_tendril=1;
    if(associates.get_dim()==0){
        associates_from_tendril=0;

        if(_good_points.get_dim()>thin_to){
            n_thin=_good_points.get_dim()/thin_to;
            if(n_thin==1){
                n_thin=2;
            }
        }
        for(i=0;i<_good_points.get_dim();i++){
            if(i%n_thin==0){
                associates.add(_good_points.get_data(i));
            }
        }
    }

    if(i_origin>=0){
        if(associates.contains(i_origin)==0){
            associates.add(i_origin);
        }
        if(associates.contains(specified)==0){
            associates.add(specified);
        }
    }

    cost_fn dchifn(_chifn,associates);
    if(_scalar_norm>0.0){
        dchifn.freeze_norm(_scalar_norm);
    }
    dchifn.set_envelope(0.25*(target()-chimin()));
    dchifn.set_bases(_basis_vectors);
    printf("in simplex_boundary_search\n");

    sprintf(log_message,"    associates %d path %d\n", associates.get_dim(),_tendril_path.get_rows());
    write_to_log(log_message);

    simplex_minimizer ffmin;
    ffmin.set_bases(_basis_vectors);
    ffmin.set_chisquared(&dchifn);
    array_1d<double> min,max;
    min.set_name("dalex_simplex_search_min");
    max.set_name("dalex_simplex_search_min");

    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,0.0);
        max.set(i,_chifn->get_characteristic_length(i));
    }

    ffmin.set_minmax(min,max);

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
    if(i_origin>=0 && associates_from_tendril==1){
        for(i=0;i<associates.get_dim();i++){
            ellipse_pts.add_row(_chifn->get_pt(associates.get_data(i)));
        }
    }

    array_1d<double> distances,distances_sorted;
    array_1d<int> distances_dex;
    distances.set_name("simp_bou_dist");
    distances_sorted.set_name("simp_bou_dist_sorted");
    distances_dex.set_name("simp_bou_dist_dex");
    if(ellipse_pts.get_rows()<2*_chifn->get_dim()){
        printf("    need to build ellipse out of nearest good points\n");
        for(i=0;i<_good_points.get_dim();i++){
            distances.add(cardinal_distance(specified,_good_points.get_data(i)));
            distances_dex.add(_good_points.get_data(i));
        }
        sort(distances,distances_sorted,distances_dex);
        ellipse_pts.reset_preserving_room();
        for(i=0;i<2*_chifn->get_dim();i++){
            ellipse_pts.add_row(_chifn->get_pt(distances_dex.get_data(i)));
        }
    }

    dummy_ellipse.build(ellipse_pts);
    write_to_log("built ellipse\n");

    int i_anchor;
    array_1d<double> base_dir;
    base_dir.set_name("dalex_simplex_boundary_base_dir");
    if(i_origin>=0){
        for(i=0;i<_chifn->get_dim();i++){
            base_dir.set(i,_chifn->get_pt(specified,i)-_chifn->get_pt(i_origin,i));
        }
    }
    else{
        for(i=0;i<_chifn->get_dim();i++){
            base_dir.set(i,_chifn->get_pt(specified,i)-_chifn->get_pt(mindex(),i));
        }
        if(_chifn->get_fn(specified)>target()){
            for(i=0;i<_chifn->get_dim();i++){
                base_dir.multiply_val(i,-1.0);
            }
        }
    }

    double base_dir_norm = base_dir.normalize();
    sprintf(log_message,"base_dir_norm %e\n",base_dir_norm);
    write_to_log(log_message);

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

    double local_target,tt;
    double dd,avg_dd;
    int bad_iterations;

    if(specified>=0){
        i_anchor=specified;
        printf("    anchord dchifn %e\n",dchifn(_chifn->get_pt(i_anchor)));
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
                i_bisect1=i_anchor;
                tt=local_target;
                bad_iterations=0;
                while(i_bisect1==i_anchor && bad_iterations<10){
                    i_bisect1=bisection(i_anchor,bisect_dir,tt,0.001);
                    tt+=0.5*(target()-chimin());
                    bad_iterations++;
                }

                if(i_bisect1==i_anchor){
                    bad_iterations=0;
                    while(i_bisect1==i_anchor){
                        for(j=0;j<_chifn->get_dim();j++){
                            bisect_dir.set(j,normal_deviate(_chifn->get_dice(),0.0,1.0));
                        }
                        i_bisect1=bisection(i_anchor,bisect_dir,local_target,0.001);
                        bad_iterations++;
                        if(bad_iterations%50==0){
                            sprintf(log_message,"    bad_iterations on i_bisect1 %d -- %e %e\n",
                            bad_iterations,_chifn->get_fn(i_anchor),local_target);
                            write_to_log(log_message);
                        }
                        if(bad_iterations>200){
                            printf("WARNING 200 bad iterations on i_bisect1\n");
                            exit(1);
                        }
                    }
                }

                for(j=0;j<_chifn->get_dim();j++){
                    bisect_dir.set(j,-1.0*dummy_ellipse.bases(i,j)+base_dir.get_data(j));
                }

                i_bisect2=i_anchor;
                tt=local_target;
                bad_iterations=0;
                while(i_bisect2==i_anchor && bad_iterations<10){
                    i_bisect2=bisection(i_anchor,bisect_dir,tt,0.001);
                    tt+=0.5*(target()-chimin());
                    bad_iterations++;
                }

                if(i_bisect2==i_anchor){
                    bad_iterations=0;
                    while(i_bisect2==i_anchor){
                        for(j=0;j<_chifn->get_dim();j++){
                            bisect_dir.set(j,normal_deviate(_chifn->get_dice(),0.0,1.0));
                        }
                        i_bisect2=bisection(i_anchor,bisect_dir,local_target,0.001);
                        bad_iterations++;
                        if(bad_iterations%50==0){
                            sprintf(log_message,"    bad_iterations on i_bisect2 %d -- %e %e\n",
                            bad_iterations,_chifn->get_fn(i_anchor), local_target);
                            write_to_log(log_message);
                        }
                        if(bad_iterations>200){
                            printf("WARNING 200 bad iterations on i_bisect2\n");
                            exit(1);
                        }
                    }
                }

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
                    printf("could not find good i_bisect\n");
                    exit(1);
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

            if(seed.get_rows()<_chifn->get_dim()){
                sprintf(log_message,"after looping dim; seed is %d\n",seed.get_rows());
                write_to_log(log_message);
                exit(1);
            }

            if(seed.get_rows()==_chifn->get_dim()){
                avg_dd=0.0;
                for(i=0;i<_chifn->get_dim();i++){
                    dd=0.0;
                    for(j=0;j<_chifn->get_dim();j++){
                        dd+=power(_chifn->get_pt(i_anchor,j)-seed.get_data(i,j),2);
                    }
                    avg_dd+=sqrt(dd);
                }
                avg_dd/=double(_chifn->get_dim());
                for(i=0;i<_chifn->get_dim();i++){
                    avg_pt.set(i,_chifn->get_pt(specified,i)+avg_dd*base_dir.get_data(i));
                }
                evaluate(avg_pt,&mu1,&i);
                printf("    fn_avg %e\n",dchifn(avg_pt));
                seed.add_row(avg_pt);
            }
        }
    }
    else{
        printf("calling _explorers.get_seed(); did not expect that\n");
        exit(1);
    }

    sprintf(log_message,"got seeds; fn anchor %e; %d %d %d\n",
    _chifn->get_fn(i_anchor),specified,i_origin,i_anchor);
    write_to_log(log_message);

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
    sprintf(log_message,"    starting from %e; %e; %e\n",start_min,mu,start_max);
    write_to_log(log_message);
    sprintf(log_message,"    starting from (delta) %e; %e; %e\n",start_min-target(),mu-target(),start_max-target());
    write_to_log(log_message);

    array_1d<double> minpt;
    minpt.set_name("dalex_simplex_search_minpt");

    // loop over points, calling dchifn so that they
    // get stored in dchifn's cache
    for(i=pt_start;i<_chifn->get_pts();i++){
        mu=dchifn(_chifn->get_pt(i));
    }

    ffmin.find_minimum(seed,minpt);

    sprintf(log_message,"    found minimum -- chimin %e\n",chimin());
    write_to_log(log_message);

    array_1d<int> path_row;
    path_row.set_name("path_row");

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
        evaluate(minpt,&mu,i_next);
    }

    write_to_end_pt_file(i_next[0]);

    int pre_fill=_chifn->get_pts();
    array_1d<double> perp_dir;
    array_1d<double> trial_dir;
    perp_dir.set_name("dalex_simp_bou_perp_dir");
    trial_dir.set_name("dalex_simp_bou_trial_dir");
    double base_norm;
    int iteration;
    double step;
    double rr;
    int n_good_fill=0;
    if(i_next[0]!=specified){
        for(i=0;i<_chifn->get_dim();i++){
            base_dir.set(i,_chifn->get_pt(i_next[0],i)-_chifn->get_pt(specified,i));
        }
        base_norm=base_dir.normalize();
        for(iteration=0;iteration<2*_chifn->get_dim();iteration++){
            for(i=0;i<_chifn->get_dim();i++){
                perp_dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
            }
            component=0.0;
            for(i=0;i<_chifn->get_dim();i++){
                component+=perp_dir.get_data(i)*base_dir.get_data(i);
            }
            for(i=0;i<_chifn->get_dim();i++){
                perp_dir.subtract_val(i,component*base_dir.get_data(i));
            }
            perp_dir.normalize();
            rr=_chifn->random_double();
            for(i=0;i<_chifn->get_dim();i++){
                trial_dir.set(i,base_dir.get_data(i)+rr*perp_dir.get_data(i));
            }
            trial_dir.normalize();
            for(step=0.1*base_norm;step<1.05*base_norm;step+=0.1*base_norm){
                for(i=0;i<_chifn->get_dim();i++){
                    trial.set(i,_chifn->get_pt(specified,i)
                               +step*trial_dir.get_data(i));
                }
                evaluate(trial,&mu,&i);
                if(mu<target()){
                    n_good_fill++;
                }
            }
        }
    }

    sprintf(log_message,"    filling called %d good %d\n",_chifn->get_pts()-pre_fill,n_good_fill);
    write_to_log(log_message);

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

    sprintf(log_message,"    adjusted %e from %e\n",
    dchifn(_chifn->get_pt(i_next[0])),_chifn->get_fn(i_next[0]));
    write_to_log(log_message);

    sprintf(log_message,"    min is %e target %e\n",chimin(),target());
    write_to_log(log_message);

    int path_start,path_end,path_mid;
    path_start=-1;
    path_end=-1;
    path_mid=-1;
    array_1d<double> path_dist,path_dist_sorted;
    array_1d<int> path_dist_dex;
    path_dist.set_name("dalex_simp_bou_path_dist");
    path_dist_sorted.set_name("dalex_simp_bou_path_dist_sorted");
    path_dist_dex.set_name("dalex_simp_bou_path_dist_dex");

    for(i=pt_start;i<_chifn->get_pts();i++){
        if(_chifn->get_fn(i)<target()){
            mu=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                mu+=(_chifn->get_pt(i,j)-_chifn->get_pt(specified,j))*base_dir.get_data(j);
            }
            path_dist.add(mu);
            path_dist_dex.add(i);
        }
    }
    if(path_dist.get_dim()>1){
        sort(path_dist, path_dist_sorted, path_dist_dex);
    }

    if(path_dist.get_dim()>0){
        path_start=path_dist_dex.get_data(0);
    }
    if(path_dist.get_dim()>1){
        path_end=path_dist_dex.get_data(path_dist.get_dim()-1);
    }
    if(path_dist.get_dim()>2){
        path_mid=path_dist_dex.get_data(path_dist.get_dim()/2);
    }

    write_to_end_pt_file(path_mid);

    double dd_start,dd_mid,dd_end;

    if(path_start>=0 || path_end>=0 || path_mid>0){
        sprintf(log_message,"going to try to add to path %d %d\n",
        path_start,_chifn->get_pts());
        write_to_log(log_message);
        sprintf(log_message,"start %d end %d\n",path_start,path_end);
        write_to_log(log_message);

        for(i=path_start;i<_chifn->get_pts();i++){
            dd_start=-1.0;
            dd_end=-1.0;
            dd_mid=-1.0;
            if(_tendril_path.get_rows()==0){
                _tendril_path.set_cols(2);
            }

            if(_chifn->get_fn(i)<_chifn->target() &&
               _chifn->get_search_type_log(i)==_type_tendril){

                path_row.set(0,i);

                if(path_start>=0){
                    dd_start=distance(i,path_start);
                }
                if(path_end>=0){
                    dd_end=distance(i,path_end);
                }
                if(path_mid>=0){
                    dd_mid=distance(i,path_mid);
                }

                if(dd_start<0.0 && dd_end<0.0 && dd_mid<0.0){
                    printf("WARNING constructing path but all dd are <0\n");
                    exit(1);
                }

                if(dd_start<0.0){
                    if(dd_end<0.0){
                        path_row.set(1,path_mid);
                    }
                    else if(dd_mid<0.0){
                        path_row.set(1,path_end);
                    }
                    else if(dd_mid<dd_end){
                        path_row.set(1,path_mid);
                    }
                    else{
                        path_row.set(1,path_end);
                    }
                }
                else if(dd_mid<0.0){
                    if(dd_start<0.0){
                        path_row.set(1,path_end);
                    }
                    else if(dd_end<0.0){
                        path_row.set(1,path_start);
                    }
                    else if(dd_start<dd_end){
                        path_row.set(1,path_start);
                    }
                    else{
                        path_row.set(1,path_end);
                    }
                }
                else if(dd_end<0.0){
                    if(dd_start<0.0){
                        path_row.set(1,path_mid);
                    }
                    else if(dd_mid<0.0){
                        path_row.set(1,path_start);
                    }
                    else if(dd_start<dd_mid){
                        path_row.set(1,path_start);
                    }
                    else{
                        path_row.set(1,path_mid);
                    }
                }
                else{
                    if(dd_start<=dd_end && dd_start<=dd_mid){
                        path_row.set(1,path_start);
                    }
                    else if(dd_mid<=dd_end && dd_mid<=dd_start){
                        path_row.set(1,path_mid);
                    }
                    else{
                        path_row.set(1,path_end);
                    }
                }

                if(path_row.get_data(1)<0){
                    printf("WARNING somehow path_row(1)<0\n");
                    exit(1);
                }
                _tendril_path.add_row(path_row);
            }
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

int dalex::_exploration_simplex(int i1, int i0, array_1d<int> &associates){

    array_2d<double> associate_pts;
    associate_pts.set_name("dalex_exp_simp_associate_pts");
    int i;
    for(i=0;i<associates.get_dim();i++){
        associate_pts.add_row(_chifn->get_pt(associates.get_data(i)));
    }
    cost_fn dchifn(_chifn, associates);
    if(_scalar_norm>0.0){
        dchifn.freeze_norm(_scalar_norm);
    }
    dchifn.set_bases(_basis_vectors);
    printf("in _exploration_simplex\n");

    array_2d<double> seed;
    seed.set_name("dalex_exp_sim_seed");
    array_1d<double> min,max;
    min.set_name("dalex_exp_simp_min");
    max.set_name("dalex_exp_simp_max");

    array_1d<double> dir;
    dir.set_name("dalex_exp_sim_dir");
    array_2d<double> dir_log;
    dir_log.set_name("dalex_exp_sim_dir_log");
    for(i=0;i<_chifn->get_dim();i++){
        dir.set(i,_chifn->get_pt(i1,i)-_chifn->get_pt(i0,i));
    }
    double dir_norm;
    dir_norm=dir.normalize();
    array_1d<double> newpt,newdir;
    newpt.set_name("dalex_exp_simp_newpt");
    newdir.set_name("dalex_exp_simp_newdir");
    for(i=0;i<_chifn->get_dim();i++){
        newpt.set(i,_chifn->get_pt(i1,i)+0.1*dir_norm*dir.get_data(i));
    }
    seed.add_row(newpt);
    dir_log.add_row(dir);
    double component;
    int j;
    while(seed.get_rows()<_chifn->get_dim()){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        for(i=0;i<dir_log.get_rows();i++){
            component=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                component+=dir.get_data(j)*dir_log.get_data(i,j);
            }
            for(j=0;j<_chifn->get_dim();j++){
                dir.subtract_val(j,component*dir_log.get_data(i,j));
            }
        }
        component=dir.normalize();
        if(component>1.0e-20){
            for(i=0;i<_chifn->get_dim();i++){
                newdir.set(i,dir.get_data(i)+dir_log.get_data(0,i));
            }
            newdir.normalize();
            for(i=0;i<_chifn->get_dim();i++){
                newpt.set(i,_chifn->get_pt(i1,i)+0.1*dir_norm*newdir.get_data(i));
            }
            seed.add_row(newpt);
            dir_log.add_row(dir);
        }
    }
    for(i=0;i<_chifn->get_dim();i++){
        newpt.set(i,0.0);
    }
    for(i=0;i<seed.get_rows();i++){
        for(j=0;j<_chifn->get_dim();j++){
            newpt.add_val(j,seed.get_data(i,j));
        }
    }
    for(i=0;i<_chifn->get_dim();i++){
        newpt.divide_val(i,double(seed.get_rows()));
    }
    seed.add_row(newpt);

    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,0.0);
        max.set(i,_chifn->get_characteristic_length(i));
    }

    simplex_minimizer ffmin;
    ffmin.set_bases(_basis_vectors);
    ffmin.set_chisquared(&dchifn);
    ffmin.set_dice(_chifn->get_dice());
    ffmin.set_minmax(min,max);
    ffmin.use_gradient();
    ffmin.find_minimum(seed,newpt);
    double mu;
    int i_found;
    evaluate(newpt,&mu,&i_found);
    printf("exp simplex found %e\n",mu);
    return i_found;

}


void dalex::find_tendril_candidates(double factor_in){

    _particle_candidates.reset_preserving_room();
    _origin_candidates.reset_preserving_room();

    char log_message[letters];

    write_to_log("finding tendril candidates\n");

    printf("finding tendril candidates\n");

    array_1d<int> associates;
    associates.set_name("find_tendrils_associates");
    array_1d<int> associates_raw;
    associates_raw.set_name("find_tendril_associates_raw");
    array_2d<double> ellipse_pts;
    ellipse_pts.set_name("find_tendrils_ellipse_pts");
    array_1d<double> center;
    center.set_name("find_tendrils_center");
    int i_center;
    int i,j;
    for(i=0;i<_chifn->get_pts();i++){
        if(_chifn->get_fn(i)<target()){
            ellipse_pts.add_row(_chifn->get_pt(i));
            associates.add(i);
            associates_raw.add(i);
        }
    }
    if(ellipse_pts.get_rows()==0){
        printf("CANNOT find tendrils; not good points\n");
        exit(1);
    }
    for(i=0;i<_chifn->get_dim();i++){
        center.set(i,_chifn->get_pt(mindex(),i));
    }
    i_center=mindex();

    int n_thin;

    if(associates.get_dim()>20000){
        n_thin=associates_raw.get_dim()/20000;
        if(n_thin==1){
            n_thin=2;
        }
        associates.reset_preserving_room();
        for(i=0;i<associates_raw.get_dim();i++){
            if(i%n_thin==0){
                associates.add(associates_raw.get_data(i));
            }
        }
    }

    ellipse good_ellipse;
    good_ellipse.use_geo_center();
    good_ellipse.build(ellipse_pts);

    cost_fn dchifn(_chifn,associates,1);

    if(_scalar_norm>0.0){
        dchifn.freeze_norm(_scalar_norm);
    }

    dchifn.set_bases(_basis_vectors);
    printf("in find_tendril_candidates\n");
    double envelope=0.25*(target()-chimin());
    if(envelope<2.0){
        envelope=2.0;
    }
    dchifn.set_envelope(envelope);

    array_1d<double> min,max;
    min.set_name("find_tendrils_min");
    max.set_name("find_tendrils_max");
    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,0.0);
        max.set(i,_chifn->get_characteristic_length(i));
    }

    simplex_minimizer ffmin;
    ffmin.set_bases(_basis_vectors);
    ffmin.set_minmax(min,max);
    ffmin.set_chisquared(&dchifn);

    array_1d<double> trial,minpt;
    trial.set_name("find_tendrils_trial");
    minpt.set_name("find_tendrils_minpt");
    double mu;
    int i_found;

    array_1d<int> particles;
    particles.set_name("find_tendrils_particles");
    array_1d<double> fn_val;
    fn_val.set_name("find_tendrils_fn_val");
    array_1d<int> fn_val_dex;
    fn_val_dex.set_name("find_tendrils_fn_val_dex");

    array_2d<double> seed;
    seed.set_name("find_tendrils_seed");

    int idim,jdim;
    double sgn,cost;
    int pt_start;
    array_1d<int> at_least_one;
    at_least_one.set_name("find_tendril_at_least_one");
    double factor;
    int try_again;
    array_1d<double> retry_pt;
    retry_pt.set_name("find_tendril_retry_pt");
    double retry_dist,retry_dist_min;
    int i_retry;

    array_1d<double> bisection_dir,delta_dir;
    bisection_dir.set_name("find_tendril_candidates_bisection_dir");
    delta_dir.set_name("find_tendril_candidates_delta_dir");
    double bisection_target;
    array_1d<int> seed_int_list;
    seed_int_list.set_name("find_tendril_candidates_seed_int_list");

    int strikes;
    int failures;
    int i_start;
    double mu_retry;
    double mu_found;

    array_1d<double> spherical_dir;
    spherical_dir.set_name("spherical_dir");

    printf("n associates %d\n",associates.get_dim());
    for(idim=0;idim<2*_chifn->get_dim();idim++){
            pt_start=_chifn->get_pts();
            seed.reset_preserving_room();
            i_found=-1;
            factor=factor_in;
            for(i=0;i<_chifn->get_dim();i++){
                spherical_dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
            }
            spherical_dir.normalize();
            while(i_found<0){
                printf("    running with factor %e\n",factor);

                for(i=0;i<_chifn->get_dim();i++){
                    trial.set(i,good_ellipse.geo_center(i));
                }
                for(i=0;i<_chifn->get_dim();i++){
                    for(j=0;j<_chifn->get_dim();j++){
                        trial.add_val(j,factor*good_ellipse.radii(i)*spherical_dir.get_data(i)*good_ellipse.bases(i,j));
                    }
                }
                evaluate(trial,&mu,&i_found);
                printf("    first point is %d\n",i_found);
                if(factor>2.0){
                    factor-=1.0;
                }
                else{
                    factor*=0.5;
                }
            }

            if(_chifn->get_fn(i_found)>chimin()+0.1*_chifn->get_deltachi()){
                // just naively using a multiple of the ellipse radius worked
                write_to_end_pt_file(i_found);
                seed.add_row(trial);
                for(jdim=0;jdim<_chifn->get_dim();jdim++){
                    if(spherical_dir.get_data(jdim)<0.0){
                        sgn=-1.0;
                    }
                    else{
                        sgn=1.0;
                    }
                    for(i=0;i<_chifn->get_dim();i++){
                       bisection_dir.set(i,sgn*good_ellipse.bases(jdim,i));
                    }
                    i=bisection(i_found,bisection_dir,_chifn->get_fn(i_found)+_chifn->get_deltachi(),0.1);
                    seed.add_row(_chifn->get_pt(i));
                }
            }
            else{
                // alternative scheme for finding seed
                write_to_log("doing alternative seed in find_tendril_candidates\n");
                seed_int_list.reset_preserving_room();
                while(seed.get_rows()<_chifn->get_dim()+1){
                    for(i=0;i<_chifn->get_dim();i++){
                        bisection_dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
                        trial.set(i,bisection_dir.get_data(i));
                    }
                    trial.normalize();
                    bisection_target=_chifn->get_fn(i_center)+10.0*_chifn->get_deltachi();
                    i_found=bisection(i_center, bisection_dir, bisection_target, 0.1*_chifn->get_deltachi());
                    if(i_found>=0 && seed_int_list.contains(i_found)==0 &&
                       _chifn->get_fn(i_found)>target()){
                        sprintf(log_message,"i_found %d chisq %e min %e\n",
                                i_found,_chifn->get_fn(i_found),chimin());
                        write_to_log(log_message);
                        if(seed.get_rows()==0){
                            write_to_end_pt_file(i_found);
                        }
                        seed.add_row(_chifn->get_pt(i_found));
                        seed_int_list.add(i_found);
                    }
                }
                write_to_log("    done doing alternative seed in find_tendril_candidates\n");
            }

            evaluate(seed(0), &mu, &i_start);
            if(i_start<0){
                printf("WARNING i_start %d\n",i_start);
                exit(1);
            }
            ffmin.find_minimum(seed,minpt);
            evaluate(minpt,&mu,&i_found);

            if(_chifn->get_fn(i_found)>target()){
                try_again=1;
            }
            else{
                try_again=0;
            }

            strikes=0;
            failures=0;
            while(try_again==1){

                if(_limit>0 && _chifn->get_pts()>_limit){
                    _chifn->write_pts();
                    return;
                }

                sprintf(log_message,"    trying again because %e > %e\n",
                        _chifn->get_fn(i_found),target());
                write_to_log(log_message);
                seed.reset_preserving_room();

                i_retry=-1;
                for(i=0;i<_good_points.get_dim();i++){
                    j=_good_points.get_data(i);
                    retry_dist=distance(i_start,j);
                    if(i==0 || retry_dist<retry_dist_min){
                        retry_dist_min=retry_dist;
                        i_retry=j;
                    }
                }
                if(i_retry<=0){
                    printf("CANNOT retry; i_retry<0\n");
                    exit(1);
                }

                for(i=0;i<_chifn->get_dim();i++){
                    retry_pt.set(i,0.5*(_chifn->get_pt(i_retry,i)+
                                        _chifn->get_pt(i_start,i)));
                }

                evaluate(retry_pt,&mu_retry,&i_retry);
                mu_found=_chifn->get_fn(i_found);

                if(mu_found<target()+_chifn->get_deltachi() ||
                   dchifn(_chifn->get_pt(i_found))<target()){

                    seed.add_row(_chifn->get_pt(i_found));
                    i_start=i_found;
                    printf("    using i_found %e\n",mu_found);

                }
                else{
                    seed.add_row(retry_pt);
                    i_start=i_retry;
                    printf("    using i_retry %e %e\n",
                    mu_retry,distance(mindex(),i_retry));
                }

                if(i_start<0){
                    printf("WARNING retry point %d\n",i_start);
                    exit(1);
                }
                for(jdim=0;jdim<_chifn->get_dim();jdim++){
                    for(i=0;i<_chifn->get_dim();i++){
                        bisection_dir.set(i,sgn*good_ellipse.bases(jdim,i));
                    }
                    i=bisection(i_found,bisection_dir,_chifn->get_fn(i_found)+_chifn->get_deltachi(),0.1);
                    seed.add_row(_chifn->get_pt(i));
                }
                ffmin.set_dice(_chifn->get_dice());
                ffmin.use_gradient();
                ffmin.find_minimum(seed,minpt);
                evaluate(minpt,&mu,&i_found);
                ffmin.do_not_use_gradient();
                if(_chifn->get_fn(i_found)<target()){
                    try_again=0;
                }
                if(_chifn->get_fn(i_found)>target() &&
                   dchifn(minpt)<target()){

                    strikes++;
                    if(strikes>=3){
                        try_again=0;
                    }
                }
                failures++;
                if(failures>_chifn->get_dim()/2){
                    try_again=0;
                }
            }

            if(_chifn->get_fn(i_found)<target()){
                at_least_one.add(i_found);
            }

            particles.add(i_found);
            write_to_end_pt_file(i_found);
            cost=dchifn(minpt);
            fn_val.add(cost);


            if(_chifn->get_dim()>9){
                sprintf(log_message,"    got %d chisq %e cost %e\n",
                        particles.get_dim(),mu,cost);
                write_to_log(log_message);

                sprintf(log_message,"    pt %e %e \n\n",
                                    minpt.get_data(6),
                                    minpt.get_data(9));
                write_to_log(log_message);
            }
            else{
                sprintf(log_message,"    got %d chisq %e cost %e\n\n",
                        particles.get_dim(),mu,cost);
                write_to_log(log_message);
            }

            j=0;
            for(i=pt_start;i<_chifn->get_pts();i++){
                if(_chifn->get_fn(i)<target()){
                    associates.add(i);
                    associates_raw.add(i);
                    j++;
                }
            }
            if(j>0){
                if(associates.get_dim()>20000){
                    n_thin=associates_raw.get_dim()/20000;
                    if(n_thin==1){
                        n_thin=2;
                    }
                    associates.reset_preserving_room();
                    for(i=0;i<at_least_one.get_dim();i++){
                        associates.add(at_least_one.get_data(i));
                    }
                    for(i=0;i<associates_raw.get_dim();i++){
                        if(i%n_thin==0 && associates.contains(associates_raw.get_data(i))==0){
                            associates.add(associates_raw.get_data(i));
                        }
                    }
                }
                dchifn.build(_chifn,associates,1);
            }

            if(_limit>0 && _chifn->get_pts()>_limit){
                _chifn->write_pts();
                return;
            }

            sprintf(log_message,"    chimin %e\n",chimin());
            write_to_log(log_message);

    }

    if(_scalar_norm<0.0){
        _scalar_norm=dchifn.scalar_norm();
    }

    for(i=0;i<fn_val.get_dim();i++){
        fn_val_dex.add(i);
    }
    array_1d<double> fn_val_sorted;
    fn_val_sorted.set_name("find_tendrils_fn_val_sorted");
    sort(fn_val,fn_val_sorted,fn_val_dex);

    for(i=0;i<particles.get_dim();i++){
        _particle_candidates.set(i,particles.get_data(i));
        _origin_candidates.set(i,-1);
    }
    sprintf(log_message,"done finding tendrils -- chimin %e\n",chimin());
    write_to_log(log_message);

    array_2d<double> simplex_seed;
    seed.set_name("find_candidates_simplex_seed");
    for(i=0;i<_chifn->get_dim();i++){
        simplex_seed.add_row(_chifn->get_pt(particles.get_data(fn_val_dex.get_data(i))));
    }
    simplex_seed.add_row(_chifn->get_pt(mindex()));

    array_1d<double> new_min_pt;

    simplex_minimizer ffmin_2;
    ffmin_2.set_bases(_basis_vectors);
    ffmin_2.set_chisquared(_chifn);
    ffmin_2.set_minmax(min,max);
    ffmin_2.set_dice(_chifn->get_dice());
    ffmin_2.use_gradient();

    double min_0=chimin();
    ffmin_2.find_minimum(simplex_seed,new_min_pt);
    if(chimin()<min_0){
        iterate_on_minimum();
    }
    _chifn->write_pts();

}


void dalex::get_new_tendril(int *particle, int *origin){
    particle[0]=-1;
    origin[0]=-1;
    int i,j,ip,io;
    int is_outside;
    int old_type;

    cost_fn dchifn;
    if(_scalar_norm>0.0){
        dchifn.freeze_norm(_scalar_norm);
    }
    dchifn.set_bases(_basis_vectors);
    printf("in get_new_tendril\n");
    array_1d<int> associates;
    array_1d<double> cost_val,cost_val_sorted;
    array_1d<int> cost_val_dex;
    associates.set_name("new_tendril_associates");
    cost_val.set_name("new_tendril_cost_val");
    cost_val_sorted.set_name("new_tendril_cost_val_sorted");
    cost_val_dex.set_name("new_tendril_cost_val_dex");

    int invalid_candidates;

    while(particle[0]<0){
        invalid_candidates = 0;
        for(i=0;i<_particle_candidates.get_dim();i++){
            if(_particle_candidates.get_data(i)<0){
                invalid_candidates++;
            }
        }
        if(invalid_candidates>=_chifn->get_dim()/2){
            old_type=_chifn->get_search_type();
            _chifn->set_search_type(_type_init_tendril);
            find_tendril_candidates(3.0);
            _chifn->set_search_type(old_type);
        }

        if(_limit>0 && _chifn->get_pts()>_limit){
            _chifn->write_pts();
            return;
        }
        associates.reset_preserving_room();
        cost_val.reset_preserving_room();
        cost_val_sorted.reset_preserving_room();
        cost_val_dex.reset_preserving_room();
        for(i=0;i<_chifn->get_pts();i++){
            if(_chifn->get_fn(i)<target()){
                if(_chifn->get_search_type_log(i)!=_type_init_tendril){
                    associates.add(i);
                }
            }
        }

        if(associates.get_dim()>0){
            dchifn.build(_chifn,associates,1);
            for(i=0;i<_particle_candidates.get_dim();i++){
                if(_particle_candidates.get_data(i)<0){
                    cost_val.add(2.0*exception_value);
                }
                else{
                    cost_val.add(dchifn(_chifn->get_pt(_particle_candidates.get_data(i))));
                }
                cost_val_dex.add(i);
            }
        }
        else{
            // if no associates, just choose particle that is
            // farthest away from distance
            for(i=0;i<_particle_candidates.get_dim();i++){
                if(_particle_candidates.get_data(i)<0){
                   cost_val.add(2.0*exception_value);
                }
                else{
                    cost_val.add(-1.0*_chifn->distance(mindex(),
                                 _particle_candidates.get_data(i)));
                }
                cost_val_dex.add(i);
            }
        }

        sort(cost_val, cost_val_sorted, cost_val_dex);

        for(i=0;i<_chifn->get_dim()/2;i++){
            ip=_particle_candidates.get_data(cost_val_dex.get_data(i));
            io=_origin_candidates.get_data(cost_val_dex.get_data(i));
            if(ip>=0){
                is_outside=1;
                for(j=0;j<_exclusion_zones.ct();j++){
                    if(_exclusion_zones(j)->contains(_chifn->get_pt(ip))==1){
                        is_outside=0;
                        break;
                    }
                }
                if(is_outside==1){
                    particle[0]=ip;
                    origin[0]=io;
                    _particle_candidates.set(cost_val_dex.get_data(i),-1);
                    _origin_candidates.set(cost_val_dex.get_data(i),-1);
                    printf("returning tendril %d %d\n",particle[0],origin[0]);
                    return;
                }
                else{
                    _particle_candidates.set(cost_val_dex.get_data(i),-1);
                    _origin_candidates.set(cost_val_dex.get_data(i),-1);
                }
            }
        }
    }
}


void dalex::init_fill(){
    int path_len=_tendril_path.get_rows();
    int old_type=_chifn->get_search_type();
    _chifn->set_search_type(_type_init);
    find_tendril_candidates(1.0);
    int i;
    for(i=0;i<_particle_candidates.get_dim();i++){
        _particle_candidates.set(i,-1);
        _origin_candidates.set(i,-1);
    }

    _chifn->set_search_type(old_type);
}

void dalex::tendril_search(){

    if(_tendril_init==0){
        init_fill();
        _tendril_init=1;
        if(_limit>0 && _chifn->get_pts()>_limit){
            _chifn->write_pts();
            return;
        }
    }

    int particle, origin;

    get_new_tendril(&particle,&origin);
    _tendril_search(particle, origin);
    _tendril_search(mindex(), particle);

}

void dalex::_tendril_search(int particle, int origin){

    char log_message[letters];
    sprintf(log_message,"running tendril search - %d\n",_chifn->get_pts());
    write_to_log(log_message);

    int pt_start=_chifn->get_pts();
    int pt_prime=_chifn->get_pts();
    array_1d<double> dir,midpt;
    int i_found,i_next;
    double mu;
    dir.set_name("tendril_dir");
    midpt.set_name("tendril_mid_pt");
    int i,j;

    array_2d<double> ellipse_pts;
    ellipse_pts.set_name("tendril_ellipse_pts");
    ellipse local_ellipse;
    int is_a_strike;

    array_1d<int> new_particles,new_origins;
    new_particles.set_name("tendril new_particles");
    new_origins.set_name("new_origins");

    int strikes=0;
    _has_struck=0;

    ellipse_pts.reset_preserving_room();
    local_ellipse.use_geo_center();

    while(strikes<3){

        is_a_strike=0;
        pt_start=_chifn->get_pts();

        is_a_strike=simplex_boundary_search(particle,origin,_exclusion_zones,&i_next,1);

        for(j=pt_start;j<_chifn->get_pts();j++){
            if(_chifn->get_fn(j)<target()){
                ellipse_pts.add_row(_chifn->get_pt(j));
            }
        }

        if(ellipse_pts.get_rows()>2*_chifn->get_dim()){
            local_ellipse.build(ellipse_pts);
        }

        if(is_a_strike==0){
            origin=particle;
            particle=i_next;
            if(_chifn->get_fn(i_next)<target()+1.0e-6){
                strikes=0;
            }
            else{
                is_a_strike=1;
                strikes++;
            }
        }
        else{
             strikes++;
             _has_struck=1;
        }

        printf("pts %d lim %d\n",_chifn->get_pts(),_limit);
        if(_chifn->get_dim()>9)printf("got to %e %e\n",_chifn->get_pt(i_next,6),_chifn->get_pt(i_next,9));
        printf("is a strike: %d; strikes %d\n",is_a_strike,strikes);
        printf("\n");

        if(ellipse_pts.get_rows()>2*_chifn->get_dim()){
            local_ellipse.build(ellipse_pts);
            _exclusion_zones.add(local_ellipse);
            ellipse_pts.reset_preserving_room();
        }

    }


    printf("\n    pts %d limit %d\n",_chifn->get_pts(),_limit);
    printf("    strikeouts %d\n",_strikeouts);
    _chifn->write_pts();
    _update_good_points();
}

void dalex::initialize_min_exploration(){
    _min_explorers.initialize();
    array_1d<double> dir,midpt;
    dir.set_name("dalex_init_min_exp_dir");
    midpt.set_name("dalex_init_min_exp_midpt");
    double mu;
    int i_found;
    int i,j;
    while(_min_explorers.get_n_particles()<2*_chifn->get_dim()){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(), 0.0, 1.0));
        }
        dir.normalize();
        i_found=bisection(mindex(),dir,target(),0.001);
        if(i_found!=mindex()){
            for(i=0;i<_chifn->get_dim();i++){
                midpt.set(i,0.5*(_chifn->get_pt(mindex(),i)+_chifn->get_pt(i_found,i)));
            }
            evaluate(midpt,&mu,&i_found);
            if(i_found!=mindex()){
                _min_explorers.add_particle(_chifn->get_pt(i_found));
            }
        }
    }

}

void dalex::min_explore(int n_particles, int n_steps){
    printf("\nmin exploring\n");
    int pt_0=_chifn->get_pts();

    if(_min_explorers.get_n_particles()==0){
        initialize_min_exploration();
    }

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

    int n_simplex=0;
    int n_start;
    double min_before_simp;
    double d_simp;

    array_1d<double> chimin_arr;
    chimin_arr.set_name("simplex_search_chimin_arr");
    int max_iter=10;
    double min_diff=1.0;

    int old_mindex=-1;
    char log_message[letters];

    double d_min=0.01*_chifn->get_deltachi();
    if(d_min>0.1){
        d_min=0.1;
    }

    int strikes=0;
    int mindex_0 = mindex();

    while(strikes<3){
        min_0=chimin();
        n_start= _chifn->get_pts();
        min_before_simp = chimin();
        simplex_search(mindex());
        n_simplex = _chifn->get_pts()-n_start;
        d_simp = chimin()-min_before_simp;
        min_1=chimin();

        if(mindex() == mindex_0){
            strikes=3;
        }
        else if(min_1<min_0-d_min){
            strikes=0;
        }
        else{
            strikes++;
        }

        sprintf(log_message,"in iterate: min1 %.5e min0 %.5e diff %.2e ",
                min_1,min_0,min_1-min_0);
        write_to_log(log_message);

        sprintf(log_message,"n_simp %d %.2e ",
                n_simplex,d_simp);
        write_to_log(log_message);

        write_to_log("\n");

        if(chimin_arr.get_dim()<max_iter){
            chimin_arr.add(chimin());
        }
        else{
            for(i=0;i<max_iter-1;i++){
                chimin_arr.set(i,chimin_arr.get_data(i+1));
            }
            chimin_arr.set(max_iter-1,chimin());
            if(chimin_arr.get_data(0)-chimin_arr.get_data(max_iter-1)<min_diff){
                break;
            }
        }

    }

    array_2d<int> tendril_cache;
    tendril_cache.set_name("tendril_cache");
    array_1d<int> good_cache;
    good_cache.set_name("good_cache");
    array_1d<int> origins;
    origins.set_name("iterate_min_origins");
    array_2d<double> ellipse_pts;
    ellipse_pts.set_name("iterate_ellipse_pts");
    ellipse local_ellipse;

    if(chimin()<_reset_chimin-_reset_threshold){

        sprintf(log_message, "finding bases: min %e\n",chimin());
        write_to_log(log_message);

        for(i=0;i<_tendril_path.get_rows();i++){
            if(_chifn->get_fn(_tendril_path.get_data(i,0))<target()+1.0e-6){
                tendril_cache.add_row(_tendril_path(i));
            }
        }

        for(i=0;i<_good_points.get_dim();i++){
            if(_chifn->get_fn(_good_points.get_data(i))<target()){
                good_cache.add(_good_points.get_data(i));
            }
        }

        _reset_chimin=chimin();
        _good_points.reset_preserving_room();
        _tendril_path.reset_preserving_room();
        _exclusion_zones.reset();

        for(i=0;i<good_cache.get_dim();i++){
            _good_points.add(good_cache.get_data(i));
        }

        for(i=0;i<tendril_cache.get_rows();i++){
            _tendril_path.add_row(tendril_cache(i));
            if(origins.contains(tendril_cache.get_data(i,1))==0){
                origins.add(tendril_cache.get_data(i,1));
            }
        }

        for(i=0;i<origins.get_dim();i++){
            ellipse_pts.reset_preserving_room();
            for(j=0;j<_tendril_path.get_rows();j++){
                if(_tendril_path.get_data(j,1)==origins.get_data(i)){
                    ellipse_pts.add_row(_chifn->get_pt(_tendril_path.get_data(j,0)));
                }
            }
            if(ellipse_pts.get_rows()>2*_chifn->get_dim()){
                local_ellipse.build(ellipse_pts);
                _exclusion_zones.add(local_ellipse);
            }
        }

        find_bases();
    }
    else if(chimin()<min_00-0.01){
        find_bases();
    }

    printf("done iterating %e %d\n",chimin(),_chifn->get_called());

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
