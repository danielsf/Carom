#include "dalex.h"

void dalex::build(chisq_wrapper *cc){

    _chifn=cc;

    _explorers.set_chifn(_chifn);

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
    int i;
    int pts_0=_chifn->get_pts();
    assess_good_points();

    int has_explored=0;

    if(mindex()!=_simplex_mindex){
        find_bases();
        explore();
        simplex_search();
        has_explored=1;
        if(chimin()<_basis_chimin-(target()-chimin())){
            find_bases();
        }
    }

    if(has_explored==0){
        explore();
    }
    tendril_search();
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
        seed.add_row(_chifn->get_pt(specified.get_data(i))[0]);
    }

    while(seed.get_rows()<_chifn->get_dim()+1){
        i=_chifn->random_int()%_explorers.get_n_particles();
        if(chosen_seed.contains(i)==0){
            _explorers.get_pt(i,trial);
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

    if(_basis_chimin-chimin()>target()-chimin()){
        find_bases();
    }

}


int dalex::bisection(int ilow, array_1d<double> &dir, double local_target, double tol){
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

    return bisection(_chifn->get_pt(ilow)[0], trial_high, local_target, tol);
}


int dalex::bisection(int ilow, int ihigh, double local_target, double tol){
    safety_check("bisection(int, int)");
    return bisection(_chifn->get_pt(ilow)[0], _chifn->get_pt(ihigh)[0],
                     local_target, tol);
}


int dalex::bisection(array_1d<double>& lowball_in, array_1d<double>& highball_in,
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

void dalex::find_bases(){

    safety_check("find_bases");

    assess_good_points();

    array_1d<double> dir;
    int i,j;
    int i_pt;
    double mu;
    int ct;
    for(ct=0;ct<4*_chifn->get_dim();ct++){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dir.normalize();
        i_pt=bisection(mindex(),dir,0.5*(target()+chimin()),0.1);
        add_charge(i_pt);

    }

    array_1d<int> associates;
    associates.set_name("find_bases_associates");

    int thin_by=0;
    if(_good_points.get_dim()>20000){
        thin_by=_good_points.get_dim()/20000;
    }
    if(thin_by==1){
        thin_by=2;
    }

    ct=0;
    for(i=0;i<_good_points.get_dim();i++){
        if(thin_by==0 || i%ct==0){
            associates.add(_good_points.get_data(i));
        }
    }

    dchi_interior_simplex dchifn(_chifn,associates);
    dchifn.copy_bases(_basis_vectors);
    for(i=0;i<_chifn->get_dim();i++){
        _basis_norm.set(i,dchifn.get_hyper_norm(i));
    }
    _basis_chimin=chimin();
}

int dalex::simplex_boundary_search(){
    return simplex_boundary_search(-1,0);
}

int dalex::simplex_boundary_search(int specified, int use_median){

    safety_check("simplex_boundary_search");
    printf("\ndoing dalex.simplex_boundary_search() %d\n",_chifn->get_pts());
    int pt_start=_chifn->get_pts();
    assess_good_points();
    assess_charges();

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
            if(specified<=0){
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

    double v0,pv0;
    array_1d<double> pvmin,pvmax,vmin,vmax;
    vmin.set_name("dalex_bou_vmin");
    vmax.set_name("dalex_bou_vmax");
    pvmin.set_name("dalex_bou_pvmin");
    pvmax.set_name("dalex_bou_pvmax");

    v0=-1.0;
    pv0=-1.0;

    for(i=0;i<associates.get_dim();i++){
        for(j=0;j<_chifn->get_dim();j++){
            mu=_chifn->get_pt(associates.get_data(i),j);

            if(i==0 || mu<vmin.get_data(j)){
                vmin.set(j,mu);
            }

            if(i==0 || mu>vmax.get_data(j)){
                vmax.set(j,mu);
            }

            mu=0.0;
            for(k=0;k<_chifn->get_dim();k++){
                mu+=_chifn->get_pt(associates.get_data(i),k)*_basis_vectors.get_data(j,k);
            }

            if(i==0 || mu<pvmin.get_data(j)){
                pvmin.set(j,mu);
            }

            if(i==0 || mu>pvmax.get_data(j)){
                pvmax.set(j,mu);
            }
        }

    }

    if(pvmin.get_dim()==_chifn->get_dim()){
        v0=1.0;
        pv0=1.0;
        for(i=0;i<_chifn->get_dim();i++){
            v0*=(vmax.get_data(i)-vmin.get_data(i));
            pv0*=(pvmax.get_data(i)-pvmin.get_data(i));
        }
    }

    printf("    associates %d path %d\n", associates.get_dim(),_tendril_path.get_rows());

    dchi_interior_simplex dchifn(_chifn,associates);

    if(use_median==1){
        dchifn.use_median();
    }

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

    array_2d<double> seed;
    seed.set_name("dalex_simplex_search_seed");

    array_2d<double> dummy_bases;

    if(specified>=0){
        dchifn.copy_bases(dummy_bases);
        seed.add_row(_chifn->get_pt(specified)[0]);
        for(i=0;i<_chifn->get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                trial.set(j,seed.get_data(0,j)+dummy_bases.get_data(i,j)*dchifn.get_hyper_norm(i)*0.01);
            }
            seed.add_row(trial);
        }
    }
    else{
        _explorers.get_seed(seed);
    }

    int i_min=-1;
    double mu_min;
    double start_min;
    int i_start_min;
    for(i=0;i<seed.get_rows();i++){
        mu=dchifn(seed(i)[0]);
        if(i==0 || mu<start_min){
            start_min=mu;
            _chifn->evaluate(seed(i)[0],&mu,&i_start_min);
        }
    }
    printf("    starting from %e\n",start_min);

    array_1d<double> minpt;
    minpt.set_name("dalex_simplex_search_minpt");

    ffmin.find_minimum(seed,minpt);

    for(i=specified;i<_chifn->get_pts();i++){
        if(i>=0 && _chifn->get_fn(i)<target()){
            i_min=i;
        }
    }

    array_1d<int> path_row;
    path_row.set_name("path_row");

    int start_path;
    if(pt_start<i_start_min){
        start_path=i_start_min;
    }
    else{
        start_path=pt_start;
    }

    for(i=start_path;i<_chifn->get_pts();i++){
        if(_tendril_path.get_rows()==0){
            _tendril_path.set_cols(2);
        }
        if(_chifn->get_fn(i)<_chifn->target()){
            path_row.set(0,i);
            if(distance(i,i_min)<distance(i,i_start_min)){
                path_row.set(1,i_min);
            }
            else{
                path_row.set(1,i_start_min);
            }
            _tendril_path.add_row(path_row);
        }

    }


    int i_good_start;

    _update_good_points(pt_start);

    if(_log!=NULL){
        _log->add(_log_dchi_simplex,i_min);
    }

    if(_chifn->get_dim()>9){
        printf("    actually found %e -- %e %e\n",
        _chifn->get_fn(i_min),_chifn->get_pt(i_min,6), _chifn->get_pt(i_min,9));
    }

    printf("    adjusted %e from %e\n",
    dchifn(_chifn->get_pt(i_min)[0]),_chifn->get_fn(i_min));

    printf("    min is %e target %e\n",chimin(),target());
    if(_chifn->get_dim()>9){
       printf("    minpt at %e %e\n",
       _chifn->get_pt(mindex(),6),
       _chifn->get_pt(mindex(),9));
    }

    printf("    v0 %e pv0 %e\n",v0,pv0);

    double v1,pv1;
    if(v0<0.0){
        return 0;
    }

    for(i=0;i<_chifn->get_dim();i++){
        mu=_chifn->get_pt(i_min,i);
        if(mu<vmin.get_data(i)){
            vmin.set(i,mu);
        }
        if(mu>vmax.get_data(i)){
            vmax.set(i,mu);
        }
        mu=0.0;
        for(j=0;j<_chifn->get_dim();j++){
            mu+=_chifn->get_pt(i_min,j)*_basis_vectors.get_data(i,j);
        }
        if(mu<pvmin.get_data(i)){
            pvmin.set(i,mu);
        }
        if(mu>pvmax.get_data(i)){
            pvmax.set(i,mu);
        }
    }

    v1=1.0;
    pv1=1.0;
    for(i=0;i<_chifn->get_dim();i++){
        v1*=(vmax.get_data(i)-vmin.get_data(i));
        pv1*=(pvmax.get_data(i)-pvmin.get_data(i));
    }

    printf("    v1 %e pv1 %e\n",v1,pv1);

    if(v1>1.1*v0 || pv1>1.1*pv0){
        return 0;
    }

    return 1;

}


void dalex::explore(){
    printf("\nexploring\n");
    int pt_0=_chifn->get_pts();

    _explorers.set_n_particles(2*_chifn->get_dim());
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

    _explorers.set_associates(associates);
    _explorers.sample(4*_chifn->get_dim());

    if(_log!=NULL){
         for(i=pt_0;i<_chifn->get_pts();i++){
             _log->add(_log_mcmc, i);
         }
    }

    _update_good_points(pt_0);
}

void dalex::get_gradient(int origin, array_1d<double> &norm, array_1d<double> &grad_out){

    array_1d<double> trial,grad;
    trial.set_name("dalex_gradient_trial");
    grad.set_name("dalex_gradient_grad");
    int i;
    for(i=0;i<_chifn->get_dim();i++){
        trial.set(i,_chifn->get_pt(origin,i));
    }
    double mu_origin;
    evaluate(trial,&mu_origin,&i);
    int ix;
    double dx,y1,y2;
    for(ix=0;ix<_chifn->get_dim();ix++){
        dx=0.001*norm.get_data(ix);
        for(i=0;i<_chifn->get_dim();i++){
            trial.set(i,_chifn->get_pt(origin,i)-dx*_basis_vectors.get_data(ix,i));
        }
        evaluate(trial,&y1,&i);
        for(i=0;i<_chifn->get_dim();i++){
            trial.set(i,_chifn->get_pt(origin,i)+dx*_basis_vectors.get_data(ix,i));
        }
        evaluate(trial,&y2,&i);
        grad.set(ix,(y2-y1)/(2.0*dx));
        grad_out.set(ix,0.0);
    }

    for(ix=0;ix<_chifn->get_dim();ix++){
        for(i=0;i<_chifn->get_dim();i++){
            grad_out.add_val(i,grad.get_data(ix)*_basis_vectors.get_data(ix,i));
        }
    }
}


void dalex::tendril_search(){

    add_charge(_chifn->mindex());

    int i,j,k;
    int pt_0=_chifn->get_pts();
    assess_good_points();
    _update_good_points();
    int n_good_0=_good_points.get_dim();

    /*array_1d<double> dir,trial;
    array_1d<double> c_min,c_max;
    double mu;
    for(i=0;i<_chifn->get_dim();i++){
        c_min.set(i,2.0*exception_value);
        c_max.set(i,-2.0*exception_value);
    }

    for(i=0;i<_end_points.get_dim();i++){
        for(j=0;j<_chifn->get_dim();j++){
            mu=0.0;
            for(k=0;k<_chifn->get_dim();k++){
                mu+=_chifn->get_pt(_end_points.get_data(i),k)*_basis_vectors.get_data(j,k);
            }
            if(mu<c_min.get_data(j)){
                c_min.set(j,mu);
            }
            if(mu>c_max.get_data(j)){
                c_max.set(j,mu);
            }
        }
    }

    for(i=0;i<_chifn->get_dim();i++){
        if(c_min.get_data(i)>exception_value || c_max.get_data(i)<-1.0*exception_value ||
           c_max.get_data(i)-c_min.get_data(i)<1.0e-10){

            c_min.set(i,0.0);
            c_max.set(i,_basis_norm.get_data(i));

        }
    }

    int n_cand=100;
    int i_cand;
    int i_found;
    int i_particle=-1;
    double dd,dd_best,dd_local_min;
    dd_best=-1.0;
    for(i_cand=0;i_cand<n_cand || i_particle<0;i_cand++){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dir.normalize();
        for(i=0;i<_chifn->get_dim();i++){
            trial.set(i,_chifn->get_pt(0.0,i));
        }
        for(i=0;i<_chifn->get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                trial.add_val(j,0.5*dir.get_data(i)*(c_max.get_data(i)-c_min.get_data(i))*_basis_vectors.get_data(i,j));
            }
        }
        i_found=bisection(mindex(),trial,target(),0.1);
        if(i_found>=0){
            dd_local_min=distance(i_found,mindex());
            for(i=0;i<_end_points.get_dim();i++){
                dd=distance(i_found,_end_points.get_data(i));
                if(dd<dd_local_min){
                    dd_local_min=dd;
                }
            }
            if(dd_local_min>dd_best){
                i_particle=i_found;
                dd_best=dd_local_min;
            }
        }
    }*/


    double mu;
    int i_found;
    array_1d<int> specified;
    array_1d<double> trial;
    specified.add(mindex());
    for(j=0;j<_chifn->get_dim();j++){
        mu=2.0*(_chifn->random_double()-0.5);
        for(i=0;i<_chifn->get_dim();i++){
            trial.set(i,_chifn->get_pt(mindex(),i)+mu*_basis_vectors.get_data(j,i)*_basis_norm.get_data(j));
        }
        _chifn->evaluate(trial,&mu,&i_found);
        specified.add(i_found);
    }
    simplex_search(specified);

    simplex_boundary_search();
    _update_good_points();

    int i_particle=_good_points.get_data(_good_points.get_dim()-1);

    if(_log!=NULL){
        _log->add(_log_dchi_simplex,i_particle);
    }

    double volume,p_volume;
    array_1d<double> min,max,min_p,max_p;
    min.set_name("dalex_tendril_min");
    max.set_name("dalex_tendril_max");
    min_p.set_name("dalex_tendril_min_p");
    max_p.set_name("dalex_tendril_max_p");

    assess_good_points();
    int ip,ix;
    ip=mindex();
    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,_chifn->get_pt(ip,i));
        max.set(i,_chifn->get_pt(ip,i));

        mu=0.0;
        for(j=0;j<_chifn->get_dim();j++){
            mu+=_chifn->get_pt(ip,j)*_basis_vectors.get_data(i,j);
        }

        min_p.set(i,mu);
        max_p.set(i,mu);
    }

    double volume_0,p_volume_0;
    volume_0=1.0;
    p_volume_0=1.0;
    for(i=0;i<_chifn->get_dim();i++){
        volume_0*=(max.get_data(i)-min.get_data(i));
        p_volume_0*=(max_p.get_data(i)-min_p.get_data(i));
    }

    printf("    volume %e %e\n",volume_0,p_volume_0);

    array_1d<double> dir1,dir2,trial_center;
    dir1.set_name("dalex_simplex_boundary_dir1");
    dir2.set_name("dalex_simplex_boundary_dir2");
    trial_center.set_name("dalex_simplex_boundary_trial_center");

    int i_origin,ct_last;

    int strikes=0;
    int iteration=0;
    int use_median=0;
    int is_a_strike;

    while(strikes<3){

        /*if(strikes>0){
            use_median=1;
        }
        else{
            use_median=0;
        }*/

        iteration++;

        printf("    strikes %d use_median %d\n",strikes,use_median);
        add_charge(_chifn->mindex());

        i_origin=i_particle;
        ct_last=_chifn->get_pts();
        is_a_strike=simplex_boundary_search(i_particle, use_median);

        i_particle=_good_points.get_data(_good_points.get_dim()-1);

        add_charge(i_particle);

        if(is_a_strike==1){
            strikes++;
            i_particle=i_origin;
        }
        else{
            strikes=0;
        }

        for(i=0;i<_chifn->get_dim();i++){
            if(_chifn->get_pt(i_particle,i)<min.get_data(i)){
                min.set(i,_chifn->get_pt(i_particle,i));
            }

            if(_chifn->get_pt(i_particle,i)>max.get_data(i)){
                max.set(i,_chifn->get_pt(i_particle,i));
            }

            mu=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                mu+=_chifn->get_pt(i_particle,j)*_basis_vectors.get_data(i,j);
            }

            if(mu<min_p.get_data(i)){
                min_p.set(i,mu);
            }

            if(mu>max_p.get_data(i)){
                max_p.set(i,mu);
            }
        }

        volume=1.0;
        p_volume=1.0;
        for(i=0;i<_chifn->get_dim();i++){
            volume*=(max.get_data(i)-min.get_data(i));
            p_volume*=(max_p.get_data(i)-min_p.get_data(i));
        }

        if(volume>volume_0*1.1 && p_volume>p_volume_0*1.1){
            volume_0=volume;
            p_volume_0=p_volume;
        }
        printf("    volume %e %e\n",volume_0,p_volume_0);

    }

}

