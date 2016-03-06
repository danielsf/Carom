#include "dalex.h"

void dalex::build(chisq_wrapper *cc){

    _chifn=cc;

    int i,j;
    for(i=0;i<_chifn->get_dim();i++){
        _propagate_bisection(i);
    }

    _basis_vectors.set_cols(_chifn->get_dim());
    for(i=0;i<_chifn->get_dim();i++){
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

    if(_particles.get_dim()==0 && _origins.get_dim()==0){
        return _chifn->get_characteristic_length(dex);
    }

    double min,max;
    int i,ip;
    for(i=0;i<_particles.get_dim()+_origins.get_dim();i++){

        if(i>=_particles.get_dim()){
            ip=_origins.get_data(i-_particles.get_dim());
        }
        else{
            ip=_particles.get_data(i);
        }

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
    for(i=0;i<_chifn->get_dim();i++){
        propagate(i);
    }

    if(mindex()!=_simplex_mindex){
        simplex_search();
        if(chimin()<_basis_chimin-(target()-chimin())){
            _explorers.reset();
            find_bases();
        }
    }

    explore();
    simplex_boundary_search();

}


void dalex::simplex_search(){
    printf("\n    doing dalex_simplex %e -- %e %e\n",chimin(),target(),_target_factor);
    safety_check("simplex_search");
    array_1d<double> min,max;
    min.set_name("dalex_simplex_min");
    max.set_name("dalex_simplex_max");
    array_2d<double> seed;
    seed.set_name("dalex_simplex_seed");
    array_1d<double> trial,grad;

    int i,j;

    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,2.0*exception_value);
        max.set(i,-2.0*exception_value);
    }

    for(i=0;i<_particles.get_dim();i++){
        for(j=0;j<_chifn->get_dim();j++){
            if(_chifn->get_pt(_particles.get_data(i),j)<min.get_data(j)){
                min.set(j,_chifn->get_pt(_particles.get_data(i),j));
            }
            if(_chifn->get_pt(_particles.get_data(i),j)>max.get_data(j)){
                max.set(j,_chifn->get_pt(_particles.get_data(i),j));
            }
        }
    }

    for(i=0;i<_chifn->get_dim();i++){
        if(max.get_data(i)-min.get_data(i)<1.0e-20){
            min.set(i,0.0);
            max.set(i,_chifn->get_characteristic_length(i));
        }
    }

    int i_found;
    double wgt;
    double mu_best,mu;
    array_1d<double> trial_best;
    for(i=0;i<_particles.get_dim();i++){
        mu_best=2.0*exception_value;
        for(wgt=0.25;wgt<1.0;wgt+=0.25){
            for(j=0;j<_chifn->get_dim();j++){
                trial.set(j,wgt*_chifn->get_pt(_origins.get_data(i),j)+(1.0-wgt)*_chifn->get_pt(_particles.get_data(i),j));
            }
            evaluate(trial,&mu,&j);
            if(mu<mu_best){
                mu_best=mu;
                for(j=0;j<_chifn->get_dim();j++){
                    trial_best.set(j,trial.get_data(j));
                }
            }
        }
        seed.add_row(trial_best);
    }

    calculate_gradient(mindex(), grad);
    for(i=0;i<_chifn->get_dim();i++){
        grad.multiply_val(i,-1.0);
    }
    double gnorm,local_target=target();
    int i_min=mindex();
    i_found=i_min;
    gnorm=grad.normalize();
    while(i_found==i_min){
        i_found = bisection(i_min,grad,local_target,0.1);

        if(i_found==i_min){
            printf("    need to increase gradient target %e %e %e\n",local_target,_chifn->get_fn(i_min),gnorm);
        }
        local_target*=2.0;
    }
    if(i_min==i_found){
        printf("cannot proceed with simplex; gradient did not move %e\n",grad.normalize());
        exit(1);
    }
    seed.add_row(_chifn->get_pt(i_found)[0]);

    simplex_minimizer ffmin;
    ffmin.set_chisquared(_chifn);
    ffmin.set_minmax(min,max);
    ffmin.set_dice(_chifn->get_dice());
    ffmin.use_gradient();
    ffmin.find_minimum(seed,trial);

    printf("    after dalex_simplex chimin %e\n",chimin());
    _simplex_mindex=mindex();
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


void dalex::calculate_gradient(int i_origin, array_1d<double> &grad){

    safety_check("gradient");

    grad.reset_preserving_room();

    array_1d<double> trial;
    trial.set_name("dalex_gradient_trial");

    int ii,jj,i_found,zero_dim,aborted;
    double step_factor,mu,dx,dy;

    array_1d<double> norm,max,min;
    norm.set_dim(_chifn->get_dim());
    for(ii=0;ii<_chifn->get_dim();ii++){
        max.set(ii, -2.0*exception_value);
        min.set(ii, 2.0*exception_value);
    }

    int ip,kk;
    for(ii=0;ii<_particles.get_dim();ii++){
        if(_particles.get_data(ii)>0){
            ip=_particles.get_data(ii);
            for(jj=0;jj<_chifn->get_dim();jj++){
                mu=0.0;
                for(kk=0;kk<_chifn->get_dim();kk++){
                    mu+=_chifn->get_pt(ip,kk)*_basis_vectors.get_data(jj,kk);
                }
                if(mu<min.get_data(jj)){
                    min.set(jj,mu);
                }
                if(mu>max.get_data(jj)){
                    max.set(jj,mu);
                }
            }
        }
    }

    for(ii=0;ii<_origins.get_dim();ii++){
        if(_origins.get_data(ii)>0){
            ip=_origins.get_data(ii);
            for(jj=0;jj<_chifn->get_dim();jj++){
                mu=0.0;
                for(kk=0;kk<_chifn->get_dim();kk++){
                    mu+=_chifn->get_pt(ip,kk)*_basis_vectors.get_data(jj,kk);
                }
                if(mu<min.get_data(jj)){
                    min.set(jj,mu);
                }
                if(mu>max.get_data(jj)){
                    max.set(jj,mu);
                }
            }
        }
    }

    for(ii=0;ii<_chifn->get_dim();ii++){
        if(max.get_data(ii)-min.get_data(ii)>1.0e-20 && fabs(max.get_data(ii)-min.get_data(ii))<exception_value){
            norm.set(ii,max.get_data(ii)-min.get_data(ii));
        }
        else{
            norm.set(ii,1.0);
        }
    }

    zero_dim=0;
    aborted=0;

    array_1d<double> grad_basis;
    grad_basis.set_name("dalex_calc_grad_grad_basis");

    for(ii=0;ii<_chifn->get_dim();ii++){
        step_factor=0.001;
        while(grad_basis.get_dim()<=ii){

            for(kk=0.0;kk<_chifn->get_dim();kk++){
                trial.set(kk,_chifn->get_pt(i_origin,kk));
                trial.add_val(kk,step_factor*norm.get_data(kk)*_basis_vectors.get_data(ii,kk));
            }

            evaluate(trial,&mu,&i_found);

            if(i_found<0 || i_found==i_origin){
                for(kk=0.0;kk<_chifn->get_dim();kk++){
                    trial.set(kk,_chifn->get_pt(i_origin,kk));
                    trial.subtract_val(kk,step_factor*norm.get_data(kk)*_basis_vectors.get_data(ii,kk));
                }
                evaluate(trial,&mu,&i_found);
            }

            if(i_found>=0 && i_found!=i_origin){
                dx=0.0;
                for(kk=0;kk<_chifn->get_dim();kk++){
                    dx+=(_chifn->get_pt(i_origin,kk)-trial.get_data(kk))*_basis_vectors.get_data(ii,kk);
                }
                dy=_chifn->get_fn(i_origin)-mu;
                grad_basis.set(ii,dy/dx);
            }
            else if(i_found==i_origin){
                step_factor*=1.8;
                aborted++;
            }
            else{
                step_factor*=0.5;
                aborted++;
            }

            if(aborted>=15){
                zero_dim++;
                grad_basis.set(ii,0.0);
            }

        }

    }

    for(ii=0;ii<_chifn->get_dim();ii++){
        grad.set(ii,0.0);
    }

    for(ii=0;ii<_chifn->get_dim();ii++){
        for(kk=0;kk<_chifn->get_dim();kk++){
            grad.add_val(kk,grad_basis.get_data(ii)*_basis_vectors.get_data(ii,kk));
        }
    }

    if(zero_dim==_chifn->get_dim()){
        printf("WARNING dalex gradient zeroed out all dimensions\n");
        exit(1);
    }

}


void dalex::propagate(int dex){
    safety_check("propagate");

    if(dex>=_particles.get_dim()){
        _propagate_bisection(dex);
    }

    int i_particle=_particles.get_data(dex);
    int i_origin=_origins.get_data(dex);

    if(i_particle<=0 || i_origin<=0 || i_particle==i_origin){
        _propagate_bisection(dex);
    }
    else{
        //_propagate_midpt(dex);
        _propagate_ricochet(dex);
    }
}


void dalex::_propagate_bisection(int dex){
    safety_check("_propagate_bisection");
    /*printf("   bisecting %d %e %d\n",dex,chimin(),_chifn->get_pts());
    if(dex<_particles.get_dim()){
        printf("    %d %d %e %e\n",
        _particles.get_data(dex),_origins.get_data(dex),
        _chifn->get_fn(_particles.get_data(dex)),
        _chifn->get_fn(_origins.get_data(dex)));
    }*/
    array_1d<double> dir;
    dir.set_name("dalex_propagate_bisection_dir");
    int i,i_found_1,i_found_2;
    i_found_1=-1;

    while(i_found_1<0){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dir.normalize();
        i_found_1=bisection(mindex(),dir,target(),0.1);
        if(_particles.contains(i_found_1)==1 || _origins.contains(i_found_1)==1){
            i_found_1=-1;
        }
    }


    i_found_2=-1;
    while(i_found_2<0){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dir.normalize();
        i_found_2=bisection(mindex(),dir,target(),0.1);
        if(_origins.contains(i_found_2)==1 || _particles.contains(i_found_2)==1 || i_found_2==i_found_1){
            i_found_2=-1;
        }
    }

    double dd,dd1,dd2;
    dd1=2.0*exception_value;
    dd2=2.0*exception_value;
    for(i=0;i<_particle_log.get_dim();i++){
        dd=_chifn->distance(i_found_1,_particle_log.get_data(i));
        if(dd<dd1){
            dd1=dd;
        }
        dd=_chifn->distance(i_found_2,_particle_log.get_data(i));
        if(dd<dd2){
            dd2=dd;
        }
    }

    if(dd1>dd2){
        _particles.set(dex,i_found_1);
        _origins.set(dex,i_found_2);
    }
    else{
        _particles.set(dex,i_found_2);
        _origins.set(dex,i_found_1);
    }

    _particle_log.add(i_found_1);
    _particle_log.add(i_found_2);

}


void dalex::_propagate_ricochet(int dex){
    safety_check("_propagate_ricochet");
    //printf("    ricocheting %d\n",dex);
    array_1d<double> dir,gradient,reflected_dir;
    dir.set_name("dalex_propagate_ricochet_dir");
    gradient.set_name("dalex_propagate_ricochet_gradient");
    reflected_dir.set_name("dalex_propagate_ricochet_relfected_dir");

    if(_chifn->get_fn(_particles.get_data(dex))>_chifn->target()){
        _propagate_bisection(dex);
        return;
    }

    int i_particle,i_origin;
    i_particle=_particles.get_data(dex);
    i_origin=_origins.get_data(dex);

    calculate_gradient(i_particle, gradient);
    gradient.normalize();
    int i;
    for(i=0;i<_chifn->get_dim();i++){
        dir.set(i,_chifn->get_pt(i_particle,i)-_chifn->get_pt(i_origin,i));
    }
    dir.normalize();

    double component=0.0;
    for(i=0;i<_chifn->get_dim();i++){
        component+=dir.get_data(i)*gradient.get_data(i);
    }

    for(i=0;i<_chifn->get_dim();i++){
        reflected_dir.set(i,dir.get_data(i)-2.0*component*gradient.get_data(i));
    }

    int new_particle;
    new_particle=bisection(i_particle,reflected_dir,target(),0.1);

    _origins.set(dex,i_particle);
    _particles.set(dex,new_particle);
    if(_log!=NULL){
        _log->add(_log_ricochet,new_particle);
    }
}


void dalex::_propagate_midpt(int dex){
    safety_check("_propagate_midpt");
    array_1d<double> midpt,dir,dir_0;
    midpt.set_name("prop_mid_mid");
    dir.set_name("prop_mid_dir");
    dir_0.set_name("prop_mid_dir_0");

    int i_particle,i_origin;
    i_particle=_particles.get_data(dex);
    i_origin=_origins.get_data(dex);

    int i;
    double mu;
    for(i=0;i<_chifn->get_dim();i++){
        midpt.set(i,0.666*_chifn->get_pt(i_particle,i)+0.333*_chifn->get_pt(i_origin,i));
        dir_0.set(i,_chifn->get_pt(i_particle,i)-_chifn->get_pt(i_origin,i));
    }
    int i_mid;
    double mu_mid;
    evaluate(midpt, &mu_mid, &i_mid);
    if(i_mid<0 || mu_mid>target() || i_particle==i_origin){
        _propagate_bisection(dex);
        return;
    }

     dir_0.normalize();
     double component=0.0;
     for(i=0;i<_chifn->get_dim();i++){
         dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
     }
     for(i=0;i<_chifn->get_dim();i++){
         component+=dir.get_data(i)*dir_0.get_data(i);
     }
     for(i=0;i<_chifn->get_dim();i++){
         dir.subtract_val(i,component*dir_0.get_data(i));
     }
     dir.normalize();
     int i_found_1,i_found_2;
     i_found_1=bisection(i_mid,dir,target(),0.1);
     if(_chifn->get_fn(i_mid)>target()){
         _propagate_bisection(dex);
         return;
     }

     for(i=0;i<_chifn->get_dim();i++){
         dir.multiply_val(i,-1.0);
     }
     i_found_2=bisection(i_mid, dir, target(), 0.1);

     if(i_found_1<0 || i_found_2<0){
         _propagate_bisection(dex);
         return;
     }

     double dd1,dd2,dd;
     dd1=2.0*exception_value;
     dd2=2.0*exception_value;
     for(i=0;i<_particle_log.get_dim();i++){
         dd=_chifn->distance(i_found_1,_particle_log.get_data(i));
         if(dd<dd1){
             dd1=dd;
         }
         dd=_chifn->distance(i_found_2,_particle_log.get_data(i));
         if(dd<dd2){
             dd2=dd;
         }
     }

     if(dd1>dd2){
         _particles.set(dex,i_found_1);
         _origins.set(dex,i_found_2);
     }
     else{
         _particles.set(dex,i_found_2);
         _origins.set(dex,i_found_1);
     }

    _particle_log.add(i_found_1);
    _particle_log.add(i_found_2);

    if(_log!=NULL){
        if(i_found_1>=0){
            _log->add(_log_ricochet, i_found_1);
        }

        if(i_found_2>=0){
            _log->add(_log_ricochet, i_found_2);
        }
    }

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
    bases_out(idim)->normalize();

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

        bases_out(ix)->normalize();

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
        bases(ix)->normalize();
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
        bases(ix)->normalize();
    }

    validate_bases(bases,"dalex_guess_bases");

    printf("validated guessed bases\n");
}



void dalex::find_bases(){

    safety_check("find_bases");

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
    while(_basis_associates.get_dim()<4*_chifn->get_dim()+n_0){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dir.normalize();
        i_pt=bisection(mindex(),dir,0.5*(target()+chimin()),0.1);

        if(fabs(_chifn->get_fn(i_pt)-0.5*(target()+chimin()))>1.0 && chimin()<500.0){
            printf("WARNING failed to get associate within tol %e %e-- %e %e\n",
            0.5*(target()+chimin()),_chifn->get_fn(i_pt),chimin(),target());
            printf("first\n");
            exit(1);
        }

        if(i_pt!=mindex() && _basis_associates.contains(i_pt)==0){
            _basis_associates.add(i_pt);
            i_pt=bisection(mindex(),dir,0.25*chimin()+0.75*target(),0.1);
            if(i_pt!=mindex() && _basis_associates.contains(i_pt)==0){
                _basis_associates.add(i_pt);

                if(fabs(_chifn->get_fn(i_pt)-(0.75*target()+0.25*chimin()))>1.0 && chimin()<500.0){
                    printf("WARNING failed to get associate within tol %e %e -- %e %e\n",
                    0.75*target()+0.25*chimin(),_chifn->get_fn(i_pt),chimin(),target());
                    printf("second\n");
                    exit(1);
                }

            }
        }
    }


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

    array_1d<int> cand_1,cand_2;
    cand_1.set_name("dalex_find_bases_cand_1");
    cand_2.set_name("dalex_find_bases_cand_2");

    if(changed_bases==1){
        for(i=0;i<_chifn->get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                dir.set(j,-1.0*_basis_vectors.get_data(i,j));
            }
            i_pt=bisection(mindex(),_basis_vectors(i)[0],target(),0.1);
            _particle_log.add(i_pt);
            if(_log!=NULL && i_pt>=0){
                _log->add(_log_compass, i_pt);
            }
            cand_1.add(i_pt);
            if(fabs(_chifn->get_fn(i_pt)-target())>0.2 && chimin()<500.0){
                printf("WARNING at end of basis %e wanted %e\n",
                _chifn->get_fn(i_pt),target());
                exit(1);
            }
            i_pt=bisection(mindex(),dir,target(),0.1);
            _particle_log.add(i_pt);
            if(_log!=NULL && i_pt>=0){
                _log->add(_log_compass, i_pt);
            }
            cand_2.add(i_pt);
            if(fabs(_chifn->get_fn(i_pt)-target())>0.2 && chimin()<500.0){
                printf("WARNING at end of basis %e wanted %e\n",
                _chifn->get_fn(i_pt),target());
                exit(1);
            }
        }
    }

    _basis_chimin=chimin();

    double dd,dd1,dd2;
    int ip1,ip2,k;
    for(i=0;i<_chifn->get_dim() && changed_bases==1;i++){
        ip1=cand_1.get_data(i);
        ip2=cand_2.get_data(i);
        dd1=2.0*exception_value;
        dd2=2.0*exception_value;
        for(j=0;j<_particle_log.get_dim();j++){
            if(_particle_log.get_data(j)!=ip1 && _particle_log.get_data(j)!=ip2){
                 dd=_chifn->distance(ip1,_particle_log.get_data(j));
                 if(dd<dd1){
                     dd1=dd;
                 }
                 dd=_chifn->distance(ip2,_particle_log.get_data(j));
                 if(dd<dd2){
                     dd2=dd;
                 }
            }
        }

        if(dd1>dd2){
            _particles.set(i,ip1);
            _origins.set(i,ip2);
        }
        else{
            _particles.set(i,ip2);
            _origins.set(i,ip1);
        }
    }


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

void dalex::simplex_boundary_search(){
    printf("\ndoing dalex.simplex_boundary_search()\n");
    int pt_start=_chifn->get_pts();

    int i_node,i_pt;
    int i,j;
    double xmin,xmax,xx;

    assess_good_points();

    gp_lin interpolator;
    interpolator.set_kd_fn(_chifn->get_tree(), _chifn->get_fn_arr());
    interpolator.set_ell_factor(1.0);

    dchi_boundary_simplex_gp dchifn(_chifn,&interpolator,_good_points);

    simplex_minimizer ffmin;
    ffmin.set_chisquared(&dchifn);
    ffmin.set_dice(_chifn->get_dice());
    array_1d<double> min,max;
    min.set_name("dalex_simplex_search_min");
    max.set_name("dalex_simplex_search_min");

    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,0.0);
        max.set(i,get_norm(i));
    }

    ffmin.set_minmax(min,max);
    ffmin.use_gradient();

    array_2d<double> seed;
    seed.set_name("dalex_simplex_search_seed");

    seed.set_cols(_chifn->get_dim());
    int iFound;
    array_1d<double> trial;
    array_1d<int> seed_dex;
    double ftrial;
    trial.set_name("dalex_simplex_search_trial");
    seed_dex.set_name("dalex_simplex_search_seed_dex");
    int i_min=-1;
    double mu_min;

    for(i=0;i<_chifn->get_dim()+1;i++){
        printf("    explorer is %d\n",_explorers.get_data(i));
        seed.add_row(_chifn->get_pt(_explorers.get_data(i))[0]);
    }

    double mu,start_min;
    for(i=0;i<seed.get_rows();i++){
        mu=dchifn(seed(i)[0]);
        if(i==0 || mu<start_min){
            start_min=mu;
        }
    }
    printf("    starting from %e\n",start_min);

    array_1d<double> minpt;
    minpt.set_name("dalex_simplex_search_minpt");

    ffmin.find_minimum(seed,minpt);

    double interp_val=interpolator(minpt);

    evaluate(minpt, &mu, &i_min);

    printf("    interp %e actual %e -- %e\n",interp_val,interpolator(minpt),mu);

    if(i_min<0){
        i_min=bisection(_chifn->get_pt(_chifn->mindex())[0],minpt,target(),0.1);
        printf("    set i_min to %d\n",i_min);
    }

    if(_log!=NULL){
        _log->add(_log_dchi_simplex,i_min);
    }


    assess_good_points();

    double tol=0.01*(target()-chimin());
    int i_bisect;
    if(_chifn->get_fn(i_min)-target()>tol){
        i_bisect=bisection(_chifn->mindex(),i_min,target(),tol);
        printf("    got i_bisect %d\n",i_bisect);
        if(i_bisect>=0 && _log!=NULL){
            _log->add(_log_dchi_simplex,i_bisect);
        }
    }

    printf("    actually found %e -- %e %e\n",
    _chifn->get_fn(i_min),_chifn->get_pt(i_min,0), _chifn->get_pt(i_min,1));

    printf("    adjusted %e\n",dchifn(_chifn->get_pt(i_min)[0]));
    printf("    interpolated %e\n",interpolator(_chifn->get_pt(i_min)[0]));

    printf("    min is %e target %e\n",chimin(),target());

}


void dalex::explore(){
    printf("\nexploring\n");

    array_1d<double> trial;
    trial.set_name("dalex_explore_trial");
    int i;
    int i_found,i_high;
    double mu;
    while(_explorers.get_dim()<_chifn->get_dim()+1){
        for(i=0;i<_chifn->get_dim();i++){
            trial.set(i,_chifn->get_min(i)+
                      _chifn->random_double()*(_chifn->get_max(i)-_chifn->get_min(i)));
        }
        evaluate(trial,&mu,&i_high);
        if(mu>target()){
            i_found=bisection(_chifn->get_pt(_chifn->mindex())[0],trial,target(),0.1);
        }
        else{
            i_found=i_high;
        }
        if(i_found!=_chifn->mindex() && _explorers.contains(i_found)==0 && i_found>=0){
            _explorers.add(i_found);
        }
    }

    printf("past initialization\n");

    array_1d<double> norm,min,max;
    norm.set_name("dalex_explore_norm");
    min.set_name("dalex_explore_min");
    max.set_name("dalex_explore_max");
    int j,k;
    for(i=0;i<_explorers.get_dim();i++){
        for(j=0;j<_chifn->get_dim();j++){
            mu=0.0;
            for(k=0;k<_chifn->get_dim();k++){
                mu+=_chifn->get_pt(_explorers.get_data(i),k)*_basis_vectors.get_data(j,k);
            }
            if(i==0 || mu<min.get_data(j)){
                min.set(j,mu);
            }
            if(i==0 || mu>max.get_data(j)){
                max.set(j,mu);
            }
        }

    }

    assess_good_points();
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
        }
    }

    for(i=0;i<_chifn->get_dim();i++){
        if(max.get_data(i)-min.get_data(i)>1.0e-20){
            norm.set(i,max.get_data(i)-min.get_data(i));
        }
        else{
            norm.set(i,1.0);
        }
    }

    printf("got norm\n");

    double roll,ratio;
    int i_step,n_steps;

    n_steps=20;

    array_2d<double> pts;
    pts.set_name("dalex_explore_pts");
    for(i=0;i<_explorers.get_dim();i++){
        pts.add_row(_chifn->get_pt(_explorers.get_data(i))[0]);
    }

    printf("set points\n");

    assess_good_points();
    gp_lin interpolator;
    interpolator.set_kd_fn(_chifn->get_tree(), _chifn->get_fn_arr());
    interpolator.set_ell_factor(1.0);

    dchi_boundary_simplex_gp dchifn(_chifn, &interpolator, _good_points);

    double rr=1.0/sqrt(double(_chifn->get_dim()));

    int ip;
    array_1d<double> dir;
    dir.set_name("dalex_explore_dir");

    int accept_it;
    array_1d<int> acceptance;
    array_1d<double> ff_val;
    acceptance.set_name("dalex_explore_acceptance");
    ff_val.set_name("dalex_explore_ff_val");
    for(i=0;i<_explorers.get_dim();i++){
        acceptance.set(i,0);
        ff_val.set(i,dchifn(pts(i)[0]));
    }

    printf("taking steps\n");

    for(i_step=0;i_step<n_steps;i_step++){
        for(ip=0;ip<_explorers.get_dim();ip++){
            for(i=0;i<_chifn->get_dim();i++){
                dir.set(i,normal_deviate(_chifn->get_dice(),0.0,rr*norm.get_data(i)));
            }

            for(i=0;i<_chifn->get_dim();i++){
                trial.set(i,pts.get_data(ip,i)+dir.get_data(i));
            }
            mu=dchifn(trial);

            accept_it=0;
            if(mu<ff_val.get_data(ip)){
                 accept_it=1;
            }
            else{
                roll=_chifn->random_double();
                ratio=exp((ff_val.get_data(ip)-mu)/_explorer_temp);
                if(roll<ratio){
                   accept_it=1;
                }
            }

            if(accept_it==1){
                acceptance.add_val(ip,1);
                ff_val.set(ip,mu);
                for(i=0;i<_chifn->get_dim();i++){
                    pts.set(ip,i,trial.get_data(i));
                }
            }
        }
    }

    printf("done with steps\n");

    double mu_min,mu_max;

    for(i=0;i<_explorers.get_dim();i++){
        evaluate(pts(i)[0],&mu,&ip);
        if(i==0 || mu<mu_min){
            mu_min=mu;
        }
        if(i==0 || mu>mu_max){
            mu_max=mu;
        }
        if(ip<0){
            j=ip;
            ip=bisection(_chifn->get_pt(_chifn->mindex())[0],pts(i)[0],target()+(target()-chimin()),0.1);
        }
        _explorers.set(i,ip);
    }

    int min_acc,max_acc;
    for(i=0;i<acceptance.get_dim();i++){
        if(i==0 || acceptance.get_data(i)<min_acc){
            min_acc=acceptance.get_data(i);
        }
        if(i==0 || acceptance.get_data(i)>max_acc){
            max_acc=acceptance.get_data(i);
        }
    }

    if(min_acc<n_steps/3){
        _explorer_temp*=10.0;
    }
    else if(max_acc>(3*n_steps)/4){
        _explorer_temp*=0.15;
    }

    printf("done exploring %e %e %e\nmin %e max %e\n",
    double(min_acc)/double(n_steps),
    double(max_acc)/double(n_steps),
    _explorer_temp,mu_min,mu_max);
}
