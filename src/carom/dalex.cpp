#include "dalex.h"

void dalex::build(chisq_wrapper *cc){

    _chifn=cc;

    int i;
    for(i=0;i<_chifn->get_dim();i++){
        _propagate_bisection(i);
    }

}

void dalex::search(){
    safety_check("search");
    int i;
    for(i=0;i<_chifn->get_dim();i++){
        propagate(i);
    }

    if(_chifn->mindex()!=_simplex_mindex){
        simplex_search();
    }

}


void dalex::simplex_search(){
    printf("\n    doing dalex_simplex %e -- %e %e\n",_chifn->chimin(),target(),_target_factor);
    safety_check("simplex_search");
    array_1d<double> min,max;
    min.set_name("dalex_simplex_min");
    max.set_name("dalex_simplex_max");
    array_2d<double> seed;
    seed.set_name("dalex_simplex_seed");
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

    double mu_min;
    int i_min=-1;
    for(i=0;i<_particles.get_dim();i++){
        seed.add_row(_chifn->get_pt(_origins.get_data(i))[0]);
        if(i_min<0 || _chifn->get_fn(_particles.get_data(i))<mu_min){
            if(_particles.get_data(i)!=_chifn->mindex()){
                mu_min=_chifn->get_fn(_particles.get_data(i));
                i_min=_particles.get_data(i);
            }
        }
    }

    array_1d<double> trial;
    for(i=0;i<_chifn->get_dim();i++){
        trial.set(i,_chifn->get_pt(_chifn->mindex(),i)+0.5*(_chifn->get_pt(_chifn->mindex(),i)-_chifn->get_pt(i_min,i)));
    }
    seed.add_row(trial);

    simplex_minimizer ffmin;
    ffmin.set_chisquared(_chifn);
    ffmin.set_minmax(min,max);
    ffmin.set_dice(_chifn->get_dice());
    ffmin.use_gradient();
    ffmin.find_minimum(seed,trial);

    printf("    after dalex_simplex chimin %e\n",_chifn->chimin());
    _simplex_mindex=_chifn->mindex();
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
        _chifn->evaluate(trial_high, &mu, &i_found);
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
    _chifn->evaluate(lowball_in, &flow, &ii);
    _chifn->evaluate(highball_in, &fhigh, &jj);

    if(flow>local_target){
        printf("WARNING flow is greater than target %e %e\n",flow,local_target);
        exit(1);
    }

    i_found=ii;

    for(ii=0;ii<_chifn->get_dim();ii++){
        lowball.set(ii,lowball_in.get_data(ii));
        highball.set(ii,highball_in.get_data(ii));
    }

    int ct;
    for(ct=0;ct<20 &&
       (ct<5 || fabs(_chifn->get_fn(i_found)-local_target)>tol); ct++){

        for(ii=0;ii<_chifn->get_dim();ii++){
            trial.set(ii, 0.5*(lowball.get_data(ii)+highball.get_data(ii)));
        }

        _chifn->evaluate(trial,&mu,&i_trial);

        if(mu<local_target){
            for(ii=0;ii<_chifn->get_dim();ii++){
                lowball.set(ii,trial.get_data(ii));
            }
        }
        else{
            for(ii=0;ii<_chifn->get_dim();ii++){
                highball.set(ii,trial.get_data(ii));
            }
        }

        if(mu<local_target && fabs(mu-local_target)<fabs(_chifn->get_fn(i_found)-local_target)){
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
    for(jj=0;jj<_chifn->get_dim();jj++){
        trial.set(jj,_chifn->get_pt(i_origin,jj));
    }

    zero_dim=0;
    aborted=0;

    for(ii=0;ii<_chifn->get_dim();ii++){
        step_factor=0.001;
        while(grad.get_dim()<=ii){

            trial.set(ii,_chifn->get_pt(i_origin,ii));

            trial.add_val(ii,step_factor*_chifn->get_characteristic_length(ii));
            _chifn->evaluate(trial,&mu,&i_found);

            if(i_found<0 || i_found==i_origin){
                trial.set(ii,_chifn->get_pt(i_origin,ii));
                trial.subtract_val(ii,step_factor*_chifn->get_characteristic_length(ii));
                _chifn->evaluate(trial,&mu,&i_found);
            }

            if(i_found>=0 && i_found!=i_origin){
                dx=_chifn->get_pt(i_origin,ii)-trial.get_data(ii);
                dy=_chifn->get_fn(i_origin)-mu;
                grad.add(dy/dx);
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
                grad.add(0.0);
            }

        }

        trial.set(ii,_chifn->get_pt(i_origin,ii));

    }

    if(zero_dim==_chifn->get_dim()){
        printf("WARNING dalex gradient zeroed out all dimensions\n");
        exit(1);
    }

}


void dalex::propagate(int dex){
    safety_check("propagate");

    if(dex>_particles.get_dim()){
        _propagate_bisection(dex);
    }

    int i_particle=_particles.get_data(dex);
    int i_origin=_origins.get_data(dex);

    if(i_particle<=0 || i_origin<=0){
        _propagate_bisection(dex);
    }
    else if(_chifn->get_fn(i_particle)>target() || i_particle==i_origin){
        _propagate_bisection(dex);
    }
    else{
        //_propagate_ricochet(dex);
        _propagate_midpt(dex);
    }
}


void dalex::_propagate_bisection(int dex){
    safety_check("_propagate_bisection");
    /*printf("   bisecting %d %e %d\n",dex,_chifn->chimin(),_chifn->get_pts());
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
        i_found_1=bisection(_chifn->mindex(),dir,target(),0.1);
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
        i_found_2=bisection(_chifn->mindex(),dir,target(),0.1);
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
    _chifn->evaluate(midpt, &mu_mid, &i_mid);
    if(i_mid<0 || mu_mid>_chifn->target() || i_particle==i_origin){
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
     i_found_1=bisection(i_mid,dir,_chifn->target(),0.1);
     for(i=0;i<_chifn->get_dim();i++){
         dir.multiply_val(i,-1.0);
     }
     i_found_2=bisection(i_mid, dir, _chifn->target(), 0.1);

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

}
