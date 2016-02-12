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
}

int dalex::bisection(int ilow, array_1d<double> &dir, double target, double tol){
    safety_check("bisection(int, arr)");
    array_1d<double> trial_high;
    trial_high.set_name("dalex_bisection_trial_high");
    double rr=1.0;
    int ii,i_found;
    double mu=-2.0*exception_value;
    for(ii=0;ii<_chifn->get_dim();ii++){
        trial_high.set(ii,_chifn->get_pt(ilow,ii));
    }

    while(mu<=target){
        for(ii=0;ii<_chifn->get_dim();ii++){
            trial_high.add_val(ii,rr*dir.get_data(ii));
        }
        _chifn->evaluate(trial_high, &mu, &i_found);
    }

    return bisection(_chifn->get_pt(ilow)[0], trial_high, target, tol);
}


int dalex::bisection(int ilow, int ihigh, double target, double tol){
    safety_check("bisection(int, int)");
    return bisection(_chifn->get_pt(ilow)[0], _chifn->get_pt(ihigh)[0],
                     target, tol);
}


int dalex::bisection(array_1d<double>& lowball_in, array_1d<double>& highball_in,
                     double target, double tol){

    safety_check("bisection(arr, arr)");
    array_1d<double> trial,lowball,highball;
    trial.set_name("dalex_bisection_trial");
    lowball.set_name("dalex_bisection_lowball");
    highball.set_name("dalex_bisection_highball");
    int i_found,ii,jj,i_trial;

    double mu,flow,fhigh;
    _chifn->evaluate(lowball_in, &flow, &ii);
    _chifn->evaluate(highball_in, &fhigh, &jj);

    i_found=ii;

    for(ii=0;ii<_chifn->get_dim();ii++){
        lowball.set(ii,lowball_in.get_data(ii));
        highball.set(ii,highball_in.get_data(ii));
    }

    int ct;
    for(ct=0;ct<20 &&
       (ct<5 || fabs(_chifn->get_fn(i_found)-target)>tol); ct++){

        for(ii=0;ii<_chifn->get_dim();ii++){
            trial.set(ii, 0.5*(lowball.get_data(ii)+highball.get_data(ii)));
        }

        _chifn->evaluate(trial,&mu,&i_trial);

        if(mu<target){
            for(ii=0;ii<_chifn->get_dim();ii++){
                lowball.set(ii,trial.get_data(ii));
            }
        }
        else{
            for(ii=0;ii<_chifn->get_dim();ii++){
                highball.set(ii,trial.get_data(ii));
            }
        }

        if(mu<target && fabs(mu-target)<fabs(_chifn->get_fn(i_found)-target)){
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
    else if(_chifn->get_fn(i_particle)>_chifn->target() || i_particle==i_origin){
        _propagate_bisection(dex);
    }
    else{
        _propagate_ricochet(dex);
    }
}


void dalex::_propagate_bisection(int dex){
    safety_check("_propagate_bisection");
    printf("   bisecting %d %e %d\n",dex,_chifn->chimin(),_chifn->get_pts());
    array_1d<double> dir;
    dir.set_name("dalex_propagate_bisection_dir");
    int i,i_particle,i_origin;
    i_particle=-1;

    while(i_particle<0){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dir.normalize();
        i_particle=bisection(_chifn->mindex(),dir,_chifn->target(),0.1);
        if(_particles.contains(i_particle)==1){
            i_particle=-1;
        }
    }

    _particles.set(dex,i_particle);

    i_origin=-1;
    while(i_origin<0){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dir.normalize();
        i_origin=bisection(_chifn->mindex(),dir,_chifn->target(),0.1);
        if(_origins.contains(i_origin)==1 || _particles.contains(i_origin)==1){
            i_origin=-1;
        }
    }
    _origins.set(dex,i_origin);

}


void dalex::_propagate_ricochet(int dex){
    safety_check("_propagate_ricochet");
    printf("    ricocheting %d\n",dex);
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
    new_particle=bisection(i_particle,reflected_dir,_chifn->target(),0.1);

    _origins.set(dex,i_particle);
    _particles.set(dex,new_particle);
}
