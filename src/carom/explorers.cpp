#include "explorers.h"

void explorers::set_bases(){

    array_1d<double> vv;
    vv.set_name("exp_set_bases_vv");
    int i,j;
    double component;

    _bases.reset_preserving_room();
    while(_bases.get_rows()<_chifn->get_dim()){
        for(i=0;i<_chifn->get_dim();i++){
            vv.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        for(i=0;i<_bases.get_rows();i++){
            component=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                component+=vv.get_data(j)*_bases.get_data(i,j);
            }
            for(j=0;j<_chifn->get_dim();j++){
                vv.subtract_val(j,component*_bases.get_data(i,j));
            }
        }
        component=vv.normalize();
        if(component>1.0e-10){
            _bases.add_row(vv);
        }
    }

    array_1d<double> min,max;
    min.set_name("exp_set_base_min");
    max.set_name("exp_set_base_max");

    int k;
    for(i=0;i<_particles.get_rows();i++){
        for(j=0;j<_chifn->get_dim();j++){
            component=0.0;
            for(k=0;k<_chifn->get_dim();k++){
                component+=_particles.get_data(i,k)*_bases.get_data(j,k);
            }
            if(j>=min.get_dim() || component<min.get_data(j)){
                min.set(j,component);
            }
            if(j>=max.get_dim() || component>max.get_data(j)){
                max.set(j,component);
            }
        }
    }

    for(i=0;i<_associates.get_dim();i++){
        for(j=0;j<_chifn->get_dim();j++){
            component=0.0;
            for(k=0;k<_chifn->get_dim();k++){
                component+=_chifn->get_pt(_associates.get_data(i),k)*_bases.get_data(j,k);
            }
            if(j>=min.get_dim() || component<min.get_data(j)){
                min.set(j,component);
            }
            if(j>=max.get_dim() || component>max.get_data(j)){
                max.set(j,component);
            }
        }
    }

    for(i=0;i<_chifn->get_dim();i++){
        _norm.set(i,0.1*(max.get_data(i)-min.get_data(i)));
    }

}

void explorers::initialize_particles(){

    _particles.reset();

    array_1d<double> dir,trial;
    dir.set_name("exp_init_dir");
    trial.set_name("exp_init_trial");
    int i,j;

    double mu;
    int i_found;

    while(_particles.get_rows()<_n_particles){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dir.normalize();
        for(i=0;i<_chifn->get_dim();i++){
            trial.set(i,_chifn->get_pt(_chifn->mindex(),i));
        }
        for(i=0;i<_chifn->get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                trial.add_val(j,10.0*_norm.get_data(i)*dir.get_data(i)*_bases.get_data(i,j));
            }
        }
        _chifn->evaluate(trial,&mu,&i_found);
        if(i_found>=0){
            _particles.add_row(trial);
        }

    }

}


void explorers::sample(int n_steps, array_2d<double> &model_bases){

    set_bases();

    initialize_particles();

    dchi_interior_simplex dchifn(_chifn, _associates, model_bases);

    _mindex=-1;

    double mu;
    array_1d<double> mu_arr;
    mu_arr.set_name("exp_sample_mu_arr");

    int i_step;
    int i_dim,ip;
    int i;
    int accept_it;
    double rr,roll,ratio;
    array_1d<double> trial;
    trial.set_name("exp_sample_trial");

    if(_accepted.get_dim()!=_n_particles){
       _accepted.set_dim(_n_particles);
       _accepted.zero();
    }

    for(ip=0;ip<_n_particles;ip++){
        mu=dchifn(_particles(ip)[0]);
        if(ip==0 || mu<_mu_min){
            _mindex=ip;
            _mu_min=mu;
        }
        mu_arr.set(ip,mu);
    }

    double needed_temp;
    double old_temp;
    array_1d<int> req_temp_dex;
    array_1d<double> req_temp_sorted;
    req_temp_dex.set_name("exp_sample_req_temp_dex");
    req_temp_sorted.set_name("exp_sample_req_temp_sorted");

    int has_been_adjusted;

    array_1d<double> acceptance_rate,acceptance_rate_sorted;
    array_1d<int> acceptance_rate_dex;
    acceptance_rate.set_name("exp_sample_acceptance_rate");
    acceptance_rate_sorted.set_name("exp_sample_acceptance_rate_sorted");
    acceptance_rate_dex.set_name("exp_sample_acceptance_rate_dex");

    double med_acc;

    printf("    starting sampling with %e\n",_mu_min);

    int scalar_acceptance=0;

    for(i_step=0;i_step<n_steps;i_step++){
        if(i_step>0 && i_step%(_chifn->get_dim())==0){
            set_bases();
        }

        for(ip=0;ip<_n_particles;ip++){
            i_dim=_chifn->random_int()%_chifn->get_dim();
            rr=normal_deviate(_chifn->get_dice(),0.0,1.0);
            for(i=0;i<_chifn->get_dim();i++){
                trial.set(i,_particles.get_data(ip,i)+
                          rr*_norm.get_data(i_dim)*_bases.get_data(i_dim,i));
            }
            accept_it=0;
            mu=dchifn(trial);
            if(mu<mu_arr.get_data(ip)){
                accept_it=1;
                needed_temp=1.0;
            }
            else{
                ratio=exp(-0.5*(mu-mu_arr.get_data(ip))/_temp);
                roll=_chifn->random_double();
                needed_temp=-0.5*(mu-mu_arr.get_data(ip))/log(roll);
                if(roll<ratio){
                    accept_it=1;
                }
            }

            _req_temp.add(needed_temp);

            if(accept_it==1){
                scalar_acceptance++;
                _accepted.add_val(ip,1);
                mu_arr.set(ip,mu);
                for(i=0;i<_chifn->get_dim();i++){
                    _particles.set(ip,i,trial.get_data(i));
                }
                if(mu<_mu_min){
                     _mindex=ip;
                     _mu_min=mu;
                }
            }

        }
        _attempted++;

        if(_attempted>0 && _attempted%(2*_chifn->get_dim())==0){
            acceptance_rate_dex.reset_preserving_room();
            for(i=0;i<_accepted.get_dim();i++){
                acceptance_rate.set(i,double(_accepted.get_data(i))/double(_attempted));
                acceptance_rate_dex.set(i,i);
            }
            sort_and_check(acceptance_rate, acceptance_rate_sorted, acceptance_rate_dex);
            med_acc=acceptance_rate_sorted.get_data(acceptance_rate_dex.get_dim()/2);
            if(med_acc>0.75 || med_acc<0.3333){
                old_temp=_temp;
                req_temp_sorted.reset_preserving_room();
                for(i=0;i<_req_temp.get_dim();i++){
                    req_temp_dex.set(i,i);
                }
                sort_and_check(_req_temp, req_temp_sorted, req_temp_dex);
                _temp=req_temp_sorted.get_data(req_temp_dex.get_dim()/2);
                if(fabs(1.0-old_temp/_temp)>0.01){
                    _req_temp.reset();
                    _accepted.zero();
                    _attempted=0;
                }
            }
        }

    }

    printf("    sampling min %e\n",_mu_min);
    printf("    temp %e\n",_temp);
    printf("    accepted %d steps %d\n",scalar_acceptance,
    _n_particles*n_steps);
}
