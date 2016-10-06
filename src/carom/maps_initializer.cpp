#include "maps_initializer.h"

void maps_initializer::initialize(){
    safety_check();
    printf("starting min is %e\n",_chifn->chimin());
    _particles.reset_preserving_room();
    _fn.reset_preserving_room();
    _dex.reset_preserving_room();
    _transform.reset_preserving_room();
    _ct_replace=0;
    int n_particles=1000;
    int i,j;
    for(i=0;i<_chifn->get_dim();i++){
        _transform.set(i,_chifn->get_characteristic_length(i));
    }

    array_1d<double> center;
    for(i=0;i<_chifn->get_dim();i++){
        center.set(i,0.5*(_chifn->get_min(i)+_chifn->get_max(i)));
    }

    _fn_max=-1.0;
    _max_dex=-1;
    array_1d<double> trial,shell,dir;
    double mu;
    int i_found;
    while(_particles.get_rows()<n_particles){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dir.normalize();
        for(i=0;i<_chifn->get_dim();i++){
            trial.set(i,center.get_data(i)+
                     0.5*(_chifn->get_max(i)-_chifn->get_min(i))*dir.get_data(i));
        }
        _chifn->evaluate(trial,&mu,&i_found);
        if(i_found>=0){
            convert_to_shell(trial,shell);
            _dex.add(i_found);
            _fn.add(mu);
            _particles.add_row(shell);
            if(_max_dex<0 || mu>_fn_max){
                _fn_max=mu;
                _max_dex=_particles.get_rows()-1;
            }
        }
    }

    set_bases();
    printf("after init %e\n",_chifn->chimin());

    /*double dd;
    for(i=0;i<_chifn->get_dim();i++){
        printf("    %e\n",_chifn->get_pt(_chifn->mindex(),i));
        dd+=power(_chifn->get_pt(_chifn->mindex(),i),2);
    }
    printf("dd %e\n",sqrt(dd));
    exit(1);*/
}


void maps_initializer::convert_to_truth(const array_1d<double> &in, array_1d<double> &out){
    int i;
    for(i=0;i<_chifn->get_dim();i++){
        out.set(i,in.get_data(i)*_transform.get_data(i));
    }
}

void maps_initializer::convert_to_shell(const array_1d<double> &in, array_1d<double> &out){
    int i;
    for(i=0;i<_chifn->get_dim();i++){
        out.set(i,in.get_data(i)/_transform.get_data(i));
    }
}

void maps_initializer::set_bases(){
    _bases.reset_preserving_room();
    _radii.reset_preserving_room();

    int dim,i,j;
    double fn_min;
    for(i=0;i<_particles.get_rows();i++){
        if(i==0 || _fn.get_data(i)<fn_min){
            fn_min=_fn.get_data(i);
            for(j=0;j<_chifn->get_dim();j++){
                _center.set(j,_particles.get_data(i,j));
            }
        }
    }
    double mu, mu_best;
    array_1d<double> dir,dir_best;

    _vol=1.0;
    while(_bases.get_rows()<_chifn->get_dim()){
        for(i=0;i<_particles.get_rows();i++){
            for(dim=0;dim<_chifn->get_dim();dim++){
                dir.set(dim,_particles.get_data(i,dim)-_center.get_data(dim));
            }

            for(j=0;j<_bases.get_rows();j++){
                mu=0.0;
                for(dim=0;dim<_chifn->get_dim();dim++){
                    mu+=dir.get_data(dim)*_bases.get_data(j,dim);
                }
                for(dim=0;dim<_chifn->get_dim();dim++){
                    dir.subtract_val(dim,mu*_bases.get_data(j,dim));
                }
            }

            mu=dir.normalize();
            if(i==0 || mu>mu_best){
                mu_best=mu;
                for(dim=0;dim<_chifn->get_dim();dim++){
                    dir_best.set(dim,dir.get_data(dim));
                }
            }
        }
        for(i=0;i<_bases.get_rows();i++){
            mu=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                mu+=_bases.get_data(i,j)*dir_best.get_data(j);
            }
            if(fabs(mu)>1.0e-4){
                printf("WARNING basis dot prod %e\n",mu);
                exit(1);
            }
        }
        mu=dir_best.normalize();
        if(fabs(mu-1.0)>1.0e-5){
            printf("WARNING basis norm %e\n",mu);
            exit(1);
        }
        _bases.add_row(dir_best);
        _radii.add(mu_best);
        _vol*=mu_best;
    }

}

void maps_initializer::evaluate(array_1d<double> &pt, double *mu, int *i_found){
    convert_to_truth(pt, _pt_true);
    _chifn->evaluate(_pt_true, mu, i_found);
}

void maps_initializer::sample(){

    int i,j;
    double rr;
    rr=-1.0;
    while(rr<0.0){
        rr=normal_deviate(_chifn->get_dice(),1.0,1.0);
    }
    for(i=0;i<_chifn->get_dim();i++){
        _dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
    }
    _dir.normalize();
    for(i=0;i<_chifn->get_dim();i++){
        _pt_shell.set(i,_center.get_data(i));
    }
    for(i=0;i<_chifn->get_dim();i++){
        for(j=0;j<_chifn->get_dim();j++){
            _pt_shell.add_val(j,rr*_radii.get_data(i)*_dir.get_data(i)*_bases.get_data(i,j));
        }
    }

    double mu;
    int i_found;
    evaluate(_pt_shell, &mu, &i_found);

    if(mu<_fn_max){
        _ct_replace++;
        _particles.set_row(_max_dex,_pt_shell);
        _fn.set(_max_dex,mu);
        _dex.set(_max_dex,i_found);
        for(i=0;i<_particles.get_rows();i++){
            if(i==0 || _fn.get_data(i)>_fn_max){
                _fn_max=_fn.get_data(i);
                _max_dex=i;
            }
        }
    }



}

void maps_initializer::search(){
    safety_check();
    initialize();
    int n_steps=200000;
    int i;
    for(i=0;i<n_steps && _vol>1.0;i++){
        sample();
        if(i%(100*_chifn->get_dim())==0 && i>0){
            printf("called %d min %e -- %d %e -- %e\n",
            _chifn->get_called(),_chifn->chimin(),_ct_replace,_fn_max,_vol);
            set_bases();
        }
    }

    array_2d<double> seed;
    array_1d<double> min,max;
    array_1d<double> fn_sorted;
    array_1d<int> fn_dex;
    for(i=0;i<_fn.get_dim();i++){
        fn_dex.add(i);
    }
    sort(_fn,fn_sorted,fn_dex);
    array_1d<double> trial;
    for(i=0;i<_chifn->get_dim()+1;i++){
        convert_to_truth(_particles(i),trial);
        seed.add_row(trial);
    }
    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,0.0);
        max.set(i,_chifn->get_characteristic_length(i));
    }

    simplex_minimizer ffmin;
    ffmin.set_chisquared(_chifn);
    ffmin.set_minmax(min,max);
    ffmin.set_dice(_chifn->get_dice());
    ffmin.use_gradient();
    ffmin.find_minimum(seed,trial);

}
