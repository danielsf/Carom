#include "maps.h"

maps::maps(){
    _write_every=3000;
    _last_wrote_log=-1;
    _last_written=0;
    _ct_simplex_boundary=0;
    _ct_simplex_min=0;
    _ct_dalex=0;
    _ct_mcmc=0;
    _init_mcmc=0;
    _last_did_min=0;
    _simplex_mindex=-1;
    _calls_to_simplex_boundary=0;
    _log.set_name("carom_log");
    sprintf(_outname,"output/carom_output.sav");
    sprintf(_timingname,"output/carom_timing.sav");
    _time_started=double(time(NULL));
    _good_points.set_name("maps_good_points");
    _duds.set_name("maps_duds");
    _duds_for_min.set_name("maps_duds_for_min");
    _failed_mins.set_name("maps_failed_mins");
    _explorers.set_name("maps_explorers");
    _explorer_temp=1.0;
}

maps::~maps(){
    write_pts();

}

double maps::evaluate(array_1d<double> &pt, int *dex){
    double mu;
    _chifn.evaluate(pt,&mu,dex);

    if(mu<_chifn.target() && dex[0]>=0){
        if(_good_points.contains(dex[0])==0){
            _good_points.add(dex[0]);
        }
    }

    return mu;
}

void maps::assess_good_points(){
    assess_good_points(0,_chifn.get_pts());
}

void maps::assess_good_points(int i_min){
    assess_good_points(i_min,_chifn.get_pts());
}

void maps::assess_good_points(int i_min, int i_max){
    int i;

    for(i=0;i<_good_points.get_dim();i++){
        if(_chifn.get_fn(_good_points.get_data(i))>_chifn.target()){
            _good_points.remove(i);
            i--;
        }
    }

    for(i=i_min;i<i_max+1 && i<_chifn.get_pts();i++){
        if(_chifn.get_fn(i)<_chifn.target()){
            if(_good_points.contains(i)==0){
                _good_points.add(i);
            }
        }
    }
}

void maps::set_dof(int dd){
    _chifn.set_dof(dd);
}

void maps::set_confidence_limit(double cc){
    _chifn.set_confidence_limit(cc);
}

void maps::set_seed(int ii){
    _chifn.set_seed(ii);
}

void maps::set_min(array_1d<double> &vv){
    _chifn.set_min(vv);
}

void maps::set_max(array_1d<double> &vv){
    _chifn.set_max(vv);
}

void maps::set_characteristic_length(int dex, double vv){
    _chifn.set_characteristic_length(dex,vv);
}

void maps::set_chisquared(chisquared *xx){
    _chifn.set_chisquared(xx);
}

void maps::set_deltachi(double xx){
    _chifn.set_deltachi(xx);
}

void maps::set_target(double tt){
    _chifn.set_target(tt);
}

void maps::set_write_every(int ww){
    _write_every=ww;
}

int maps::get_called(){
    return _chifn.get_called();
}

double maps::get_chimin(){
    return _chifn.chimin();
}

void maps::set_outname(char *nn){
    int i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _outname[i]=nn[i];
    }
    _outname[i]=0;
}

void maps::set_timingname(char *nn){
    int i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _timingname[i]=nn[i];
    }
    _timingname[i]=0;
}

void maps::initialize(int npts){
    _chifn.initialize(npts);
    _cloud.build(&_chifn);
    _cloud.set_log(&_log);
    _simplex_mindex=_chifn.mindex();
    assess_good_points(0);
    _interpolator.set_kd_fn(_chifn.get_tree(), _chifn.get_fn_arr());
    _interpolator.set_ell_factor(1.0);

    array_1d<double> trial;
    int i,i_pt;
    double mu;
    while(_explorers.get_dim()<2*_chifn.get_dim()){
        for(i=0;i<_chifn.get_dim();i++){
            trial.set(i,_chifn.get_min(i)+
                        _chifn.random_double()*(_chifn.get_max(i)-_chifn.get_min(i)));
        }
        mu=evaluate(trial, &i_pt);
        if(i_pt>=0 && mu<exception_value){
            _explorers.add(i_pt);
        }
    }

    write_pts();
}


void maps::write_log(){

    FILE *output;
    char log_name[2*letters],suffix[letters];
    array_1d<int> types;
    types.set_name("carom_write_log_types");
    types.add(_log_ricochet);
    types.add(_log_mcmc);
    types.add(_log_dchi_simplex);
    types.add(_log_simplex);
    types.add(_log_compass);
    types.add(_log_swarm);
    int i,j,ii;

    for(ii=0;ii<types.get_dim();ii++){
        if(types.get_data(ii)==_log_ricochet){
            sprintf(log_name,"%s_ricochet_log.txt",_outname);
        }
        else if(types.get_data(ii)==_log_mcmc){
            sprintf(log_name,"%s_mcmc_log.txt",_outname);
        }
        else if(types.get_data(ii)==_log_dchi_simplex){
            sprintf(log_name,"%s_dchi_simplex_log.txt",_outname);
        }
        else if(types.get_data(ii)==_log_simplex){
           sprintf(log_name,"%s_simplex_log.txt",_outname);
        }
        else if(types.get_data(ii)==_log_compass){
            sprintf(log_name,"%s_compass_log.txt",_outname);
        }
        else if(types.get_data(ii)==_log_swarm){
            sprintf(log_name,"%s_swarm_log.txt",_outname);
        }
        else{
            printf("WARNING asked for unknown log type %d\n",types.get_data(ii));
            exit(1);
        }

        if(_last_wrote_log<0){
            output=fopen(log_name,"w");
        }
        else{
            output=fopen(log_name,"a");
        }

        for(i=0;i<_log.get_cols(types.get_data(ii));i++){
            for(j=0;j<_chifn.get_dim();j++){
                fprintf(output,"%e ",_chifn.get_pt(_log.get_data(types.get_data(ii),i),j));
            }
            fprintf(output,"%e %d\n",
            _chifn.get_fn(_log.get_data(types.get_data(ii),i)),
            _log.get_data(types.get_data(ii),i));
        }

        fclose(output);
    }

    _last_wrote_log=_chifn.get_pts();
    _log.reset_preserving_room();

}


void maps::write_pts(){
    FILE *output;

    int i,j;
    output=fopen(_outname,"w");
    fprintf(output,"# ");
    for(i=0;i<_chifn.get_dim();i++){
        fprintf(output,"p%d ",i);
    }
    fprintf(output,"chisq mu sig ling\n");
    for(i=0;i<_chifn.get_pts();i++){
        for(j=0;j<_chifn.get_dim();j++){
            fprintf(output,"%.18e ",_chifn.get_pt(i,j));
        }
        fprintf(output,"%.18e 0 0 0\n",_chifn.get_fn(i));
    }
    fclose(output);

    if(_last_written==0){
        output=fopen(_timingname,"w");
        fprintf(output,"#seed %d\n",_chifn.get_seed());
    }
    else{
        output=fopen(_timingname,"a");
    }

    fprintf(output,"%d %d min %.4e target %.4e -- timing -- %.4e %.4e -- %.4e %.4e -- overhead %.4e -- %d -- ",
        _chifn.get_pts(),
        _chifn.get_called(),
        _chifn.chimin(),
        _chifn.target(),
        double(time(NULL))-_time_started,
        (double(time(NULL))-_time_started)/double(_chifn.get_pts()),
        _chifn.get_time_spent(),
        _chifn.get_time_spent()/double(_chifn.get_pts()),
        (double(time(NULL))-_time_started-_chifn.get_time_spent())/double(_chifn.get_pts()),
        _calls_to_simplex_boundary);

    fprintf(output,"\n");
    fclose(output);

    write_log();
    _last_written=_chifn.get_pts();

}

void maps::simplex_min_search(){
    printf("\ndoing maps.simplex_min_search -- min %e %d\n",
    _chifn.chimin(),_chifn.get_pts());
    int pt_start=_chifn.get_pts();
    int mindex_0=_chifn.mindex();

    assess_good_points();
    int i,j;

    array_1d<int> empty;

    dchi_multimodal_simplex dchifn(&_chifn, _good_points);

    array_2d<double> seed,old_dir;
    seed.set_name("carom_simplex_min_search_seed");
    old_dir.set_name("carom_simplex_min_search_old_dir");
    array_1d<double> dir,trial,trial_dir,origin;
    dir.set_name("carom_simplex_min_search_dir");
    trial.set_name("carom_simplex_min_search_trial");
    trial_dir.set_name("carom_simplex_min_searc_trial_dir");
    double rr,mu_best,dd;

    array_1d<int> chosen_exp;

    while(seed.get_rows()<_chifn.get_dim()+1){
        i=_chifn.random_int()%_explorers.get_dim();
        if(chosen_exp.contains(i)==0){
            seed.add_row(_chifn.get_pt(_explorers.get_data(i))[0]);
            chosen_exp.add(i);
        }
    }

    printf("got seed %e\n",_chifn.chimin());

    int i_min_0=_chifn.mindex();

    simplex_minimizer ffmin;
    ffmin.set_chisquared(&dchifn);
    ffmin.set_dice(_chifn.get_dice());
    array_1d<double> min,max;
    min.set_name("carom_simplex_search_min");
    max.set_name("carom_simplex_search_min");

    for(i=0;i<_chifn.get_dim();i++){
        min.set(i,0.0);
        max.set(i,_chifn.get_characteristic_length(i));
    }

    ffmin.set_minmax(min,max);
    ffmin.use_gradient();

    array_1d<double> minpt;
    ffmin.find_minimum(seed,minpt);

    assess_good_points();

    _last_did_min=_chifn.get_pts();

    double mu_found;
    mu_found=evaluate(minpt,&i);

    if(i>=0){
        _log.add(_log_simplex, i);
    }

    if(i>=0 && _chifn.mindex()==mindex_0){
        _failed_mins.add(i);
    }

    _ct_simplex_min+=_chifn.get_pts()-pt_start;
    printf("    min %e -- %e\n",_chifn.chimin(),mu_found);

    if(_chifn.mindex()!=i_min_0){
        printf("    traveled %e\n",
        _chifn.distance(_chifn.mindex(),i_min_0));
        for(i=0;i<_chifn.get_dim();i++){
            printf("    %e %e -- %e\n",
            _chifn.get_pt(_chifn.mindex(),i),_chifn.get_pt(i_min_0,i),
            _chifn.get_characteristic_length(i));
        }
    }

    _simplex_mindex=_chifn.mindex();
}


void maps::mcmc_search(){
    printf("\ndoing maps.mcmc_search %d %e %d %d\n",
    _ct_mcmc,_chifn.chimin(),_chifn.get_called(),_chifn.get_pts());
    int pt_start=_chifn.get_pts();

    int i;
    double delta=_chifn.target()-_chifn.chimin();

    if(_init_mcmc==0){
        _mcmc.initialize(4,99,&_chifn);
        _mcmc.set_name_root(_outname);
        for(i=0;i<_chifn.get_dim();i++){
            _mcmc.set_min(i,_chifn.get_min(i));
            _mcmc.set_max(i,_chifn.get_max(i));
        }
        _mcmc.set_burnin(1000);
        _mcmc.guess_bases(_chifn.get_pt(_chifn.mindex())[0], delta, 1);
        _init_mcmc=1;
        _mcmc_basis_min=_chifn.chimin();
    }

    if(_mcmc_basis_min-_chifn.chimin()>0.1*delta){
        _mcmc.guess_bases(_chifn.get_pt(_chifn.mindex())[0], delta, 1);
        _mcmc_basis_min=_chifn.chimin();
    }

    _mcmc.sample(100*_chifn.get_dim());

    assess_good_points(pt_start);

    _ct_mcmc+=_chifn.get_pts()-pt_start;
    printf("min %e target %e\n",_chifn.chimin(),_chifn.target());
}

void maps::explore(){

    array_1d<double> min,max,norm;
    min.set_name("maps_explore_min");
    max.set_name("maps_explore_max");
    norm.set_name("maps_explore_norm");

    int ip,i,j;
    double mu;
    for(ip=0;ip<_explorers.get_dim();ip++){
        for(i=0;i<_chifn.get_dim();i++){
            mu=0.0;
            for(j=0;j<_chifn.get_dim();j++){
                mu+=_chifn.get_pt(_explorers.get_data(ip),j)*_cloud.get_basis(i,j);
            }

            if(i>=min.get_dim() || mu<min.get_data(i)){
                min.set(i,mu);
            }

            if(i>=max.get_dim() || mu>max.get_data(i)){
                max.set(i,mu);
            }
        }
    }

    for(i=0;i<_chifn.get_dim();i++){
        if(max.get_data(i)-min.get_data(i)>1.0e-20){
            norm.set(i,max.get_data(i)-min.get_data(i));
        }
        else{
            norm.set(i,1.0);
        }
    }

    array_1d<double> f_val,f_val_sorted;
    array_1d<int> f_val_dex;
    f_val.set_name("maps_explore_fval");
    f_val_sorted.set_name("maps_explore_fval_sorted");
    f_val_dex.set_name("maps_explore_f_val_dex");
    for(i=0;i<_explorers.get_dim();i++){
        f_val.add(_chifn.get_fn(_explorers.get_data(i)));
        f_val_dex.add(i);
    }
    sort_and_check(f_val,f_val_sorted,f_val_dex);

    double f_val_median=f_val_sorted.get_data(f_val_dex.get_dim()/2);

    array_1d<double> dir,trial;
    dir.set_name("maps_explore_dir");
    trial.set_name("maps_explore_trial");

    array_1d<int> accepted;
    for(i=0;i<_explorers.get_dim();i++){
        accepted.set(i,0);
    }

    double roll,ratio;
    int i_pt,accept_it;

    int n_steps=100;
    int i_step;
    double rr;
    for(i_step=0;i_step<n_steps;i_step++){
        for(ip=0;ip<_explorers.get_dim();ip++){
            for(i=0;i<_chifn.get_dim();i++){
                dir.set(i,normal_deviate(_chifn.get_dice(),0.0,1.0));
            }
            dir.normalize();
            rr=fabs(normal_deviate(_chifn.get_dice(),0.0,0.1));
            for(i=0;i<_chifn.get_dim();i++){
                trial.set(i,_chifn.get_pt(_explorers.get_data(ip),i));
            }
            for(i=0;i<_chifn.get_dim();i++){
                for(j=0;j<_chifn.get_dim();j++){
                    trial.add_val(j,rr*norm.get_data(i)*dir.get_data(i)*_cloud.get_basis(i,j));
                }
            }
            mu=evaluate(trial, &i_pt);
            accept_it=0;
            if(mu<_chifn.get_fn(_explorers.get_data(ip))){
                accept_it=1;
            }
            else{
                roll=_chifn.random_double();
                ratio=exp((_chifn.get_fn(_explorers.get_data(ip))-mu)/_explorer_temp);
                if(roll<ratio){
                    accept_it=1;
                }
            }

            if(accept_it==1){
                accepted.add_val(ip,1);
                _explorers.set(ip,i_pt);
            }
        }
    }

    int min_acc,max_acc;
    for(i=0;i<accepted.get_dim();i++){
        if(i==0 || accepted.get_data(i)<min_acc){
            min_acc=accepted.get_data(i);
        }
        if(i==0 || accepted.get_data(i)>max_acc){
            max_acc=accepted.get_data(i);
        }
    }

    if(min_acc<n_steps/3){
        _explorer_temp*=10.0;
    }
    else if(max_acc>(3*n_steps)/4){
        _explorer_temp*=0.15;
    }

    for(i=0;i<_explorers.get_dim();i++){
        ip=_explorers.get_data(i);
        if(i==0 || _chifn.get_fn(ip)<mu){
            mu=_chifn.get_fn(ip);
        }
    }

    printf("min_acc %e max_acc %e min %e -- temp %e\n",
    double(min_acc)/double(n_steps),double(max_acc)/double(n_steps),mu,_explorer_temp);

}

void maps::nested_simplex_init(){

    int n_steps=10;
    int n_internal_steps=20000;
    int n_pts=2*_chifn.get_dim();

    array_2d<double> pts;
    array_1d<int> active;
    array_2d<double> seed;
    array_1d<double> norm,smin,smax;

    int i;
    for(i=0;i<_chifn.get_dim();i++){
        norm.set(i,_chifn.get_characteristic_length(i));
        smin.set(i,0.0);
        smax.set(i,_chifn.get_characteristic_length(i));
    }

    array_1d<double> trial;
    double mu;
    while(pts.get_rows()<n_pts){
        for(i=0;i<_chifn.get_dim();i++){
            trial.set(i,_chifn.get_min(i)+_chifn.random_double()*
                       (_chifn.get_max(i)-_chifn.get_min(i)));
        }
        mu=evaluate(trial,&i);
        if(mu<exception_value && i>=0){
            pts.add_row(trial);
        }
    }


    simplex_minimizer ffmin;
    ffmin.set_chisquared(&_chifn);
    ffmin.set_minmax(smin,smax);
    ffmin.set_dice(_chifn.get_dice());
    ffmin.use_gradient();
    ffmin.set_limit(n_internal_steps);

    int i_step;
    int j;
    for(i_step=0;i_step<n_steps;i_step++){
        seed.reset_preserving_room();
        active.reset_preserving_room();
        while(active.get_dim()!=_chifn.get_dim()+1){
            i=_chifn.random_int()%pts.get_rows();
            if(active.contains(i)==0){
                active.add(i);
                seed.add_row(pts(i)[0]);
            }
        }
        ffmin.find_minimum(seed, trial);
        for(i=0;i<active.get_dim();i++){
            ffmin.get_pt(i,trial);
            for(j=0;j<_chifn.get_dim();j++){
                pts.set(active.get_data(i),j,trial.get_data(j));
            }
        }
        printf("    chimin %e\n",_chifn.chimin());
    }

    ffmin.set_limit(-1);
    seed.reset_preserving_room();
    array_1d<double> min_val,min_val_sorted;
    array_1d<int> min_val_dex;
    for(i=0;i<pts.get_rows();i++){
        min_val_dex.add(i);
        mu=evaluate(pts(i)[0],&j);
        min_val.add(mu);
    }
    sort_and_check(min_val, min_val_sorted, min_val_dex);
    for(i=0;i<_chifn.get_dim()+1;i++){
        seed.add_row(pts(min_val_dex.get_data(i))[0]);
    }
    ffmin.find_minimum(seed, trial);

}

void maps::mcmc_init(){
    int total_per=1000;
    int adjust_every=50;
    int n_groups=3;
    int n_particles=n_groups*(_chifn.get_dim()+1);

    array_1d<int> abs_min_pt;
    array_1d<int> local_min_pt;
    array_1d<int> since_min;
    array_1d<int> particles;
    array_1d<int> accepted,accepted_sorted,accepted_dex;
    array_1d<int> total_accepted;
    abs_min_pt.set_name("mcmc_init_abs_minpt");
    local_min_pt.set_name("mcmc_init_local_minpt");
    since_min.set_name("mcmc_init_since_min");
    particles.set_name("mcmc_init_particles");
    accepted.set_name("mcmc_init_accepted");
    accepted_sorted.set_name("mcmc_init_acc_sorted");
    accepted_dex.set_name("mcmc_init_acc_dex");
    total_accepted.set_name("total_accepted");

    double _temp=1.0;

    array_1d<double> trial,dir,norm;
    trial.set_name("mcmc_init_trial");
    dir.set_name("mcmc_init_dir");
    norm.set_name("mcmc_init_norm");

    double rr,re_norm;

    re_norm=1.0;

    int ip,i,j,i_step,i_found;

    for(i=0;i<_chifn.get_dim();i++){
        norm.set(i,_chifn.get_characteristic_length(i));
    }

    int adjusted=0;
    double mu,roll,ratio;
    int accept_it;
    int min_acc,max_acc,med_acc;

    array_1d<double> min_vals,min_val_sorted;
    array_1d<int> min_dexes;
    min_vals.set_name("min_vals");
    min_val_sorted.set_name("min_val_sorted");
    min_dexes.set_name("min_dexes");

    array_1d<double> dd,dd_sorted;
    array_1d<int> dd_dexes,current_particles;
    double dd_term;
    dd.set_name("mcmc_init_dd");
    dd_sorted.set_name("mcmc_init_dd_sorted");
    dd_dexes.set_name("mcmc_init_dd_dexes");
    current_particles.set_name("mcmc_init_current_particles");

    asymm_array_2d<int> trails;
    trails.set_name("mcmc_init_trails");

    for(i=0;i<n_particles;i++){
        accepted.set(i,0);
        total_accepted.set(i,0);
    }

    for(ip=0;ip<n_particles;ip++){
       i_found=-1;
        while(i_found<0){
            for(i=0;i<_chifn.get_dim();i++){
                trial.set(i,_chifn.get_min(i)+
                           _chifn.random_double()*(_chifn.get_max(i)-_chifn.get_min(i)));
            }
            mu=evaluate(trial,&i_found);

        }
        particles.set(ip,i_found);
        trails.set(ip,0,i_found);
        abs_min_pt.set(ip,i_found);
        local_min_pt.set(ip,i_found);
        since_min.set(ip,0);
    }

    double needed_temp;
    array_1d<double> needed_temp_arr,needed_temp_sorted;
    array_1d<int> needed_temp_dex;
    double old_temp;
    int has_been_adjusted;
    int step_ct=0;
    int needs_adjustment;

    array_1d<double> geo_center,local_min,local_max;
    int min_pt_connected;

    printf("starting steps with min %e\n",_chifn.chimin());
    for(i_step=0;i_step<total_per;i_step++){
        for(ip=0;ip<n_particles;ip++){

            rr=-1.0;
            while(rr<1.0e-10){
                rr=fabs(normal_deviate(_chifn.get_dice(),re_norm,0.5*re_norm));
            }

            for(i=0;i<_chifn.get_dim();i++){
                dir.set(i,normal_deviate(_chifn.get_dice(),0.0,1.0));
            }
            dir.normalize();
            for(i=0;i<_chifn.get_dim();i++){
                trial.set(i,_chifn.get_pt(particles.get_data(ip),i)+
                            rr*norm.get_data(i)*dir.get_data(i));
            }
            mu=evaluate(trial,&i_found);


            accept_it=0;

            if(ip>=abs_min_pt.get_dim() || mu<_chifn.get_fn(abs_min_pt.get_data(ip))){
                abs_min_pt.set(ip,i_found);
            }

            if(ip>=local_min_pt.get_dim() || mu<_chifn.get_fn(local_min_pt.get_data(ip))){
                local_min_pt.set(ip,i_found);
                since_min.set(ip,0);
            }
            else{
                since_min.add_val(ip,1);
            }

            if(ip>=particles.get_dim() || mu<_chifn.get_fn(particles.get_data(ip))){
                accept_it=1;
                needed_temp=1.0;
            }
            else{
                roll=_chifn.random_double();
                ratio=exp(-0.5*(mu-_chifn.get_fn(particles.get_data(ip)))/_temp);
                needed_temp=-0.5*(mu-_chifn.get_fn(particles.get_data(ip)))/log(roll);
                if(roll<ratio){
                    accept_it=1;
                }
            }

            needed_temp_dex.add(needed_temp_arr.get_dim());
            if(needed_temp<1.0){
                needed_temp_arr.add(1.0);
            }
            else{
                needed_temp_arr.add(needed_temp);
            }

            if(accept_it==1){
                particles.set(ip,i_found);
                trails.add(ip,i_found);
                accepted.add_val(ip,1);
                total_accepted.add_val(ip,1);
            }
            else{
                trails.add(ip,particles.get_data(ip));
            }

        }
        step_ct++;

        if(i_step>0 && i_step%adjust_every==0){

            accepted_sorted.reset_preserving_room();
            accepted_dex.reset_preserving_room();
            for(i=0;i<n_particles;i++){
                accepted_dex.set(i,i);
            }
            sort_and_check(accepted,accepted_sorted,accepted_dex);

            med_acc=accepted_sorted.get_data(accepted_dex.get_dim()/2);
            min_acc=accepted_sorted.get_data(0);
            max_acc=accepted_sorted.get_data(accepted_dex.get_dim()-1);

            has_been_adjusted=0;

            needs_adjustment=0;
            if(med_acc<step_ct/3){
                needs_adjustment=1;
            }
            else if(med_acc>(3*step_ct)/4){
                needs_adjustment=-1;
            }

            if(needs_adjustment!=0){
                sort_and_check(needed_temp_arr, needed_temp_sorted, needed_temp_dex);
                old_temp=_temp;
                _temp=needed_temp_sorted.get_data(needed_temp_dex.get_dim()/2);

                if(fabs(1.0-(_temp/old_temp))>0.01){
                    has_been_adjusted=1;
                }
            }
            else{
                re_norm+=0.1;
            }

            local_min.reset_preserving_room();
            local_max.reset_preserving_room();

            for(ip=0;ip<particles.get_dim();ip++){
                current_particles.set(ip,particles.get_data(ip));
                for(i=0;i<_chifn.get_dim();i++){
                    if(ip==0 || _chifn.get_pt(particles.get_data(ip),i)<local_min.get_data(i)){
                        local_min.set(i,_chifn.get_pt(particles.get_data(ip),i));
                    }
                    if(ip==0 || _chifn.get_pt(particles.get_data(ip),i)>local_max.get_data(i)){
                        local_max.set(i,_chifn.get_pt(particles.get_data(ip),i));
                    }
                }
            }

            for(ip=0;ip<particles.get_dim();ip++){
                if(since_min.get_data(ip)>=adjust_every){
                    min_pt_connected=0;
                    for(i=0;i<_chifn.get_dim();i++){
                        trial.set(i,0.5*(_chifn.get_pt(local_min_pt.get_data(ip),i)+
                                         _chifn.get_pt(_chifn.mindex(),i)));
                    }

                    mu=evaluate(trial, &i);
                    if(mu<_chifn.get_fn(abs_min_pt.get_data(ip))){
                        abs_min_pt.set(ip,i);
                    }
                    if(mu<_chifn.get_fn(local_min_pt.get_data(ip))){
                        min_pt_connected=1;
                        local_min_pt.set(ip,i);
                    }
                    if(local_min_pt.get_data(ip)==_chifn.mindex()){
                        min_pt_connected=1;
                    }

                    if(min_pt_connected==0){
                        for(i=0;i<_chifn.get_dim();i++){
                            mu=local_max.get_data(i)-local_min.get_data(i);
                            trial.set(i,local_min.get_data(i)+_chifn.random_double()*mu);
                        }
                    }
                    else{
                        for(i=0;i<_chifn.get_dim();i++){
                            geo_center.set(i,0.0);
                        }

                        for(i=0;i<particles.get_dim();i++){
                            if(i!=ip){
                                for(j=0;j<_chifn.get_dim();j++){
                                    geo_center.add_val(j,_chifn.get_pt(particles.get_data(i),j));
                                }
                            }
                        }

                        for(i=0;i<_chifn.get_dim();i++){
                            geo_center.divide_val(i,double(particles.get_dim()-1));
                        }

                        for(i=0;i<_chifn.get_dim();i++){
                            trial.set(i, 3.0*geo_center.get_data(i)
                                         -2.0*_chifn.get_pt(particles.get_data(ip),i));
                        }
                    }

                    mu=evaluate(trial,&i_found);
                    if(i_found>=0){
                        particles.set(ip,i_found);
                        trails.add(ip,i_found);
                        local_min_pt.set(ip,i_found);
                        since_min.set(ip,0);
                        if(mu<_chifn.get_fn(abs_min_pt.get_data(ip))){
                            abs_min_pt.set(ip,i_found);
                        }
                    }
                }
            }



            needed_temp_sorted.reset_preserving_room();

            //printf("    acc %d %d %d out of %d temp %e re_norm %e min %e\n",
            //min_acc,med_acc,max_acc,step_ct,_temp, re_norm, _chifn.chimin());

            if(has_been_adjusted==1){
                needed_temp_arr.reset_preserving_room();
                needed_temp_dex.reset_preserving_room();
               for(i=0;i<n_particles;i++){
                   accepted.set(i,0);
                }
                step_ct=0;
            }

            if(has_been_adjusted==1){
                adjusted++;
            }
        }
    }

    array_1d<int> connected;
    double mu_min;
    int dex_min;

    for(i=0;i<abs_min_pt.get_dim();i++){
        if(i==0 || _chifn.get_fn(abs_min_pt.get_data(i))<mu_min){
             mu_min=_chifn.get_fn(abs_min_pt.get_data(i));
             dex_min=abs_min_pt.get_data(i);
        }
    }

    double min_disconnected=2.0*exception_value;
    int n_disconnected=0;

    for(i=0;i<abs_min_pt.get_dim();i++){
        if(abs_min_pt.get_data(i)==dex_min){
            connected.set(i,1);
        }
        else{
            for(j=0;j<_chifn.get_dim();j++){
                trial.set(j,0.5*(_chifn.get_pt(dex_min,j)+
                            _chifn.get_pt(abs_min_pt.get_data(i),j)));
            }
            mu=evaluate(trial,&j);
            if(mu<_chifn.get_fn(abs_min_pt.get_data(i))){
                connected.set(i,1);
            }
            else{
                connected.set(i,0);
                n_disconnected++;
                if(_chifn.get_fn(abs_min_pt.get_data(i))<min_disconnected){
                    min_disconnected=_chifn.get_fn(abs_min_pt.get_data(i));
                }
            }
        }
    }

    for(i=0;i<n_particles;i++){
        printf("min %e %d - %e - %d\n",
        _chifn.get_fn(abs_min_pt.get_data(i)),
        total_accepted.get_data(i),
        _chifn.get_fn(local_min_pt.get_data(i)),
        connected.get_data(i));
    }
    printf("called %d -- %e\n",_chifn.get_pts(),_chifn.chimin());
    printf("min disconnected %e - %d\n",min_disconnected,n_disconnected);

    array_1d<double> smin,smax;
    for(i=0;i<_chifn.get_dim();i++){
        smin.set(i,0.0);
        smax.set(i,_chifn.get_characteristic_length(i));
    }
    simplex_minimizer ffmin;
    ffmin.set_chisquared(&_chifn);
    ffmin.set_minmax(smin,smax);
    ffmin.set_dice(_chifn.get_dice());
    ffmin.use_gradient();

    min_vals.reset();
    min_val_sorted.reset();
    min_dexes.reset();
    for(i=0;i<abs_min_pt.get_dim();i++){
        min_vals.add(_chifn.get_fn(abs_min_pt.get_data(i)));
        min_dexes.add(abs_min_pt.get_data(i));
    }
    sort_and_check(min_vals, min_val_sorted, min_dexes);
    array_2d<double> seed;
    array_1d<int> chosen;
    for(i=0;i<_chifn.get_dim()+1;i++){
        seed.add_row(_chifn.get_pt(min_dexes.get_data(i))[0]);
        chosen.add(min_dexes.get_data(i));
    }
    ffmin.find_minimum(seed,trial);

    seed.reset_preserving_room();


    min_vals.reset();
    min_val_sorted.reset();
    min_dexes.reset();
    for(i=0;i<abs_min_pt.get_dim();i++){
        if(connected.get_data(i)==0){
            min_vals.add(_chifn.get_fn(abs_min_pt.get_data(i)));
            min_dexes.add(abs_min_pt.get_data(i));
        }
    }
    sort_and_check(min_vals, min_val_sorted, min_dexes);
    for(i=0;i<min_dexes.get_dim() && seed.get_rows()!=_chifn.get_dim()+1;i++){
        seed.add_row(_chifn.get_pt(min_dexes.get_data(i))[0]);
        chosen.add(min_dexes.get_data(i));
    }

    while(seed.get_rows()!=_chifn.get_dim()+1){
        i=_chifn.random_int()%abs_min_pt.get_dim();
        if(chosen.contains(abs_min_pt.get_data(i))==0 && abs_min_pt.get_data(i)!=dex_min){
            seed.add_row(_chifn.get_pt(abs_min_pt.get_data(i))[0]);
            chosen.add(abs_min_pt.get_data(i));
        }
    }

    ffmin.find_minimum(seed,trial);

}

void maps::simplex_init(){

    int i;
    array_1d<double> smin,smax;
    for(i=0;i<_chifn.get_dim();i++){
        smin.set(i,0.0);
        smax.set(i,_chifn.get_characteristic_length(i));
    }

    simplex_minimizer ffmin;
    ffmin.set_chisquared(&_chifn);
    ffmin.set_minmax(smin,smax);

    ffmin.set_dice(_chifn.get_dice());
    ffmin.use_gradient();

    array_2d<double> seed;
    array_1d<double> min_pt;
    array_1d<int> abs_min_pt;
    array_1d<double> trial;

    int ip,j,i_found;
    double mu;
    for(ip=0;ip<2*_chifn.get_dim();ip++){
        seed.reset_preserving_room();
        for(i=0;i<_chifn.get_dim()+1;i++){
            for(j=0;j<_chifn.get_dim();j++){
                trial.set(j,_chifn.get_min(j)+
                            _chifn.random_double()*
                            (_chifn.get_max(j)-_chifn.get_min(j)));
            }
            seed.add_row(trial);
        }

        ffmin.find_minimum(seed, min_pt);
        mu=evaluate(min_pt,&i_found);
        if(i_found>=0 && abs_min_pt.contains(i_found)==0){
            abs_min_pt.add(i_found);
        }
        else{
            ip--;
        }
    }



    array_1d<double> min_vals, min_val_sorted;
    array_1d<int> min_dexes;

    for(i=0;i<abs_min_pt.get_dim();i++){
        min_vals.add(_chifn.get_fn(abs_min_pt.get_data(i)));
        min_dexes.add(abs_min_pt.get_data(i));
    }
    sort_and_check(min_vals, min_val_sorted, min_dexes);

    array_2d<double> final_seed;
    for(i=0;i<_chifn.get_dim()+1;i++){
        final_seed.add_row(_chifn.get_pt(min_dexes.get_data(i))[0]);
    }
    ffmin.find_minimum(final_seed,trial);


}



void maps::search(int limit){
    int pt_start;

    double min0=_chifn.chimin();
    printf("before init min %e\n",_chifn.chimin());
    mcmc_init();
    printf("min now %e -> %e\n",min0,_chifn.chimin());
    printf("called %d\n",_chifn.get_pts());
    exit(1);

    while(_chifn.get_pts()<limit){

        pt_start=_chifn.get_pts();
        _cloud.search();
        _ct_dalex+=_chifn.get_pts()-pt_start;

        if(_chifn.get_pts()-_last_written>_write_every){
            write_pts();
        }
    }

    int i;
    printf("minpt -- %e\n",_chifn.get_fn(_chifn.mindex()));
    for(i=0;i<_chifn.get_dim();i++){
        printf("    %.3e\n",_chifn.get_pt(_chifn.mindex(),i));
    }
    printf("\n\n");
    write_pts();
}


int maps::bisection(int ilow, int ihigh, double target, double tol){
    return bisection(_chifn.get_pt(ilow)[0], _chifn.get_pt(ihigh)[0], target, tol);
}


int maps::bisection(int iOrigin, array_1d<double> &dir, double target, double tol){

    array_1d<double> highball;
    highball.set_name("maps_bisection_dir_highball");
    double norm=1.0;
    double fhigh=-2.0*exception_value;
    int i;
    for(i=0;i<_chifn.get_dim();i++){
        highball.set(i,_chifn.get_pt(iOrigin,i));
    }

    while(fhigh<target+tol){
        for(i=0;i<_chifn.get_dim();i++){
            highball.add_val(i,norm*dir.get_data(i));
        }
        fhigh=evaluate(highball,&i);
        norm*=2.0;
    }

    return bisection(_chifn.get_pt(iOrigin)[0], highball, target, tol);
}

int maps::bisection(array_1d<double> &low_in, array_1d<double> &high_in,
                    double target, double tol){


    int i;
    array_1d<double> trial,lowball,highball;
    trial.set_name("maps_bisection_trial");
    lowball.set_name("maps_bisection_trial");
    highball.set_name("maps_bisection_trial");

    for(i=0;i<_chifn.get_dim();i++){
        lowball.set(i,low_in.get_data(i));
        highball.set(i,high_in.get_data(i));
    }

    int iFound,iTrial;
    double mu,dmu,muFound,dmuFound;

    mu=evaluate(highball,&iFound);
    muFound=evaluate(lowball, &iFound);

    dmuFound=fabs(muFound-target);

    int ct;

    for(ct=0;(ct<5 || dmuFound>tol) && ct<40;ct++){
        for(i=0;i<_chifn.get_dim();i++){
            trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
        }
        mu=evaluate(trial,&iTrial);
        if(mu<target){
            for(i=0;i<_chifn.get_dim();i++){
                lowball.set(i,trial.get_data(i));
            }
        }
        else{
            for(i=0;i<_chifn.get_dim();i++){
                highball.set(i,trial.get_data(i));
            }
        }

        dmu=fabs(mu-target);
        if(dmu<dmuFound){
            iFound=iTrial;
            muFound=mu;
            dmuFound=dmu;
        }
    }

    return iFound;

}
