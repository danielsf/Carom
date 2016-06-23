#include "maps.h"

maps::maps(){
    _write_every=3000;
    _last_wrote_log=-1;
    _last_written=0;
    _ct_simplex_boundary=0;
    _ct_simplex_min=0;
    _ct_dalex=0;
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

void maps::mcmc_init(){
    maps_initializer initializer;
    initializer.set_chifn(&_chifn);
    initializer.search();
}

void maps::search(int limit){
    int pt_start;

    double min0=_chifn.chimin();
    printf("before init min %e\n",_chifn.chimin());
    mcmc_init();
    printf("min now %e -> %e\n",min0,_chifn.chimin());
    printf("called %d\n",_chifn.get_pts());

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
