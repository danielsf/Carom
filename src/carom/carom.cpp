#include "carom.h"

carom::carom(){
    _write_every=3000;
    _last_wrote_log=-1;
    _last_written=0;
    _ct_simplex=0;
    _ct_node=0;
    _unique_nodes=0;
    _calls_to_simplex=0;
    _log.set_name("carom_log");
    sprintf(_outname,"output/carom_output.sav");
    sprintf(_timingname,"output/carom_timing.sav");
    _time_started=double(time(NULL));
}

carom::~carom(){
    write_pts();

}

void carom::set_dof(int dd){
    _chifn.set_dof(dd);
}

void carom::set_confidence_limit(double cc){
    _chifn.set_confidence_limit(cc);
}

void carom::set_seed(int ii){
    _chifn.set_seed(ii);
}

void carom::set_min(array_1d<double> &vv){
    _chifn.set_min(vv);
}

void carom::set_max(array_1d<double> &vv){
    _chifn.set_max(vv);
}

void carom::set_characteristic_length(int dex, double vv){
    _chifn.set_characteristic_length(dex,vv);
}

void carom::set_chisquared(chisquared *xx){
    _chifn.set_chisquared(xx);
}

void carom::set_deltachi(double xx){
    _chifn.set_deltachi(xx);
}

void carom::set_target(double tt){
    _chifn.set_target(tt);
}

void carom::set_write_every(int ww){
    _write_every=ww;
}

int carom::get_called(){
    return _chifn.get_called();
}

int carom::active_nodes(){
    int ix,i=0;
    for(ix=0;ix<_nodes.get_dim();ix++){
        if(_nodes(ix)->get_activity()==1)i++;
    }

    return i;
}

void carom::set_outname(char *nn){
    int i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _outname[i]=nn[i];
    }
    _outname[i]=0;
}

void carom::set_timingname(char *nn){
    int i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _timingname[i]=nn[i];
    }
    _timingname[i]=0;
}

void carom::initialize(int npts){
    _chifn.initialize(npts);
    write_pts();
}


void carom::write_log(){

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


void carom::write_pts(){
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
        fprintf(output,"seed %d\n",_chifn.get_seed());
    }
    else{
        output=fopen(_timingname,"a");
    }

    fprintf(output,"%d min %.4e -- timing -- %.4e %.4e -- %.4e %.4e -- overhead %.4e -- %d %d %d -- ",
        _chifn.get_pts(),
        _chifn.chimin(),
        double(time(NULL))-_time_started,
        (double(time(NULL))-_time_started)/double(_chifn.get_pts()),
        _chifn.get_time_spent(),
        _chifn.get_time_spent()/double(_chifn.get_pts()),
        (double(time(NULL))-_time_started-_chifn.get_time_spent())/double(_chifn.get_pts()),
        _calls_to_simplex,_nodes.get_dim(),_unique_nodes);
    for(i=0;i<_nodes.get_dim();i++){
        fprintf(output,"%.4e %.4e %d %d -- %.4e %.4e %.4e %.4e -- convergence %d swarm expand %d",
        _nodes(i)->projected_volume(),
        _nodes(i)->volume(),
        _nodes(i)->get_n_particles(),
        _chifn.get_pts(),
        _nodes(i)->get_ricochet_growth(),
        _nodes(i)->get_mcmc_growth(),
        _nodes(i)->get_swarm_growth(),
        _nodes(i)->get_simplex_growth(),
        _nodes(i)->get_convergence_ct(),
        _nodes(i)->get_swarm_expand());

        fprintf(output," trimmed %d",_nodes(i)->get_total_trimmed());

        fprintf(output,"; ");
    }

    fprintf(output,"\n");
    fclose(output);

    printf("\nNODE CENTERS\n");
    for(i=0;i<_nodes.get_dim();i++){
        printf("    %e -- %d %e %e\n",
        _chifn.get_fn(_nodes(i)->get_center()),_nodes(i)->get_center(),
        _chifn.get_pt(_nodes(i)->get_center(),0),
        _chifn.get_pt(_nodes(i)->get_center(),1));

        if(_nodes(i)->get_total_bisections()>0){
            printf("    bisections %d %d %d\n",
            _nodes(i)->get_total_bisections(),
            _nodes(i)->get_bisection_calls(),
            _nodes(i)->get_bisection_calls()/_nodes(i)->get_total_bisections());
        }

        if(_nodes(i)->get_total_ricochets()>0){
            printf("    ricochets %d %d %d\n",
            _nodes(i)->get_total_ricochets(),
            _nodes(i)->get_ricochet_calls(),
            _nodes(i)->get_ricochet_calls()/_nodes(i)->get_total_ricochets());

            printf("    ricochet bisections %d %d %d\n",
            _nodes(i)->get_ricochet_bisections(),
            _nodes(i)->get_ricochet_bisection_calls(),
            _nodes(i)->get_ricochet_bisection_calls()/_nodes(i)->get_ricochet_bisections());

            printf("    gradient calls %d\n",
            _nodes(i)->get_gradient_calls());

            printf("    highball calls %d\n",
            _nodes(i)->get_highball_calls());
        }
        printf("    transform\n");
        for(j=0;j<_chifn.get_dim();j++){
            printf("    %e\n",_nodes(i)->get_transform(j));
        }
        printf("\n");
    }

    write_log();
    _last_written=_chifn.get_called();

}

void carom::simplex_search(){
    _calls_to_simplex++;
    int ibefore=_chifn.get_called();

    array_1d<int> local_associates;
    array_1d<double> norm;
    local_associates.set_name("carom_simplex_local_associates");
    norm.set_name("carom_simplex_norm");

    int i_node,i_pt;
    int i,j;
    double xmin,xmax,xx;

    if(_nodes.get_dim()>0){
        for(i_node=0;i_node<_nodes.get_dim();i_node++){
            i_pt=_nodes(i_node)->get_center();
            if(local_associates.contains(i_pt)==0){
                local_associates.add(i_pt);
            }

            for(i=0;i<_nodes(i_node)->get_n_boundary();i++){
                i_pt=_nodes(i_node)->get_boundary(i);
                if(local_associates.contains(i_pt)==0){
                    local_associates.add(i_pt);
                }
            }

            for(i=0;i<_nodes(i_node)->get_n_associates();i+=3){
                i_pt=_nodes(i_node)->get_associate(i);
                if(local_associates.contains(i_pt)==0){
                    local_associates.add(i_pt);
                }
            }
        }

        if(_nodes.get_dim()>1){
            for(i_node=0;i_node<_nodes.get_dim();i_node++){
                for(i=0;i<_chifn.get_dim();i++){
                    i_pt=_nodes(i_node)->get_center();
                    xmin=_chifn.get_pt(i_pt,i);
                    xmax=_chifn.get_pt(i_pt,i);
                    for(j=0;j<_nodes(i_node)->get_n_boundary();j++){
                        i_pt=_nodes(i_node)->get_boundary(j);
                        xx=_chifn.get_pt(i_pt,i);

                        if(xx<xmin){
                            xmin=xx;
                        }

                        if(xx>xmax){
                            xmax=xx;
                        }
                    }

                    if(i>=norm.get_dim() || xmax-xmin>norm.get_data(i)){
                        norm.set(i,xmax-xmin);
                    }
                }
            }
        }
    }

    dchi_multimodal_simplex dchifn(&_chifn,local_associates);

    if(norm.get_dim()>0){
        dchifn.set_norm(norm);
    }

    simplex_minimizer ffmin;
    ffmin.set_chisquared(&dchifn);
    ffmin.set_dice(_chifn.get_dice());
    array_1d<double> min,max;
    min.set_name("carom_simplex_search_min");
    max.set_name("carom_simplex_search_min");
    _chifn.get_min(min);
    _chifn.get_max(max);
    ffmin.set_minmax(min,max);
    ffmin.use_gradient();

    array_1d<double> minpt;
    minpt.set_name("carom_simplex_search_minpt");

    array_2d<double> seed;
    seed.set_name("carom_simplex_search_seed");

    seed.set_cols(_chifn.get_dim());
    int iFound;
    array_1d<double> trial;
    array_1d<int> seed_dex;
    double ftrial;
    trial.set_name("carom_simplex_search_trial");
    seed_dex.set_name("carom_simplex_search_seed_dex");

    while(seed.get_rows()<_chifn.get_dim()+1){
        for(i=0;i<_chifn.get_dim();i++){
            trial.set(i,min.get_data(i)+_chifn.random_double()*(max.get_data(i)-min.get_data(i)));
        }
        _chifn.evaluate(trial,&ftrial,&iFound);

        if(ftrial<exception_value){
            j=1;
            for(i=0;i<seed_dex.get_dim() && j==1;i++){
                if(iFound==seed_dex.get_data(i))j=0;
            }

            if(j==1){
                seed_dex.add(iFound);
                seed.add_row(trial);
            }
        }
    }

    gp_cost cost_fn;
    chisq_wrapper cost_chi;

    if(_nodes.get_dim()>0){
        cost_chi.copy(_chifn);
        cost_fn.set_chifn(&cost_chi);
        ffmin.set_cost(&cost_fn);
    }

    ffmin.find_minimum(seed,minpt);

    array_1d<int> neigh;
    array_1d<double> dd;
    neigh.set_name("carom_simplex_neigh");
    dd.set_name("carom_simplex_dd");
    _chifn.nn_srch(minpt,1,neigh,dd);

    _ct_simplex+=_chifn.get_called()-ibefore;

    _log.add(_log_simplex,neigh.get_data(0));

    assess_node(neigh.get_data(0));

    printf("done simplex searching; called cost %d _nodes %d\n",cost_fn.get_called(),_nodes.get_dim());

    write_pts();
}

void carom::search(int limit){
    int before,i;
    int active_nodes=0,goon=1,dosimplex;

    while(goon==1){
        _nodes.cull();

        if(_calls_to_simplex>_unique_nodes+2 &&
           _nodes.get_dim()>0){
            dosimplex=0;
        }
        else{
            if(active_nodes==0 || _ct_node>_ct_simplex){
                dosimplex=1;
            }
            else{
                dosimplex=0;
            }
        }

        if(dosimplex==1){
            simplex_search();
        }
        else{
            for(i=0;i<_nodes.get_dim();i++){
                if(_nodes(i)->get_activity()==1){
                    _nodes(i)->search();
                }
            }
        }

        if(_chifn.get_called()-_last_written>_write_every){
            write_pts();
        }

        active_nodes=0;
        _ct_node=0;
        for(i=0;i<_nodes.get_dim();i++){
            if(_nodes(i)->get_activity()==1){
                active_nodes++;
            }
            _ct_node+=_nodes(i)->get_ct_ricochet();
        }

        if(active_nodes==0 &&
           _calls_to_simplex>_unique_nodes+2 &&
           _nodes.get_dim()>0){
            goon=0;
            write_pts();
        }

        if(limit>0 && _chifn.get_pts()>limit){
            goon=0;
            write_pts();
        }
    }
}

void carom::assess_node(int dex){
    if(dex<0 || dex>_chifn.get_pts()){
        printf("WARNING asking to assess node %d but only have %d pts\n",
        dex,_chifn.get_pts());

        exit(1);
    }

    if(_chifn.get_fn(dex)>_chifn.target()){
        return;
    }

    if(_nodes.get_dim()==0 && _chifn.get_fn(dex)<_chifn.target()){
        _nodes.add(dex,&_chifn,&_log);
        _unique_nodes++;
        return;
    }

    int keep_it,i,ix,iFound,isAssociate,is_unique;
    double ftrial;
    array_1d<double> trial;
    trial.set_name("carom_assess_node_trial");

    keep_it=1;
    for(ix=0;ix<_nodes.get_dim() && keep_it==1;ix++){

        isAssociate=_nodes(ix)->is_this_an_associate_gross(dex);

        if(isAssociate==1){
            keep_it=0;
        }

        if(keep_it==0){
            if(_chifn.get_fn(dex)<_chifn.get_fn(_nodes(ix)->get_center())){
                _nodes(ix)->set_center(dex);
            }
        }
    }

    if(keep_it==1){
        is_unique=1;
        for(ix=0;ix<_nodes.get_dim() && is_unique==1;ix++){
            if(_nodes(ix)->is_this_an_associate(dex)==1){
                is_unique=0;
            }
        }
        _nodes.add(dex,&_chifn,&_log);
        if(is_unique==1){
            _unique_nodes++;
        }
    }

}


///////////////////////gp_cost

gp_cost::gp_cost(){
    _called=0;
    _ell=-1.0;
    _covarin.set_name("gp_cost_covarin");
    _covar.set_name("gp_cost_covar");
    _neigh.set_name("gp_cost_neigh");
    _neigh_buff.set_name("gp_cost_neigh");
    _qq.set_name("gp_cost_qq");
    _dd.set_name("gp_cost_dd");
    _chifn=NULL;
}

gp_cost::~gp_cost(){}

void gp_cost::set_chifn(chisq_wrapper *cc){
    _chifn=cc;
}

void gp_cost::is_it_safe(char *word){
    if(_chifn==NULL){
        printf("WARNING in gp_cost::%s\n",word);
        printf("_chifn is null\n");
        exit(1);
    }
}

int gp_cost::get_called(){
    return _called;
}

double gp_cost::operator()(array_1d<double> &pt){
    is_it_safe("operator()");
    _called++;

    int npts,dosrch=0;

    npts=_chifn->get_dim();
    _chifn->nn_srch(pt,npts,_neigh_buff,_dd);

    if(_ell<=0.0){
        dosrch=1;
    }

    if(_neigh.get_dim()!=_neigh_buff.get_dim()){
        dosrch=1;
    }

    if(_covarin.get_rows()!=_neigh_buff.get_dim() || _covarin.get_cols()!=_neigh_buff.get_dim()){
        dosrch=1;
    }

    if(_qq.get_dim()!=_neigh_buff.get_dim()){
        dosrch=1;
    }

    int i,j;
    double nugget=1.0e-4;
    if(dosrch==0){
       for(i=0;i<_neigh_buff.get_dim() && dosrch==0;i++){
           dosrch=1;
           for(j=0;j<_neigh.get_dim();j++){
               if(_neigh.get_data(j)==_neigh_buff.get_data(i)){
                   dosrch=0;
               }
           }
       }
    }

    array_1d<double> dd,dd_sorted;
    array_1d<int> dexes;
    int ct;

    dd.set_name("gp_cost_operator_dd");
    dd_sorted.set_name("gp_cost_operator_dd_sorted");
    dexes.set_name("gp_cost_operator_dexes");

    if(dosrch==1){
        ct=0;
        for(i=0;i<_neigh_buff.get_dim();i++){
            for(j=i+1;j<_neigh_buff.get_dim();j++){
                dd.set(ct,_chifn->distance(_neigh_buff.get_data(i),_neigh_buff.get_data(j)));
                dexes.set(ct,ct);
                ct++;
            }
        }

        sort_and_check(dd,dd_sorted,dexes);
        _ell=dd_sorted.get_data(ct/2);

        _fbar=0.0;
        for(i=0;i<_neigh_buff.get_dim();i++){
            _neigh.set(i,_neigh_buff.get_data(i));
            _fbar+=_chifn->get_fn(_neigh_buff.get_data(i));
        }
        _fbar=_fbar/double(_neigh_buff.get_dim());

        _covar.set_cols(_neigh_buff.get_dim());

        for(i=0;i<_neigh_buff.get_dim();i++){
            for(j=i;j<_neigh_buff.get_dim();j++){
                _covar.set(i,j,exp(-0.5*power(_chifn->distance(_neigh_buff.get_data(i),_neigh_buff.get_data(j))/_ell,2)));
                if(i==j){
                    _covar.add_val(i,j,nugget);
                }
                else{
                    _covar.set(j,i,_covar.get_data(i,j));
                }
            }

        }

        invert_lapack(_covar,_covarin,0);
    }

    for(i=0;i<_neigh.get_dim();i++){
        _qq.set(i,exp(-0.5*power(_chifn->distance(pt,_neigh.get_data(i))/_ell,2)));
    }

    double mu=-1.0*_fbar;
    for(i=0;i<_neigh.get_dim();i++){
        for(j=0;j<_neigh.get_dim();j++){
            mu-=_qq.get_data(i)*_covarin.get_data(i,j)*(_chifn->get_fn(_neigh.get_data(i))-_fbar);
        }
    }

    return mu;
}
