#include "mcmc.h"

void mcmc_sampler::write_chains(){
    int i,j;
    char out_name[letters];
    FILE *output;
    if(_has_written==0){
        for(i=0;i<_n_chains;i++){
            sprintf(out_name,"%s_%d.txt",_name_root,i);
            output=fopen(out_name,"w");
            fprintf(output,"# seed %d\n",_seed);
            fprintf(output,"# degen chisq ");
            for(j=0;j<_dim;j++){
                fprintf(output,"p%d ",j);
            }
            fprintf(output,"\n");
            fclose(output);
        }
    }

    _has_written=1;

    array_2d<double> local_points;
    local_points.set_name("mcmc_write_chains_local_points");
    array_1d<double> local_chisq;
    local_chisq.set_name("mcmc_write_chains_local_chisq");
    array_1d<int> local_degen;
    local_degen.set_name("mcmc_write_chains_local_degen");
    int dex;
    int i_chain;
    int k;

    array_2d<double> start_pts;
    start_pts.set_name("mcmc_write_chains_start_pts");
    array_1d<double> start_chisq;
    start_chisq.set_name("mcmc_write_chains_start_chisq");
    array_1d<int> start_degen;
    start_degen.set_name("mcmc_write_chains_start_degen");

    for(i_chain=0;i_chain<_n_chains;i_chain++){
        local_degen.reset_preserving_room();
        local_points.reset_preserving_room();
        local_chisq.reset_preserving_room();

        for(i=0;i<_point_dexes.get_cols(i_chain);i++){
            dex=_point_dexes.get_data(i_chain,i);
            local_points.add_row(_points(dex));
            local_chisq.add(_chisq_arr.get_data(dex));
            j=1;
            for(k=i+1;
                k<_point_dexes.get_cols(i_chain) &&
                _point_dexes.get_data(i_chain,k)==dex;
                k++){

                j++;

            }
            local_degen.add(j);
            i+=j-1;
        }

        sprintf(out_name,"%s_%d.txt",_name_root,i_chain);
        output=fopen(out_name,"a");

        for(i=0;i<local_points.get_rows()-1;i++){
            fprintf(output,"%d %e",local_degen.get_data(i),local_chisq.get_data(i));
            for(j=0;j<_dim;j++){
                fprintf(output," %e",local_points.get_data(i,j));
            }
            fprintf(output,"\n");
        }
        fclose(output);
        start_pts.add_row(local_points(i));
        start_chisq.add(local_chisq.get_data(i));
        start_degen.add(local_degen.get_data(i));
    }

    _points.reset_preserving_room();
    _point_dexes.reset_preserving_room();
    _chisq_arr.reset_preserving_room();

    for(i=0;i<_n_chains;i++){
        _points.add_row(start_pts(i));
        _chisq_arr.add(start_chisq.get_data(i));
        for(j=0;j<start_degen.get_data(i);j++){
            _point_dexes.set(i,j,_points.get_rows()-1);
        }
    }

}

void mcmc_sampler::sample(int steps_per_chain){
    if(_chisq_fn==NULL){
        printf("WARNING; _chisq_fn is NULL in sample\n");
        exit(1);
    }

    if(_dice==NULL){
        printf("WARNING; _dice is NULL in sample\n");
        exit(1);
    }

    if(_points.get_rows()==0){
        printf("WARNING; sampling, but haven't set start points\n");
        exit(1);
    }

    double chi_old,chi_new;
    int i_step;
    int i_chain;
    int i_basis;
    double rr;
    array_1d<double> new_pt;
    new_pt.set_name("mcmc_sample_new_pt");
    int old_dex;
    int total_ct=0;
    int i,j,accept_it;
    array_1d<double> dir;
    dir.set_name("mcmc_sample_dir");
    for(i_step=0;i_step<steps_per_chain;i_step++){
        for(i_chain=0;i_chain<_n_chains;i_chain++){
            total_ct++;
            accept_it=0;
            old_dex=_point_dexes.get_data(i_chain,_point_dexes.get_cols(i_chain)-1);
            chi_old=_chisq_arr.get_data(old_dex);
            i_basis=_dice->int32()%_dim;
            rr=fabs(normal_deviate(_dice,1.0,0.25));
            for(i=0;i<_dim;i++){
                dir.set(i,normal_deviate(_dice,0.0,1.0));
            }
            dir.normalize();
            for(i=0;i<_dim;i++){
                new_pt.set(i,_points.get_data(old_dex,i));
            }
            for(i=0;i<_dim;i++){
                 for(j=0;j<_dim;j++){
                     new_pt.add_val(j,rr*dir.get_data(i)*_bases.get_data(i,j));
                 }
            }
            chi_new=_chisq_fn[0](new_pt);
            if(chi_new<chi_old){
                accept_it=1;
            }
            else{
                rr=_dice->doub();
                if(rr<exp(-0.5*(chi_new-chi_old))){
                    accept_it=1;
                }
            }

            if(accept_it==1){
                _points.add_row(new_pt);
                _chisq_arr.add(chi_new);
                _point_dexes.add(i_chain,_points.get_rows()-1);
            }
            else{
                _point_dexes.add(i_chain,old_dex);
            }
        }

        if(total_ct>=_write_every){
            write_chains();
        }
    }
}
