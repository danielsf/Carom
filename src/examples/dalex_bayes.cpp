#include "containers.h"
#include "goto_tools.h"
#include "kd.h"

int main(int iargc, char *argv[]){

    int i,j,dim;
    char in_name[letters];
    char out_name[letters];
    in_name[0]=0;
    out_name[0]=0;
    dim=-1;
    for(i=1;i<iargc;i++){
        if(argv[i][0]=='-'){
            switch(argv[i][1]){
                case 'i':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        in_name[j]=argv[i][j];
                    }
                    in_name[j]=0;
                    break;
                case 'o':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        out_name[j]=argv[i][j];
                    }
                    out_name[j]=0;
                    break;
                case 'd':
                    i++;
                    dim=atoi(argv[i]);
                    break;
            }
        }
    }

    if(dim<0){
        printf("need to specify dim\n");
        exit(1);
    }
    if(in_name[0]==0){
        printf("need to specify in_name\n");
        exit(1);
    }
    if(out_name[0]==0){
        printf("need to specify out_name\n");
        exit(1);
    }

    int n_cols=0;
    char word[letters];
    word[0]=0;
    FILE *in_file;
    in_file = fopen(in_name, "r");
    while(compare_char("log", word)==0){
        fscanf(in_file,"%s",word);
        if(compare_char("#",word)==0){
            n_cols++;
        }
    }

    array_2d<double> dalex_pts;
    array_1d<double> dalex_chisq;
    array_1d<double> pt;
    double xx;
    dalex_pts.set_name("dalex_pts");
    dalex_chisq.set_name("dalex_chisq");
    pt.set_name("pt");

    printf("n_cols %d\n",n_cols);
    while(fscanf(in_file,"%le",&xx)>0){
        pt.set(0,xx);
        for(i=1;i<dim;i++){
            fscanf(in_file,"%le",&xx);
            pt.set(i,xx);
        }
        dalex_pts.add_row(pt);
        fscanf(in_file,"%le",&xx);
        dalex_chisq.add(xx);
        for(i=dim+1;i<n_cols;i++){
            fscanf(in_file,"%le",&xx);
        }
    }

    fclose(in_file);

    printf("dalex_pts %d %d; dalex_chisq %d\n",
    dalex_pts.get_rows(),dalex_pts.get_cols(),dalex_chisq.get_dim());

    if(dalex_pts.get_rows()!=dalex_chisq.get_dim()){
        printf("somehow got %d pts but %d chisq\n",
        dalex_pts.get_rows(),dalex_chisq.get_dim());

        exit(1);
    }

    double chisq_min=exception_value;
    for(i=0;i<dalex_chisq.get_dim();i++){
        if(dalex_chisq.get_data(i)<chisq_min){
            chisq_min=dalex_chisq.get_data(i);
        }
    }

    array_1d<double> xmin,xmax;
    xmin.set_name("xmin");
    xmax.set_name("xmax");
    for(i=0;i<dalex_pts.get_rows();i++){
        if(dalex_chisq.get_data(i)<chisq_min+21.03){
            for(j=0;j<dim;j++){
                if(j>=xmin.get_dim() || dalex_pts.get_data(i,j)<xmin.get_data(j)){
                    xmin.set(j,dalex_pts.get_data(i,j));
                }
                if(j>=xmax.get_dim() || dalex_pts.get_data(i,j)>xmax.get_data(j)){
                    xmax.set(j,dalex_pts.get_data(i,j));
                }
            }
        }
    }

    kd_tree dalex_tree(dalex_pts, xmin, xmax);

    array_1d<int> neigh;
    array_1d<double> dist;
    neigh.set_name("neigh");
    dist.set_name("dist");

    array_1d<double> delta;
    array_1d<double> local_min,local_max;
    array_1d<int> local_min_set,local_max_set;
    array_2d<double> box_min,box_max;
    local_min.set_name("local_min");
    local_max.set_name("local_max");
    local_min_set.set_name("local_min_set");
    local_max_set.set_name("local_max_set");
    box_min.set_name("box_min");
    box_max.set_name("box_max");
    delta.set_name("delta");

    int n_neigh_start=3*dim+1;
    int n_neigh;
    int n_neigh_validate=20*dim;
    int mins_set;
    int maxes_set;
    int k;
    int dim_dex;
    int neigh_dex;
    int is_valid;
    int n_pass;
    int max_n_neigh=-1;
    int is_inside;
    double vol,min_vol=2.0*exception_value,mean_vol=0.0,var_vol=0.0;

    array_1d<double> sorted_delta;
    array_1d<int> sorted_delta_dex;
    sorted_delta.set_name("sorted_delta");
    sorted_delta_dex.set_name("sorted_delta_dex");

    double t_start = double(time(NULL));

    for(i=0;i<dalex_pts.get_rows();i++){
        for(j=0;j<dim;j++){
            local_min_set.set(j,0);
            local_max_set.set(j,0);
            local_min.set(j,dalex_pts.get_data(i,j));
            local_max.set(j,dalex_pts.get_data(i,j));
        }
        mins_set=0;
        maxes_set=0;
        n_neigh=n_neigh_validate;
        is_valid=0;
        while(is_valid==0){
            n_pass=0;
            dalex_tree.nn_srch(i,n_neigh,neigh,dist);
            if(neigh.get_data(0)!=i){
                printf("neighbor search did not find self\n");
                exit(1);
            }
            delta.reset_preserving_room();
            while(is_valid==0 && n_pass<2){
                for(j=1;j<n_neigh && mins_set<dim && maxes_set<dim; j++){
                    neigh_dex = neigh.get_data(j);
                    sorted_delta.reset_preserving_room();
                    sorted_delta_dex.reset_preserving_room();
                    for(k=0;k<dim;k++){
                        delta.set(k,fabs(dalex_pts.get_data(i,k)-dalex_pts.get_data(neigh_dex,k))/(xmax.get_data(k)-xmin.get_data(k)));
                        sorted_delta_dex.set(k,k);
                    }
                    sort(delta,sorted_delta,sorted_delta_dex);

                    for(k=dim-1;k>=0;k--){
                        dim_dex=sorted_delta_dex.get_data(k);
                        xx = dalex_pts.get_data(neigh_dex,dim_dex);
                        if(xx>dalex_pts.get_data(i,dim_dex)
                           && local_max_set.get_data(dim_dex)==0){

                            local_max.set(dim_dex, 0.5*(xx+dalex_pts.get_data(i,dim_dex)));
                            maxes_set++;
                            local_max_set.set(dim_dex,1);
                            if(n_pass==0){
                                break;
                            }
                        }
                        else if(xx<dalex_pts.get_data(i,dim_dex)
                                && local_min_set.get_data(dim_dex)==0){

                            local_min.set(dim_dex, 0.5*(xx+dalex_pts.get_data(i,dim_dex)));
                            mins_set++;
                            local_min_set.set(dim_dex,1);
                            if(n_pass==0){
                                break;
                            }
                        }
                    }
                }
                is_valid=1;
                for(j=0;j<dim;j++){
                    if(local_min_set.get_data(j)==0 && local_max_set.get_data(j)==0){
                        is_valid=0;
                        n_pass++;
                        if(n_pass>=2){
                            n_neigh+=2*dim;
                        }
                        break;
                    }
                }
            }
        }

        if(n_neigh>n_neigh_validate){
            n_neigh=n_neigh+10*dim;
            dalex_tree.nn_srch(i,n_neigh,neigh,dist);
        }
        is_valid=0;
        delta.reset_preserving_room();
        sorted_delta.reset_preserving_room();
        sorted_delta_dex.reset_preserving_room();
        while(is_valid==0){
            is_valid=1;
            for(j=1;j<n_neigh;j++){
                is_inside=1;
                for(k=0;k<dim;k++){
                    delta.set(k,2.0*exception_value);
                    sorted_delta_dex.set(k,k);
                    if(!(dalex_pts.get_data(neigh.get_data(j),k)>local_min.get_data(k) &&
                         dalex_pts.get_data(neigh.get_data(j),k)<local_max.get_data(k))){
                        is_inside=0;
                    }
                    else{
                        if(dalex_pts.get_data(neigh.get_data(j),k)<dalex_pts.get_data(i,k)){
                            delta.set(k,dalex_pts.get_data(neigh.get_data(j),k)-local_min.get_data(k));
                        }
                        else{
                            delta.set(k,local_max.get_data(k)-dalex_pts.get_data(neigh.get_data(j),k));
                        }
                    }
                }

                if(is_inside==1){
                    is_valid=0;
                    sort(delta,sorted_delta,sorted_delta_dex);
                    dim_dex = sorted_delta_dex.get_data(0);
                    xx = dalex_pts.get_data(neigh.get_data(j),dim_dex);
                    if(xx<dalex_pts.get_data(i,dim_dex)){
                        local_min.set(dim_dex,0.5*(xx+dalex_pts.get_data(i,dim_dex)));
                    }
                    else{
                        local_max.set(dim_dex,0.5*(xx+dalex_pts.get_data(i,dim_dex)));
                    }
                }

            }
        }

        if(dalex_chisq.get_data(i)<chisq_min+200.0){
            for(j=0;j<dim;j++){
                if(local_min_set.get_data(j)==0 && local_max_set.get_data(j)==0){
                    printf("%d %e somehow, volume will be zero\n",i,dalex_chisq.get_data(i));
                    printf("%e\n",dalex_pts.get_data(i,j));
                    for(k=0;k<n_neigh;k++){
                         printf("    %e %e\n",
                         dalex_pts.get_data(neigh.get_data(k),j)-dalex_pts.get_data(i,j),
                         dalex_chisq.get_data(i)-dalex_chisq.get_data(neigh.get_data(k)));
                    }
                    exit(1);
                }
                if(local_min.get_data(j)>dalex_pts.get_data(i,j) ||
                   local_max.get_data(j)<dalex_pts.get_data(i,j)){

                   printf("%d somehow %e < %e < %e, which is not true\n",
                   i, local_min.get_data(j),dalex_pts.get_data(i,j),
                   local_max.get_data(j));
               }
            }

            for(j=1;j<n_neigh;j++){
                neigh_dex=neigh.get_data(j);
                is_inside=1;
                for(k=0;k<dim;k++){
                    xx=dalex_pts.get_data(neigh_dex,k);
                    if(!(xx>local_min.get_data(k) && xx<local_max.get_data(k))){
                        is_inside=0;
                        break;
                    }
                }

                if(is_inside==1){
                    printf("a neighbor is inside the box %d\n",j);
                    for(k=0;k<dim;k++){
                        printf("%e ?< %e ?< %e; %e\n",
                        local_min.get_data(k),
                        dalex_pts.get_data(neigh_dex,k),
                        local_max.get_data(k),
                        dalex_pts.get_data(i,k));
                    }
                    exit(1);
                }
            }

        }
        if(n_neigh>max_n_neigh){
            max_n_neigh=n_neigh;
        }
        vol=1.0;
        for(k=0;k<dim;k++){
            vol*=local_max.get_data(k)-local_min.get_data(k);
        }
        if(vol<min_vol){
            min_vol=vol;
        }
        mean_vol+=vol;
        var_vol+=vol*vol;
        box_min.add_row(local_min);
        box_max.add_row(local_max);
        if(box_min.get_rows()%10000==0){
            mean_vol=mean_vol/10000.0;
            var_vol = var_vol/10000.0-mean_vol*mean_vol;
            printf("%d %e %d -- %e %e %e\n",
            box_min.get_rows(),double(time(NULL))-t_start,max_n_neigh,min_vol,
            mean_vol,sqrt(var_vol));
            min_vol=2.0*exception_value;
            mean_vol=0.0;
            var_vol=0.0;
        }
    }

    if(box_min.get_rows()!=dalex_chisq.get_dim()){
        printf("somehow got %d box_min and %d chisq\n",
        box_min.get_rows(), dalex_chisq.get_dim());
        exit(1);
    }

    if(box_max.get_rows()!=dalex_chisq.get_dim()){
        printf("somehow got %d box_max %d chisq\n",
        box_max.get_rows(),dalex_chisq.get_dim());
        exit(1);
    }

    FILE *out_file;
    out_file=fopen("volume_plot.txt", "w");
    double ln_vol;
    for(i=0;i<dalex_chisq.get_dim();i++){
        vol=1.0;
        ln_vol=0.0;
        for(j=0;j<dim;j++){
            vol*=box_max.get_data(i,j)-box_min.get_data(i,j);
            ln_vol+=log(box_max.get_data(i,j)-box_min.get_data(i,j));
        }
        fprintf(out_file,"%e %e %e\n",dalex_chisq.get_data(i),vol,ln_vol);
    }
    fclose(out_file);

    double dx;
    array_1d<double> log_volume;
    log_volume.set_name("log_volume");
    for(i=0;i<dalex_chisq.get_dim();i++){
        ln_vol=0.0;
        for(j=0;j<dim;j++){
            dx=box_max.get_data(i,j)-box_min.get_data(i,j);
            if(dx>1.0e-30){
                ln_vol += log(dx);
            }
            else{
                ln_vol=-69.0*dim;
                break;
            }
        }
        log_volume.add(ln_vol);
    }

    if(log_volume.get_dim()!=dalex_chisq.get_dim()){
        printf("somehow got %d volumes and %d chisq\n",
        log_volume.get_dim(),dalex_chisq.get_dim());
    }

    array_1d<double> log_posterior;
    log_posterior.set_name("log_posterior");
    array_1d<int> log_posterior_dex;
    log_posterior_dex.set_name("log_posterior_dex");
    for(i=0;i<dalex_chisq.get_dim();i++){
        log_posterior.set(i,-0.5*(dalex_chisq.get_data(i)-chisq_min)+log_volume.get_data(i));
        log_posterior_dex.set(i,i);
    }
    if(log_posterior.get_dim()!=dalex_chisq.get_dim()){
        printf("somehow got %d log_posteriors, but %d chisq\n",
        log_posterior.get_dim(),dalex_chisq.get_dim());
        exit(1);
    }
    array_1d<double> log_posterior_sorted;
    log_posterior_sorted.set_name("log_posterior_sorted");
    sort(log_posterior,log_posterior_sorted,log_posterior_dex);

    double total_posterior=0.0;
    for(i=0;i<log_posterior.get_dim();i++){
        j=log_posterior_dex.get_data(i);
        total_posterior+=exp(log_posterior.get_data(j));
    }

    array_1d<int> chisq_dex;
    chisq_dex.set_name("chisq_dex");
    array_1d<double> chisq_sorted;
    chisq_sorted.set_name("chisq_sorted");

    for(i=0;i<dalex_chisq.get_dim();i++){
        chisq_dex.set(i,i);
    }

    sort(dalex_chisq,chisq_sorted,chisq_dex);

    double local_prob=0.0;
    out_file = fopen(out_name, "w");
    fprintf(out_file,"# prob/total chisq log_vol min0 max0 min1 max1...\n");
    for(i=0;i<dalex_chisq.get_dim();i++){
        j=chisq_dex.get_data(i);
        local_prob+=exp(log_posterior.get_data(j));

        fprintf(out_file,"%le %le %le ",
        local_prob/total_posterior,dalex_chisq.get_data(j),log_volume.get_data(j));

        for(k=0;k<dim;k++){
            fprintf(out_file,"%le %le ",box_min.get_data(j,k),box_max.get_data(j,k));
        }
        fprintf(out_file,"\n");

    }
    fclose(out_file);

}
