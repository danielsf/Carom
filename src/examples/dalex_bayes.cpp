#include "containers.h"
#include "goto_tools.h"
#include "kd.h"

int main(int iargc, char *argv[]){

    int i,j,dim;
    char in_name[letters];
    char out_name[letters];
    array_1d<int> xdexes,ydexes;
    xdexes.set_name("xdexes");
    ydexes.set_name("ydexes");
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
                case 'x':
                    i++;
                    xdexes.add(atoi(argv[i]));
                    i++;
                    ydexes.add(atoi(argv[i]));
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
    if(xdexes.get_dim()==0){
        printf("need to specify dimensions to plot\n");
        exit(1);
    }
    if(xdexes.get_dim()!=ydexes.get_dim()){
        printf("somehow got %d xes but %d ys\n",
        xdexes.get_dim(),ydexes.get_dim());
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

    int n_neigh=4*dim+1;
    int mins_set;
    int maxes_set;
    int k;
    int dim_dex;
    int neigh_dex;
    int is_valid;
    int n_pass;

    array_1d<double> sorted_delta;
    array_1d<int> sorted_delta_dex;
    sorted_delta.set_name("sorted_delta");
    sorted_delta_dex.set_name("sorted_delta_dex");

    double t_start = double(time(NULL));

    for(i=0;i<dalex_pts.get_rows();i++){
        dalex_tree.nn_srch(i,n_neigh,neigh,dist);
        if(neigh.get_data(0)!=i){
            printf("neighbor search did not find self\n");
            exit(1);
        }

        for(j=0;j<dim;j++){
            local_min_set.set(j,0);
            local_max_set.set(j,0);
            local_min.set(j,dalex_pts.get_data(i,j));
            local_max.set(j,dalex_pts.get_data(i,j));
        }
        mins_set=0;
        maxes_set=0;
        is_valid=0;
        n_pass=0;
        while(is_valid==0 && n_pass<2){
            for(j=1;j<n_neigh && mins_set<dim && maxes_set<dim; j++){
                neigh_dex = neigh.get_data(j);
                sorted_delta.reset_preserving_room();
                sorted_delta_dex.reset_preserving_room();
                for(k=0;k<dim;k++){
                    delta.set(k,fabs(dalex_pts.get_data(i,k)-dalex_pts.get_data(neigh_dex,k)));
                    sorted_delta_dex.set(k,k);
                }
                sort(delta,sorted_delta,sorted_delta_dex);

                for(k=dim-1;k>=0;k--){
                    dim_dex=sorted_delta_dex.get_data(k);
                    if(dalex_pts.get_data(neigh_dex,dim_dex)>dalex_pts.get_data(i,dim_dex)
                       && local_max_set.get_data(dim_dex)==0){

                        local_max.set(dim_dex, dalex_pts.get_data(neigh_dex,dim_dex));
                        maxes_set++;
                        local_max_set.set(dim_dex,1);
                        if(n_pass==0){
                            break;
                        }
                    }
                    else if(dalex_pts.get_data(neigh_dex,dim_dex)<dalex_pts.get_data(i,dim_dex)
                            && local_min_set.get_data(dim_dex)==0){

                        local_min.set(dim_dex, dalex_pts.get_data(neigh_dex,dim_dex));
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
                    break;
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
        }
        box_min.add_row(local_min);
        box_max.add_row(local_max);
        if(box_min.get_rows()%1000==0){
            printf("%d %e\n",box_min.get_rows(),double(time(NULL))-t_start);
        }
    }
}
