#include "hyperbox.h"

void hyperbox::_split_on_val(array_2d<double> &pts1,
                             array_1d<double> &min1,
                             array_1d<double> &max1,
                             array_2d<double> &pts2,
                             array_1d<double> &min2,
                             array_1d<double> &max2,
                             int dim_split, double x_split){

    pts1.reset_preserving_room();
    min1.reset_preserving_room();
    max1.reset_preserving_room();
    pts2.reset_preserving_room();
    min2.reset_preserving_room();
    max2.reset_preserving_room();

    int i;
    for(i=0;i<dim();i++){
        min1.set(i,_min.get_data(i));
        min2.set(i,_min.get_data(i));
        max1.set(i,_max.get_data(i));
        max2.set(i,_max.get_data(i));
    }

    max1.set(dim_split,x_split);
    min2.set(dim_split,x_split);

    for(i=0;i<n_pts();i++){
        if(_pts.get_data(i,dim_split)<x_split){
            pts1.add_row(_pts(i));
        }
        else{
            pts2.add_row(_pts(i));
        }
    }

    if(pts1.get_rows()==0 || pts2.get_rows()==0){
        printf("hyperbox::_split_on_val gave bad division %d %d\n",
        pts1.get_rows(),pts2.get_rows());
        exit(1);
    }

    for(i=0;i<dim();i++){
        if(min1.get_data(i)>max1.get_data(i)){
            printf("hyperbox::_split_on_val gave bad min1 %e max1 %e\n",
            min1.get_data(i),max1.get_data(i));
            exit(1);
        }
        if(min2.get_data(i)>max2.get_data(i)){
            printf("hyperbox::_split_on_val gave bad min2 %e max2 %e\n",
            min2.get_data(i),max2.get_data(i));
            exit(1);
        }
    }

    for(i=0;i<pts1.get_rows();i++){
        if(pts1.get_data(i,dim_split)<min1.get_data(dim_split) ||
           pts1.get_data(i,dim_split)>max1.get_data(dim_split)){

            printf("hyperbox::_split_on_val failed %e %e %e\n",
            min1.get_data(dim_split),
            pts1.get_data(i,dim_split),
            max1.get_data(dim_split));

            exit(1);
        }
    }

    for(i=0;i<pts2.get_rows();i++){
        if(pts2.get_data(i,dim_split)<min2.get_data(dim_split) ||
           pts2.get_data(i,dim_split)>max2.get_data(dim_split)){

            printf("hyperbox::_split_on_val failed %e %e %e\n",
            min2.get_data(dim_split),
            pts2.get_data(i,dim_split),
            max2.get_data(dim_split));

            exit(1);
        }
    }
}

void hyperbox::split(array_2d<double> &pts1,
                     array_1d<double> &min1,
                     array_1d<double> &max1,
                     array_2d<double> &pts2,
                     array_1d<double> &min2,
                     array_1d<double> &max2){


    if(_pts.get_rows()<2){
        printf("trying to split a hyperbox with %d points\n",
        _pts.get_rows());
        exit(1);
    }

    int i_dim;
    array_1d<int> n_lower,n_upper;
    n_lower.set_name("hyperbox_split_n_lower");
    n_upper.set_name("hyperbox_split_n_upper");
    int i;
    double xmid;
    n_lower.set_dim(dim());
    n_upper.set_dim(dim());
    for(i_dim=0;i_dim<dim();i_dim++){
        n_lower.set(i_dim,0);
        n_upper.set(i_dim,0);
        xmid=0.5*(_min.get_data(i_dim)+_max.get_data(i_dim));
        for(i=0;i<_pts.get_rows();i++){
            if(_pts.get_data(i,i_dim)<xmid){
                n_lower.add_val(i_dim,1);
            }
            else{
                n_upper.add_val(i_dim,1);
            }
        }
    }

    int n_best;
    int dim_best=-1;
    for(i_dim=0;i_dim<dim();i_dim++){
        if(n_lower.get_dim()!=0 && n_lower.get_dim()!=n_pts()){
            if(dim_best<0 ||
               abs(n_lower.get_data(i_dim)-n_pts()/2)<abs(n_best-n_pts()/2)){

                dim_best=i_dim;
                n_best=n_lower.get_data(i_dim);
            }
        }
    }

    if(dim_best>=0){
        _split_on_val(pts1,min1,max1,pts2,min2,max2,dim_best,
                      0.5*(_min.get_data(dim_best)+_max.get_data(dim_best)));

        return;
    }

    array_1d<int> sorted_dex;
    array_1d<double> x_val,x_val_sorted;
    sorted_dex.set_name("hyperbox_split_sorted_dex");
    x_val.set_name("hyperbox_split_x_val");
    x_val_sorted.set_name("hyperbox_split_x_val_sorted");

    x_val.set_dim(n_pts());
    x_val_sorted.set_dim(n_pts());
    sorted_dex.set_dim(n_pts());

    double xmid_best;
    double dist_best;
    double dist;
    int n1;

    for(i_dim=0;i_dim<dim();i_dim++){
        x_val.reset_preserving_room();
        sorted_dex.reset_preserving_room();
        x_val_sorted.reset_preserving_room();
        for(i=0;i<n_pts();i++){
            x_val.set(i,_pts.get_data(i,i_dim));
            sorted_dex.set(i,i);
        }
        sort(x_val,x_val_sorted,sorted_dex);

        xmid=x_val_sorted.get_data(sorted_dex.get_dim()/2);

        n1=0;
        for(i=0;i<n_pts();i++){
            if(_pts.get_data(i,i_dim)<xmid){
                n1++;
            }
        }

        if(n1==n_pts() || n1==0){
            continue;
        }

        dist=fabs(1.0-xmid/0.5*(_min.get_data(i_dim)+_max.get_data(i_dim)));
        if(dim_best<0 || dist<dist_best){
            xmid_best=xmid;
            dim_best=i_dim;
        }
    }

    if(dim_best<0){
       printf("splitting on the median did not help in hypberbox:split");
       exit(1);
    }

    _split_on_val(pts1,min1,max1,pts2,min2,max2,dim_best,xmid_best);

}
