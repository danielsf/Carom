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

double hyperbox::_median_metric(int i_dim, double xmid, array_2d<double> &pts){
    int i;
    int n1=0;
    for(i=0;i<pts.get_rows();i++){
        if(pts.get_data(i,i_dim)<xmid){
            n1++;
        }
    }

    if(n1==0 || n1==pts.get_rows()){
        return -1.0;
    }

    double ell=_max.get_data(i_dim)-_min.get_data(i_dim);
    return fabs(xmid-0.5*(_max.get_data(i_dim)+_min.get_data(i_dim)))/ell;
}

double hyperbox::_n_metric(int i_dim, double xmid, array_2d<double> &pts){
    int i;
    int n1=0;
    for(i=0;i<pts.get_rows();i++){
        if(pts.get_data(i,i_dim)<xmid){
            n1++;
        }
    }

    if(n1==0 || n1==pts.get_rows()){
        return -1.0;
    }

    return fabs(n1-pts.get_rows()/2.0);
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
    int dim_best;
    array_1d<double> local_chi;
    array_1d<double> local_chi_sorted;
    array_1d<int> local_chi_dex;
    local_chi.set_name("split_local_chi");
    local_chi_sorted.set_name("split_local_chi_sorted");
    local_chi_dex.set_name("split_local_chi_dex");
    double metric;
    double metric_best;
    double xx_best;
    double xx;
    dim_best=-1;
    int i;
    int i_mid;
    for(i=0;i<_pts.get_rows();i++){
        local_chi.add(_pts.get_data(i,dim()));
        local_chi_dex.set(i,i);
    }
    sort(local_chi,local_chi_sorted,local_chi_dex);
    i_mid=local_chi_dex.get_data(local_chi.get_dim()/2);
    for(i_dim=0;i_dim<dim();i_dim++){
        xx=_pts.get_data(i_mid,i_dim);
        metric=_median_metric(i_dim,xx,_pts);
        if(metric>-0.1){
            if(dim_best<0 || metric<metric_best){
                dim_best=i_dim;
                xx_best=xx;
            }
        }
    }

    if(dim_best>=0){
        _split_on_val(pts1,min1,max1,pts2,min2,max2,dim_best,xx_best);
        return;
    }

    for(i_dim=0;i_dim<dim();i_dim++){
        xx=0.5*(_pts.get_data(local_chi_dex.get_data(0),i_dim)+
                _pts.get_data(local_chi_dex.get_data(local_chi_dex.get_dim()-1),i_dim));
        metric=_n_metric(i_dim,xx,_pts);
        if(metric>-0.1){
            if(dim_best<0 || metric<metric_best){
                dim_best=i_dim;
                xx_best=xx;
            }
        }
    }

    if(dim_best>=0){
        _split_on_val(pts1,min1,max1,pts2,min2,max2,dim_best,xx_best);
        return;
    }

    printf("all splits failed %e %e %e %d\n",local_chi_sorted.get_data(0),
    local_chi_sorted.get_data(local_chi_sorted.get_dim()/2),
    local_chi_sorted.get_data(local_chi_sorted.get_dim()-1),
    local_chi_sorted.get_dim());
    exit(1);
}
