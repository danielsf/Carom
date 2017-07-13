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

    int i;
    int i_dim;
    double xmid;
    int dim_best=-1;
    double metric;
    double metric_best;

    for(i_dim=0;i_dim<dim();i_dim++){
        xmid=0.5*(_max.get_data(i_dim)+_min.get_data(i_dim));
        metric=_n_metric(i_dim,xmid,_pts);
        if(metric>-0.1){
            if(dim_best<0 || metric<metric_best){
                dim_best=i_dim;
                metric_best=metric;
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
    array_1d<double> pt_xmid;
    pt_xmid.set_name("hypberbox_split_pt_xmid");

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
        pt_xmid.set(i_dim,0.5*(x_val_sorted.get_data(0)+x_val_sorted.get_data(sorted_dex.get_dim()-1)));

        metric=_n_metric(i_dim,xmid,_pts);
        if(metric>-0.1){
            if(dim_best<0 || metric<metric_best){
                dim_best=i_dim;
                xmid_best=xmid;
                metric_best=metric;
            }
        }
    }

    if(dim_best>=0){
        _split_on_val(pts1,min1,max1,pts2,min2,max2,dim_best,xmid_best);
        return;
    }

    for(i_dim=0;i_dim<dim();i_dim++){
        metric=_n_metric(i_dim,pt_xmid.get_data(i_dim),_pts);
        if(metric>-0.1){
            if(dim_best<0 || metric<metric_best){
                dim_best=i_dim;
                metric_best=metric;
            }
        }
    }

    if(dim_best<0){
        printf("exhausted all splitting searches\n");
        exit(1);
    }
    _split_on_val(pts1,min1,max1,pts2,min2,max2,dim_best,pt_xmid.get_data(dim_best));

}
