#include "ellipse.h"

int ellipse::contains(const array_1d<double> &pt){
    int i,j;
    double sum=0.0;
    double component;
    for(i=0;i<_bases.get_rows();i++){
        component=0.0;
        for(j=0;j<_bases.get_rows();j++){
            component+=(pt.get_data(j)-_center.get_data(j))*_bases.get_data(i,j);
        }
        sum+=power(component/_radii.get_data(i),2);
        if(sum>1.0){
            return 0;
        }
    }
    return 1;
}

void ellipse::build(array_2d<double> &pts_in){

    int n_pts = pts_in.get_rows();
    int dim = pts_in.get_cols();

    if(n_pts<dim){
        printf("WARNING cannot construct ellipse from %d pts in %d dimensions\n",
        n_pts,dim);
        exit(1);
    }

    _bases.reset_preserving_room();
    _radii.reset_preserving_room();

    array_1d<double> dir,dir_max,min,max;
    dir.set_name("ellipse_build_dir");
    dir_max.set_name("ellipse_build_dir_max");
    min.set_name("ellipse_build_min");
    max.set_name("ellipse_build_max");
    int i,j,k;
    _center.set_dim(dim);
    _center.zero();
    for(i=0;i<n_pts;i++){
        for(j=0;j<dim;j++){
            if(j>=min.get_dim() || pts_in.get_data(i,j)<min.get_data(j)){
                min.set(j,pts_in.get_data(i,j));
            }
            if(j>=max.get_dim() || pts_in.get_data(i,j)>max.get_data(j)){
                max.set(j,pts_in.get_data(i,j));
            }
        }
    }

    for(i=0;i<dim;i++){
        _center.set(i,0.5*(min.get_data(i)+max.get_data(i)));
    }

    printf("got center\n");

    double component,norm,norm_max;

    while(_bases.get_rows()!=dim){
        for(i=0;i<n_pts;i++){
            for(j=0;j<dim;j++){
                dir.set(j,pts_in.get_data(i,j)-_center.get_data(j));
            }
            for(j=0;j<_bases.get_rows();j++){
                component=0.0;
                for(k=0;k<dim;k++){
                    component+=dir.get_data(k)*_bases.get_data(j,k);
                }
                for(k=0;k<dim;k++){
                    dir.subtract_val(k,component*_bases.get_data(j,k));
                }
            }
            norm=dir.normalize();
            if(i==0 || norm>norm_max){
                norm_max=norm;
                for(j=0;j<dim;j++){
                    dir_max.set(j,dir.get_data(j));
                }
            }
        }
        if(norm_max<1.0e-20){
            printf("WARNING ellipse basis dir has norm %e\n",norm_max);
            exit(1);
        }
        _bases.add_row(dir_max);
        _radii.add(norm_max);
    }

    printf("got bases\n");

    int is_valid=0;

    while(is_valid==0){
        is_valid=1;
        for(i=0;i<n_pts;i++){
            if(contains(pts_in(i))==0){
                is_valid=0;
                break;
            }
        }

        if(is_valid==0){
            _set_radii(pts_in);
        }
    }

}


void ellipse::_set_radii(array_2d<double> &pts_in){

    array_1d<int> bad_dexes;
    int i;
    for(i=0;i<pts_in.get_rows();i++){
        if(contains(pts_in(i))==0){
            bad_dexes.add(i);
        }
    }

    int ipt;
    array_1d<double> norm,norm_sorted;
    array_1d<int> norm_dex,norm_ct;
    for(i=0;i<pts_in.get_cols();i++){
        norm_ct.set(i,0);
    }

    int j,k;
    double component;
    for(i=0;i<bad_dexes.get_dim();i++){
        ipt=bad_dexes.get_data(i);
        for(j=0;j<pts_in.get_cols();j++){
            component=0.0;
            for(k=0;k<pts_in.get_cols();k++){
                component+=(pts_in.get_data(ipt,k)-_center.get_data(k))*_bases.get_data(j,k);
            }
            norm.set(j,fabs(component)/_radii.get_data(j));
            norm_dex.set(j,j);
        }
        sort(norm,norm_sorted,norm_dex);
        j=norm_dex.get_data(pts_in.get_cols()-1);
        norm_ct.add_val(j,1);
    }

    norm_dex.reset_preserving_room();
    for(i=0;i<pts_in.get_cols();i++){
        norm_dex.set(i,i);
    }
    array_1d<int> norm_ct_sorted;
    sort(norm_ct,norm_ct_sorted,norm_dex);
    i=norm_dex.get_data(pts_in.get_cols()-1);

    _radii.multiply_val(i,1.1);

}
