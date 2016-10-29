#include "ellipse.h"

int ellipse::contains(array_1d<double> &pt){
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

    array_1d<double> bad_pt,pt,pt_proj,bad_pt_proj;
    array_1d<double> pt_norm,pt_norm_sorted;
    array_1d<int> pt_norm_dex,adjust_dexes;
    double max_remainder,local_remainder;
    int is_valid=0;

    bad_pt.set_name("ell_build_bad_pt");
    pt.set_name("ell_build_pt");
    pt_proj.set_name("ell_build_pt_proj");
    bad_pt_proj.set_name("ell_build_bad_pt_proj");
    pt_norm.set_name("ell_build_pt_norm");
    pt_norm_sorted.set_name("ell_build_pt_norm_sorted");
    pt_norm_dex.set_name("ell_build_pt_norm_dex");
    adjust_dexes.set_name("ell_build_adjust_dexes");

    while(is_valid==0){
        for(i=0;i<n_pts;i++){
            local_remainder=0.0;
            for(j=0;j<dim;j++){
                component=0.0;
                for(k=0;k<dim;k++){
                    component+=(pts_in.get_data(i,k)-_center.get_data(k))*_bases.get_data(j,k);
                }
                local_remainder+=power(component/_radii.get_data(j),2);
                pt.set(j,pts_in.get_data(i,j));
                pt_proj.set(j,fabs(component));
            }
            if(i==0 || local_remainder>max_remainder){
                max_remainder=local_remainder;
                for(j=0;j<dim;j++){
                    bad_pt.set(j,pt.get_data(j));
                    bad_pt_proj.set(j,pt_proj.get_data(j));
                }
            }
        }
        if(max_remainder<=1.0){
            is_valid=1;
        }
        else{
             adjust_dexes.reset_preserving_room();
             for(i=0;i<dim;i++){
                 pt_norm.set(i,bad_pt_proj.get_data(i)/_radii.get_data(i));
                 pt_norm_dex.set(i,i);
             }
             sort(pt_norm, pt_norm_sorted, pt_norm_dex);
             component=0.0;
             for(i=0;i<dim;i++){
                 component+=pt_norm_sorted.get_data(i);
                 if(component>0.5){
                     adjust_dexes.add(pt_norm_dex.get_data(i));
                 }
             }
             while(contains(bad_pt)==0){
                 for(i=0;i<adjust_dexes.get_dim();i++){
                     j=adjust_dexes.get_dim()-1-i;
                     _radii.multiply_val(adjust_dexes.get_data(j),1.1);
                     if(contains(bad_pt)==1){
                         break;
                     }
                 }
             }
        }
    }

}
