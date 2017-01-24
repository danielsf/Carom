#include "ellipse.h"


void recursive_ellipse(array_2d<double> &pts_in, ellipse_list &list_out){

    ellipse singleton;
    singleton.build(pts_in);

    if(pts_in.get_rows()<2*pts_in.get_cols()){
        list_out.add(singleton);
        return;
    }

    array_2d<double> first_half,second_half;
    first_half.set_name("first_half");
    second_half.set_name("second_half");

    int i;
    for(i=0;i<pts_in.get_rows()/2;i++){
        first_half.add_row(pts_in(i));
    }

    for(;i<pts_in.get_rows();i++){
        second_half.add_row(pts_in(i));
    }

    ellipse first_ellipse;
    ellipse second_ellipse;

    first_ellipse.build(first_half);
    second_ellipse.build(second_half);

    double vol_0,vol_1,vol_2;
    vol_0=1.0;
    vol_1=1.0;
    vol_2=1.0;
    for(i=0;i<pts_in.get_cols();i++){
        vol_0*=singleton.radii(i);
        vol_1*=first_ellipse.radii(i);
        vol_2*=second_ellipse.radii(i);
    }

    printf("volumes %e %e %e %e\n",vol_0,vol_1,vol_2,vol_1+vol_2);

    double n_pts=double(pts_in.get_rows());
    double dim = double(pts_in.get_cols());

    double bic_0=log(vol_0);
    double bic_1=log(vol_1+vol_2)+(2)*log(dim);

    printf("bic %e %e -- %e %e %e\n",bic_0,bic_1,log(vol_1+vol_2), log(n_pts), log(vol_0));

    if(bic_0<bic_1){
        list_out.add(singleton);
        return;
    }

    first_half.reset();
    second_half.reset();

    for(i=0;i<pts_in.get_rows()/2;i++){
        first_half.add_row(pts_in(i));
    }
    recursive_ellipse(first_half, list_out);

    first_half.reset();

    for(;i<pts_in.get_rows();i++){
        first_half.add_row(pts_in(i));
    }
    recursive_ellipse(first_half, list_out);

}


int main(int iargc, char **argv){

    Ran dice(88);
    array_1d<double> dir,epsilon;
    dir.set_name("dir");
    epsilon.set_name("epsilon");
    
    int dim=5;
    int i;
    for(i=0;i<dim;i++){
        dir.set(i,normal_deviate(&dice,0.0,1.0));
    }
    dir.normalize();
    
    array_1d<double> mean,trial;
    array_2d<double> pts;
    mean.set_name("mean");
    trial.set_name("trial");
    pts.set_name("pts");
    
    for(i=0;i<dim;i++){
        mean.set(i,0.0);
    }
    
    int j;
    for(i=0;i<1000;i++){
        if(i==250){
            for(j=0;j<dim;j++){
                dir.set(j,normal_deviate(&dice,0.0,1.0));
            }
            dir.normalize();
        }
        for(j=0;j<dim;j++){
            epsilon.set(j,normal_deviate(&dice,0.0,1.0));
        }
        epsilon.normalize();
        for(j=0;j<dim;j++){
            trial.set(j,mean.get_data(j)+1.0*epsilon.get_data(j));
        }
        pts.add_row(trial);
        for(j=0;j<dim;j++){
            mean.add_val(j,dir.get_data(j)*0.1);
        }
    }

    FILE *output;
    output=fopen("pts_test.txt", "w");
    for(i=0;i<pts.get_rows();i++){
        for(j=0;j<pts.get_cols();j++){
            fprintf(output,"%e ",pts.get_data(i,j));
        }
        fprintf(output,"\n");
    }
    fclose(output);
    ellipse_list e_list;

    recursive_ellipse(pts, e_list);
    printf("ellipses: %d\n",e_list.ct());

    double component;
    int ix;
    char out_name[100];
    double sign;
    for(ix=0;ix<e_list.ct();ix++){
        sprintf(out_name,"ellipse_%d.txt",ix);
        output=fopen(out_name,"w");
        for(i=0;i<dim;i++){
            for(sign=-1.0;sign<1.5;sign+=2.0){
                for(j=0;j<dim;j++){
                    fprintf(output,"%e ",e_list(ix)->center(j)+sign*e_list(ix)->bases(i,j)*e_list(ix)->radii(i));
                }
                fprintf(output,"\n");
            }
        }
        
        fclose(output);
    }

}
