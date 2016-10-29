#include "containers.h"
#include "goto_tools.h"
#include "ellipse.h"

int main(){

    int dim=5;
    int n_pts=10000;
    int seed1=99;
    int seed2=342;

    array_1d<double> center;
    array_2d<double> bases;
    array_2d<double> pts;
    Ran dice(seed1);

    int i,j,k;
    array_1d<double> dir;
    double component,norm;
    while(bases.get_rows()!=dim){
        for(i=0;i<dim;i++){
            dir.set(i,normal_deviate(&dice, 0.0, 1.0));
        }
        for(i=0;i<bases.get_rows();i++){
            component=0.0;
            for(j=0;j<dim;j++){
                component+=dir.get_data(j)*bases.get_data(i,j);
            }
            for(j=0;j<dim;j++){
                dir.subtract_val(j,component*bases.get_data(i,j));
            }
        }
        norm=dir.normalize();
        if(norm>1.0e-20){
            bases.add_row(dir);
        }
    }

    array_1d<double> radii;
    for(i=0;i<dim;i++){
        center.set(i,normal_deviate(&dice, 0.0, 10.0));
        radii.set(i,0.1+dice.doub()*2.0);
    }

    ellipse_sampler sampler;
    sampler.initialize(dim, seed2);
    array_1d<double> base_pt,true_pt;
    for(i=0;i<n_pts;i++){
        sampler.get_pt(base_pt);
        for(j=0;j<dim;j++){
            true_pt.set(j,base_pt.get_data(j)*radii.get_data(j)+center.get_data(j));
        }
        pts.add_row(true_pt);
    }

    ellipse ell;
    ell.build(pts);
    
    double dd=0.0;
    for(i=0;i<dim;i++){
        dd+=power(ell.center(i)-center.get_data(i),2);
    }
    printf("center distance %e\n",sqrt(dd));
    
    for(i=0;i<n_pts;i++){
        if(ell.contains(pts(i))==0){
            printf("WARNING does not contain a pt\n");
            exit(1);
        }
    }


    double dot;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            dot=0.0;
            for(k=0;k<dim;k++){
                dot+=bases.get_data(i,k)*ell.bases(j,k);
            }
            printf("%d %d -- %e -- %e %e\n",
            i,j,dot,radii.get_data(i),ell.radii(j));
        }
        printf("\n");
    }

}
