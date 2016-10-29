#include "containers.h"
#include "goto_tools.h"
#include "ellipse.h"

int main(){

    array_1d<double> center;
    array_2d<double> bases;
    array_2d<double> pts;
    Ran dice(99);
    int dim=4;
    int n_pts=1000;

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
    sampler.initialize(dim, 112);
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
    
    


}
