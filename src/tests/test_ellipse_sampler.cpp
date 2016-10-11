#include "goto_tools.h"

int main(){

    Ran dice(99.0);

    int dim=4;
    int i_dim;
    array_1d<double> radii;
    for(i_dim=0;i_dim<dim;i_dim++){
        radii.set(i_dim,dice.doub()*2.0);
    }

    int n_samples=10000;
    int i_sample;
    array_1d<double> pt;
    ellipse_sampler sampler(dim,3215);
    double remainder;
    FILE *output;
    output=fopen("ellipse_samples.txt", "w");
    for(i_sample=0;i_sample<n_samples;i_sample++){
        sampler.get_pt(pt);
        for(i_dim=0;i_dim<dim;i_dim++){
            pt.multiply_val(i_dim,radii.get_data(i_dim));
            fprintf(output,"%e ",pt.get_data(i_dim));
        }
        fprintf(output," %e\n",remainder);
        remainder=0.0;
        for(i_dim=0;i_dim<dim;i_dim++){
            remainder+=power(pt.get_data(i_dim)/radii.get_data(i_dim),2);
        }
        remainder=sqrt(remainder);
        if(remainder>1.01){
            printf("WARNING remainder %e\n",remainder);
            exit(1);
        }
    }
    fclose(output);

    for(i_dim=0;i_dim<dim;i_dim++){
      printf("radius %d %e\n",i_dim,radii.get_data(i_dim));
    }

}
