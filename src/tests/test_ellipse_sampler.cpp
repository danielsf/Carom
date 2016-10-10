#include "goto_tools.h"

int main(){

    Ran dice(99);

    array_2d<double> cos_n_grid;
    int dim=4;
    int steps=1000;
    cos_n_grid.set_dim(dim+1,steps);
    int i_dim,i_step;
    double d_x=1.0/double(steps);
    double theta;
    for(i_dim=0;i_dim<cos_n_grid.get_rows();i_dim++){
        for(i_step=0;i_step<steps;i_step++){
            theta=asin(i_step*d_x);
            cos_n_grid.set(i_dim,i_step,integrate_cos_n(0.0,theta,i_dim));
            /*if(i_dim==4 && i_step>0){
                printf("%d %d %e\n",i_dim,i_step,
                cos_n_grid.get_data(i_dim,i_step)-cos_n_grid.get_data(i_dim,i_step-1));
            }*/
        }
    }

    array_1d<double> radii;
    for(i_dim=0;i_dim<dim;i_dim++){
        radii.set(i_dim,dice.doub()*2.0);
    }

    int n_samples=500;
    int i_sample;
    int i_max_vol;
    int dim_used;
    array_1d<int> dim_record;
    double remainder;
    double sqrt_remainder;
    double coord,roll,max_vol,cos_term;
    int dim_left;
    array_1d<double> pt;
    FILE *output;
    output=fopen("ellipse_samples.txt", "w");
    for(i_sample=0;i_sample<n_samples;i_sample++){
        remainder=1.0;
        sqrt_remainder=1.0;
        dim_record.reset_preserving_room();
        for(dim_used=0;dim_used<dim;dim_used++){

            dim_left=dim-dim_used-1;
            i_dim=-1;
            while(i_dim<0 || dim_record.contains(i_dim)==1){
                i_dim=dice.int32()%dim;
            }
            dim_record.add(i_dim);

            i_max_vol=int(sqrt_remainder/d_x);
            if(i_max_vol==cos_n_grid.get_cols()){
                i_max_vol--;
            }

            max_vol=cos_n_grid.get_data(dim_left,i_max_vol);
            roll=2.0*(dice.doub()-0.5)*max_vol;
            cos_term=0.0;
            for(i_step=0;i_step<cos_n_grid.get_cols()-1 && cos_term<fabs(roll);i_step++){
                cos_term=cos_n_grid.get_data(dim_left,i_step);
            }
            if(i_step!=0){
               i_step--;
            }
            coord=i_step*d_x;
            if(roll<0.0){
               coord*=-1.0;
            }
            /*if(i_dim==0){
                printf("roll %e i_step %d coord %e\n",roll,i_step,coord);
            }*/
            pt.set(i_dim,coord*radii.get_data(i_dim));
            remainder-=coord*coord;
            if(remainder<0.0){
                printf("WARNING remainder %e %d %e %e\n",remainder,i_dim,max_vol,coord);
                exit(1);
            }
            sqrt_remainder=sqrt(remainder);
        }
        for(i_dim=0;i_dim<dim;i_dim++){
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
