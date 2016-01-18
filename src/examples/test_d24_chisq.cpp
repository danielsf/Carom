#include "exampleLikelihoods.h"

int main(){

    gaussianJellyBean24 chisq;

    int dim=chisq.get_dim();

    array_1d<double> pt_projected,params,pt;
    int i,j;
    for(i=0;i<dim;i++){

            pt_projected.set(i,chisq.get_center(0,i));
        
    }


    double chimin,mu;
    array_1d<double> min_pt;

    chimin=2.0*exception_value;

    double xx,yy;
    for(xx=-40.0;xx<40.0;xx+=0.04){
        for(yy=-40.0;yy<40.0;yy+=0.04){
            pt_projected.set(0,xx);
            pt_projected.set(0,yy);
            for(i=0;i<dim;i++){
                pt.set(i,0.0);
            }
            for(i=0;i<dim;i++){
                for(j=0;j<dim;j++){
                    pt.add_val(i,pt_projected.get_data(j)*chisq.get_basis(j,i));
                }
            }
            
            chisq.convert_params(pt, params, 0);
            for(i=2;i<dim;i++){
                if(i%4==2){
                    if(fabs(params.get_data(i)-1.0)>0.1){
                        printf("pt_projected %e\n",pt_projected.get_data(i));
                        printf("WARNING pp %d %e\n",i,params.get_data(i));
                        exit(1);
                    }
                }
                else{
                    if(fabs(params.get_data(i))>0.1){
                        printf("WARNING pp %d %e\n",i,params.get_data(i));
                    }
                }
            }
            mu=chisq(pt);
            if(mu<chimin){
                chimin=mu;
                for(i=0;i<dim;i++){
                    min_pt.set(i,pt.get_data(i));
                }
                printf("chimin %e\n",chimin);
                for(i=0;i<dim;i++){
                    printf("    %.3e\n",min_pt.get_data(i));
                }
            }
        }
    }

    printf("min pt -- %e\n",chimin);
    for(i=0;i<dim;i++){
        printf("%e\n",min_pt.get_data(i));
    }

}
