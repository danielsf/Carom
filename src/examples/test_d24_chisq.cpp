#include "exampleLikelihoods.h"

int main(){

    gaussianJellyBean12 chisq;

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
    for(xx=-40.0;xx<40.0;xx+=0.1){
        for(yy=-40.0;yy<40.0;yy+=0.1){
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

    min_pt.set(0,1.585e+00);
    min_pt.set(1,1.588e+01);
    min_pt.set(2,-5.004e+00);
    min_pt.set(3,1.515e+01);
    min_pt.set(4,-7.553e+00);
    min_pt.set(5,-4.807e+00);
    min_pt.set(6,3.819e+00);
    min_pt.set(7,2.069e+01);
    min_pt.set(8,-7.932e+00);
    min_pt.set(9,8.844e+00);
    min_pt.set(10,6.206e+00);
    min_pt.set(11,7.611e-01);

    printf("chisq at by hand min %e\n",chisq(min_pt));

    char word[100];
    double chisq_val;
    double junk;
    FILE *input;
    input=fopen("../Multinest_v3.9/chains/gaussianJellyBean_d12_s99_n300_carom.sav","r");
    for(i=0;i<dim+5;i++){
        fscanf(input,"%s",word);
    }
    printf("word %s\n",word);

    array_1d<double> trial;
    trial.set_name("trial");
    int n_connected=0;
    int n_disconnected=0;
    double target=21.0+chisq(min_pt);
    while(fscanf(input,"%le",&mu)==1){
        pt.set(0,mu);
        for(i=1;i<dim;i++){
            fscanf(input,"%le",&mu);
            pt.set(i,mu);
        }
        fscanf(input,"%le",&chisq_val);
        for(i=0;i<3;i++){
            fscanf(input,"%le",&mu);
        }
        if(chisq_val<target){
            for(i=0;i<dim;i++){
                trial.set(i,0.5*(pt.get_data(i)+min_pt.get_data(i)));
            }
            mu=chisq(trial);
            if(mu>target){
                n_disconnected++;
            }
            else{
                n_connected++;
                if(n_connected%1000==0){
                    printf("n_connected %d dis %d\n",n_connected,n_disconnected);
                }
            }
        }
    }

    fclose(input);
    printf("n_connected %d n_dis %d\n",n_connected,n_disconnected);

}
