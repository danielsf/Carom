#include "exampleLikelihoods.h"
#include "gp.h"
#include "gp_lin.h"

int main(){

    Ran dice(42);

    gaussianJellyBean24 chifn;

    array_2d<double> points;
    array_1d<double> fn,max,min,gmax,gmin,trial;
    int dim=24;
    points.set_cols(dim);

    int i,j;
    for(i=0;i<dim;i++){
        min.set(i,-40.0);
        max.set(i,40.0);
        gmin.set(i,0.0);
        gmax.set(i,1.0);
    }
    for(i=0;i<10000;i++){
        for(j=0;j<dim;j++){
            trial.set(j,min.get_data(j)+dice.doub()*(max.get_data(j)-min.get_data(j)));
        }
        fn.add(chifn(trial));
        points.add_row(trial);
    }

    gp gg;
    gp_lin gg_lin;
    gg.build(points,fn,min,max);
    gg_lin.build(points,fn,min,max);

    array_2d<double> test_set;
    array_1d<double> test_ff;
    test_set.set_cols(dim);
    for(i=0;i<3000;i++){
        for(j=0;j<dim;j++){
            trial.set(j,min.get_data(j)+dice.doub()*(max.get_data(j)-min.get_data(j)));
        }
        test_ff.add(chifn(trial));
        test_set.add_row(trial);
    }
    gg.optimize(test_set,test_ff);
    printf("\nnow lin\n");
    gg_lin.optimize(test_set,test_ff);

    double gp_guess,lin_guess,truth;

    int lin_better=0,gp_better=0;
    double lin_bias=0.0,gp_bias=0.0;
    double lin_norm_bias=0.0,gp_norm_bias=0.0;
    double lin_sq=0.0,gp_sq=0.0;


    for(i=0;i<10000;i++){
        for(j=0;j<dim;j++){
            trial.set(j,min.get_data(j)+dice.doub()*(max.get_data(j)-min.get_data(j)));
        }
        truth=chifn(trial);
        gp_guess=gg(trial);
        lin_guess=gg_lin(trial);
        //gg.add_pt(trial,truth);
        //gg_lin.add_pt(trial,truth);
        if(fabs(truth-gp_guess)<fabs(truth-lin_guess)){
            gp_better++;
        }
        else{
            lin_better++;
        }
        lin_bias+=lin_guess-truth;
        gp_bias+=gp_guess-truth;
        lin_norm_bias+=(lin_guess-truth)/fabs(truth);
        gp_norm_bias+=(gp_guess-truth)/fabs(truth);
        lin_sq+=power(lin_guess-truth,2);
        gp_sq+=power(gp_guess-truth,2);
    }

    gp_bias=gp_bias/double(i);
    lin_bias=lin_bias/double(i);
    gp_norm_bias=gp_norm_bias/double(i);
    lin_norm_bias=lin_norm_bias/double(i);
    gp_sq=gp_sq/double(i);
    lin_sq=gp_sq/double(i);

    printf("gp better %d lin better %d\n",gp_better, lin_better);
    printf("gp bias %e lin bias %e\n",gp_bias,lin_bias);
    printf("gp norm bias %e lin norm bias %e\n",gp_norm_bias, lin_norm_bias);
    printf("gp var %e lin var %e\n",gp_sq-gp_bias*gp_bias, lin_sq-lin_bias*lin_bias);

}
