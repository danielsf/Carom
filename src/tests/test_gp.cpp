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
    exit(1);

    double guess,lin_guess,truth;

    int lin_better=0,gp_better=0;

    for(i=0;i<10000;i++){
        for(j=0;j<dim;j++){
            trial.set(j,min.get_data(j)+dice.doub()*(max.get_data(j)-min.get_data(j)));
        }
        truth=chifn(trial);
        guess=exp(gg(trial));
        lin_guess=exp(gg_lin(trial));
        gg.add_pt(trial,log(truth));
        gg_lin.add_pt(trial,log(truth));
        if(fabs(truth-guess)<fabs(truth-lin_guess)){
            gp_better++;
        }
        else{
            lin_better++;
        }
        printf("truth %e guess %e lin_guess %e\n",truth,guess,lin_guess);
        
    }

    printf("gp better %d lin better %d\n",gp_better, lin_better);

}
