#include "exampleLikelihoods.h"
#include "gp.h"

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
        fn.add(log(chifn(trial)));
        points.add_row(trial);
    }

    gp gg;
    gg.build(points,fn,min,max);

    double guess,truth;

    for(i=0;i<10;i++){
        for(j=0;j<dim;j++){
            trial.set(j,min.get_data(j)+dice.doub()*(max.get_data(j)-min.get_data(j)));
        }
        truth=chifn(trial);
        guess=gg(trial);
        gg.add_pt(trial,log(truth));
        printf("truth %e guess %e\n",log(truth),guess);
        
    }

}
