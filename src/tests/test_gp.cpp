#include "exampleLikelihoods.h"
#include "gp.h"
#include "gp_lin.h"

int main(){

    Ran dice(42);

    gaussianJellyBean12 chifn;

    array_2d<double> points;
    array_1d<double> fn,max,min,gmax,gmin,trial;
    int dim=chifn.get_dim();
    points.set_cols(dim);

    int i,j;
    for(i=0;i<dim;i++){
        min.set(i,-40.0);
        max.set(i,40.0);
        gmin.set(i,0.0);
        gmax.set(i,1.0);
    }
    for(i=0;i<5000;i++){
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

    printf("\nbuilt\n");

    array_2d<double> test_set;
    array_1d<double> test_ff;
    test_set.set_cols(dim);
    for(i=0;i<5000;i++){
        for(j=0;j<dim;j++){
            trial.set(j,min.get_data(j)+dice.doub()*(max.get_data(j)-min.get_data(j)));
        }
        test_ff.add(chifn(trial));
        test_set.add_row(trial);
    }
    //gg.optimize(test_set,test_ff);
    printf("\nnow lin\n");
    gg_lin.optimize(test_set,test_ff);
    gg.set_ell_factor(5.623414e4);
    //gg_lin.set_ell_factor(10.0);

    double gp_guess,lin_guess,truth;

    int lin_better=0,gp_better=0;
    double lin_bias=0.0,gp_bias=0.0;
    double lin_norm_bias=0.0,gp_norm_bias=0.0;
    double lin_sq=0.0,gp_sq=0.0;

    array_1d<double> gp_norm_bias_list, gp_fabs_bias_list;
    array_1d<double> lin_norm_bias_list, lin_fabs_bias_list;
    array_1d<double> gp_norm_bias_sorted, gp_fabs_bias_sorted;
    array_1d<double> lin_norm_bias_sorted, lin_fabs_bias_sorted;
    array_1d<int> gp_norm_bias_dex, gp_fabs_bias_dex;
    array_1d<int> lin_norm_bias_dex, lin_fabs_bias_dex;

    double gp_norm_bias_term,lin_norm_bias_term;

    double t_0,t_gp,t_lin;
    t_gp=0.0;
    t_lin=0.0;

    for(i=0;i<10000;i++){
        if(i%1000==0){
            printf("    ii %d\n",i);
        }
        for(j=0;j<dim;j++){
            trial.set(j,min.get_data(j)+dice.doub()*(max.get_data(j)-min.get_data(j)));
        }
        truth=chifn(trial);

        t_0=double(time(NULL));
        gp_guess=gg(trial);
        t_gp+=double(time(NULL))-t_0;

        t_0=double(time(NULL));
        lin_guess=gg_lin(trial);
        t_lin+=double(time(NULL))-t_0;

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

        gp_norm_bias_term=(gp_guess-truth)/fabs(truth);
        lin_norm_bias_term=(lin_guess-truth)/fabs(truth);

        lin_norm_bias+=lin_norm_bias_term;
        gp_norm_bias+=gp_norm_bias_term;

        lin_sq+=power(lin_guess-truth,2);
        gp_sq+=power(gp_guess-truth,2);

        gp_norm_bias_list.add(gp_norm_bias_term);
        gp_fabs_bias_list.add(fabs(gp_norm_bias_term));
        gp_norm_bias_dex.add(i);
        gp_fabs_bias_dex.add(i);

        lin_norm_bias_list.add(lin_norm_bias_term);
        lin_fabs_bias_list.add(fabs(lin_norm_bias_term));
        lin_norm_bias_dex.add(i);
        lin_fabs_bias_dex.add(i);
    }

    gp_bias=gp_bias/double(i);
    lin_bias=lin_bias/double(i);
    gp_norm_bias=gp_norm_bias/double(i);
    lin_norm_bias=lin_norm_bias/double(i);
    gp_sq=gp_sq/double(i);
    lin_sq=lin_sq/double(i);

    printf("gp better %d lin better %d\n",gp_better, lin_better);
    printf("gp bias %e lin bias %e\n",gp_bias,lin_bias);
    printf("gp norm bias %e lin norm bias %e\n",gp_norm_bias, lin_norm_bias);
    printf("gp var %e lin var %e\n",gp_sq-gp_bias*gp_bias, lin_sq-lin_bias*lin_bias);

    sort_and_check(gp_norm_bias_list, gp_norm_bias_sorted, gp_norm_bias_dex);
    sort_and_check(gp_fabs_bias_list, gp_fabs_bias_sorted, gp_fabs_bias_dex);
    sort_and_check(lin_norm_bias_list, lin_norm_bias_sorted, lin_norm_bias_dex);
    sort_and_check(lin_fabs_bias_list, lin_fabs_bias_sorted, lin_fabs_bias_dex);

    i=lin_fabs_bias_list.get_dim();

    printf("\n");

    printf("quartiles gp norm %e %e %e\n",
    gp_norm_bias_sorted.get_data(i/4),
    gp_norm_bias_sorted.get_data(i/2),
    gp_norm_bias_sorted.get_data((3*i)/4));

    printf("quartiles lin norm %e %e %e\n",
    lin_norm_bias_sorted.get_data(i/4),
    lin_norm_bias_sorted.get_data(i/2),
    lin_norm_bias_sorted.get_data((3*i)/4));

    printf("\n");

    printf("quartiles gp fabs %e %e %e\n",
    gp_fabs_bias_sorted.get_data(i/4),
    gp_fabs_bias_sorted.get_data(i/2),
    gp_fabs_bias_sorted.get_data((3*i)/4));

    printf("quartiles lin fabs %e %e %e\n",
    lin_fabs_bias_sorted.get_data(i/4),
    lin_fabs_bias_sorted.get_data(i/2),
    lin_fabs_bias_sorted.get_data((3*i)/4));

    printf("\n");
    printf("t_gp %e %e\n",t_gp,t_gp/double(i));
    printf("t_lin %e %e\n",t_lin,t_lin/double(i));

}
