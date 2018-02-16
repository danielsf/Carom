#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "jellyBean.h"
#include "exampleLikelihoods.h"
#include "dalex_driver.h"

class quadratic_fitter : public function_wrapper{

    public:

        quadratic_fitter(){
            _pts.set_name("quadratic_fitter_pts");
            _pts_raw.set_name("quadratic_fitter_pts_raw");
            _fn.set_name("quadratic_fitter_fn");
            _fn_raw.set_name("quadrattic_fitter_fn_raw");
            _work_v.set_name("quadratic_fitter_work_v");
            _work_rot.set_name("quadratic_fitter_work_rot");
            _dotprod.set_name("quadratic_fitter_dot_prod");
            _sigma.set_name("quadratic_fitter_sigma");
            _random_bases.set_name("quadratic_fitter_random_bases");
            _called=0;
            chaos = new Ran(182311);
        }

        ~quadratic_fitter(){delete chaos;}

        virtual int get_called(){
            return _called;
        }

        void set_data(array_2d<double>&, array_1d<double>&,
                      array_1d<double>&);

        void set_prior_dir(array_2d<double>&, array_1d<double>&);

        virtual double operator()(const array_1d<double>&);

        double get_aa(){
            return _best_aa;
        }

        double get_v(int i){
            return _best_v.get_data(i);
        }

    private:
        array_2d<double> _pts,_pts_raw;
        array_1d<double> _fn,_fn_raw;
        array_1d<double> _work_v;
        array_1d<double> _work_rot;
        array_1d<double> _dotprod;
        array_1d<double> _sigma;
        array_1d<double> _best_v;
        array_2d<double> _random_bases;
        Ran *chaos;
        double _fn_min;
        int _min_dex;
        int _called;
        int _dim;
        double _best_err;
        double _best_aa;
};


void quadratic_fitter::set_data(array_2d<double> &pts,
                                array_1d<double> &fn,
                                array_1d<double> &sigma){

    _pts.reset();
    _pts_raw.reset();
    _fn.reset_preserving_room();
    _fn_raw.reset_preserving_room();
    _called=0;
    _dim = pts.get_cols();
    _dotprod.reset();
    _dotprod.set_dim(_pts.get_rows());
    _random_bases.reset();
    _random_bases.set_dim(_dim,_dim);
    _best_err=2.0*exception_value;

    int i;
    for(i=0;i<pts.get_rows();i++){
        _dotprod.set(i,0.0);
        _pts_raw.add_row(pts(i));
        _pts.add_row(pts(i));
        _fn.add(fn.get_data(i));
        _fn_raw.add(fn.get_data(i));
        _sigma.add(sigma.get_data(i));
        if(i==0 || fn.get_data(i)<_fn_min){
            _min_dex=i;
            _fn_min=fn.get_data(i);
        }
    }

    array_1d<double> min_pt;
    min_pt.set_name("quadratic_fitter_set_data_min_pt");
    for(i=0;i<_dim;i++){
        min_pt.set(i,_pts.get_data(_min_dex,i));
    }

    int j;

    for(i=0;i<_pts.get_rows();i++){
        for(j=0;j<_dim;j++){
            _pts_raw.subtract_val(i,j,min_pt.get_data(j));
            _pts.subtract_val(i,j,min_pt.get_data(j));
        }
        _fn.subtract_val(i,_fn_min);
        _fn_raw.subtract_val(i,_fn_min);
    }

    for(i=0;i<_dim;i++){
        for(j=0;j<_dim;j++){
            _random_bases.set(i,j,0.0);
        }
        _random_bases.set(i,i,1.0);
    }
}

void quadratic_fitter::set_prior_dir(array_2d<double> &prior_dir, array_1d<double> &aa){
    _best_err=2.0*exception_value;
    _work_rot.reset();
    _work_v.reset();
    _pts.reset();
    _fn.reset();
    _random_bases.reset();

    int i,j,k;
    for(i=0;i<_pts_raw.get_rows();i++){
        _pts.add_row(_pts_raw(i));
        _fn.add(_fn_raw.get_data(i));
    }

    double mu;
    for(i=0;i<_pts.get_rows();i++){
        for(j=0;j<prior_dir.get_rows();j++){
            mu=0.0;
            for(k=0;k<_dim;k++){
                mu+=prior_dir.get_data(j,k)*_pts.get_data(i,k);
            }
            for(k=0;k<_dim;k++){
                _pts.subtract_val(i,k,mu*prior_dir.get_data(j,k));
            }
            _fn.subtract_val(i,mu*mu*aa.get_data(j));
        }
    }

    array_1d<double> trial_basis;
    trial_basis.set_name("quadratic_fitter_trial_basis");

    _random_bases.set_cols(_dim);
    while(_random_bases.get_rows()<_dim-prior_dir.get_rows()){
        for(i=0;i<_dim;i++){
            trial_basis.set(i,normal_deviate(chaos,0.0,1.0));
        }
        for(i=0;i<prior_dir.get_rows();i++){
            mu=0.0;
            for(j=0;j<_dim;j++){
                mu+=trial_basis.get_data(j)*prior_dir.get_data(i,j);
            }
            for(j=0;j<_dim;j++){
                trial_basis.subtract_val(j,mu*prior_dir.get_data(i,j));
            }
        }
        for(i=0;i<_random_bases.get_rows();i++){
            mu=0.0;
            for(j=0;j<_dim;j++){
                mu+=trial_basis.get_data(j)*_random_bases.get_data(i,j);
            }
            for(j=0;j<_dim;j++){
                trial_basis.subtract_val(j,mu*_random_bases.get_data(i,j));
            }
        }
        mu=trial_basis.normalize();
        if(mu>1.0e-20){
            _random_bases.add_row(trial_basis);
        }
    }

}

double quadratic_fitter::operator()(const array_1d<double> &theta){
    _called++;
    if(_work_v.get_dim()!=_dim){
        _work_v.reset();
        _work_v.set_dim(_dim);
    }
    if(_work_rot.get_dim()!=theta.get_dim()+1){
       _work_rot.reset();
       _work_rot.set_dim(theta.get_dim()+1);
    }

    int i,j;
    for(i=0;i<theta.get_dim()+1;i++){
        _work_rot.set(i,1.0);
    }

    double sin_theta;
    double cos_theta;

    for(i=0;i<theta.get_dim();i++){
        sin_theta=sin(theta.get_data(i));
        cos_theta=cos(theta.get_data(i));
        _work_rot.multiply_val(i,sin_theta);
        for(j=i+1;j<theta.get_dim()+1;j++){
            _work_rot.multiply_val(j,cos_theta);
        }
    }

    for(i=0;i<_dim;i++){
        _work_v.set(i,0.0);
        for(j=0;j<_random_bases.get_rows();j++){
            _work_v.add_val(i,_work_rot.get_data(j)*_random_bases.get_data(j,i));
        }
    }

    double mu=0.0;

    for(i=0;i<_dim;i++){
        mu+=power(_work_v.get_data(i),2);
    }
    if(fabs(mu-1.0)>1.0e-5){
       printf("WARNING _work_v**2 %e\n",mu);
       exit(1);
    }

    for(i=0;i<_pts.get_rows();i++){
        mu=0.0;
        for(j=0;j<_dim;j++){
            mu+=_pts.get_data(i,j)*_work_v.get_data(j);
        }
        _dotprod.set(i,mu*mu);
    }

    double aa=0.0;
    for(i=0;i<_dotprod.get_dim();i++){
        aa+=_fn.get_data(i)*_dotprod.get_data(i)/power(_sigma.get_data(i),2);
    }
    double denom=0.0;
    for(i=0;i<_dotprod.get_dim();i++){
        denom+=_dotprod.get_data(i)*_dotprod.get_data(i)/power(_sigma.get_data(i),2);
    }
    aa/=denom;
    double err=0.0;
    for(i=0;i<_dotprod.get_dim();i++){
        err+=power((_fn.get_data(i)-aa*_dotprod.get_data(i))/_sigma.get_data(i),2);
    }

    if(err<_best_err){
        _best_err=err;
        _best_aa=aa;
        for(i=0;i<_dim;i++){
            _best_v.set(i,_work_v.get_data(i));
        }
    }

    return err;

}



int main(int iargc, char *argv[]){
    gaussianJellyBean12 chifn;
    int i,j;

    array_2d<double> out_dir;
    array_1d<double> aa;

    char data_name[letters];
    sprintf(data_name,"output/test_180214/output_midpt.txt_end_pts.txt");

    int dim=12;
    out_dir.set_cols(dim);

    array_2d<double> pts;
    array_1d<double> row;
    array_1d<double> fn,min,max;
    array_1d<double> sigma;

    double known_min=95.3;
    double deltachi=21.03;

    for(i=0;i<dim;i++){
        min.set(i,2.0*exception_value);
        max.set(i,-2.0*exception_value);
    }

    double mu;

    FILE *in_file;
    in_file=fopen(data_name, "r");
    while(fscanf(in_file,"%le",&mu)>0){
        row.set(0,mu);
        if(mu<min.get_data(0)){
            min.set(0,mu);
        }
        if(mu>max.get_data(0)){
            max.set(0,mu);
        }
        for(i=1;i<dim;i++){
            fscanf(in_file,"%le",&mu);
            row.set(i,mu);
            if(mu<min.get_data(i)){
                min.set(i,mu);
            }
            if(mu>max.get_data(i)){
                max.set(i,mu);
            }
        }
        pts.add_row(row);
        fscanf(in_file,"%le",&mu);
        fn.add(mu);
        if(mu>known_min+deltachi){
            sigma.add(1.0+(mu-known_min)/5.0);
        }
        else{
            sigma.add(1.0);
        }
        fscanf(in_file,"%d",&j);
    }
    fclose(in_file);

    int min_dex;
    double fn_min;

    for(i=0;i<fn.get_dim();i++){
        if(i==0 || fn.get_data(i)<fn_min){
            fn_min=fn.get_data(i);
            min_dex=i;
        }
    }

    /*array_1d<double> dd,dd_sorted;
    array_1d<int> dd_dex;
    for(i=0;i<pts.get_rows();i++){
        mu=0.0;
        for(j=0;j<dim;j++){
            mu+=power(pts.get_data(i,j)-pts.get_data(min_dex,j),2);
        }
        dd.add(mu);
        dd_dex.add(i);
    }
    sort(dd,dd_sorted,dd_dex);
    for(i=0;i<pts.get_rows();i++){
        if(dd.get_data(i)<dd_sorted.get_data(10)){
            sigma.set(i,1.0);
        }
    }*/
    printf("penultimate %e\n",sigma.get_data(sigma.get_dim()-2));


    quadratic_fitter q_fit;
    q_fit.set_data(pts,fn,sigma);

    array_1d<double> th_min,th_max;
    array_2d<double> seed;



    Ran chaos(99);

    int i_dir;

    int active_dim;

    array_1d<double> minpt;

    simplex_minimizer *ffmin;
    ffmin = NULL;

    active_dim=dim;
    for(i_dir=0;i_dir<dim-1;i_dir++){
        active_dim--;
        printf("active dim %d\n",active_dim);

        if(ffmin != NULL){
            delete ffmin;
        }

        if(out_dir.get_rows()>0){
            q_fit.set_prior_dir(out_dir,aa);
        }

        th_min.reset();
        th_max.reset();
        for(i=0;i<active_dim;i++){
            th_min.set(i,0.0);
            th_max.set(i,2.0*pi);
        }

        ffmin = new simplex_minimizer;
        ffmin->set_chisquared(&q_fit);
        ffmin->set_minmax(th_min,th_max);
        ffmin->set_dice(&chaos);
        ffmin->use_gradient();
        ffmin->set_abort_max_factor(10000);

        minpt.reset_preserving_room();
        seed.reset();
        seed.set_dim(active_dim+1,active_dim);

        for(i=0;i<seed.get_rows();i++){
            for(j=0;j<seed.get_cols();j++){
                seed.set(i,j,pi+2.0*(chaos.doub()-0.5)*pi);
            }
        }

        ffmin->find_minimum(seed,minpt);

        for(i=0;i<dim;i++){
            out_dir.set(i_dir,i,q_fit.get_v(i));
        }

        aa.add(q_fit.get_aa());
    }

    if(out_dir.get_rows()!=dim-1){
        printf("WARNING only got %d dir\n",out_dir.get_rows());
        exit(1);
    }

    array_1d<double> trial;
    for(i=0;i<dim;i++){
        trial.set(i,normal_deviate(&chaos,0.0,1.0));
    }
    for(i=0;i<out_dir.get_rows();i++){
        mu=0.0;
        for(j=0;j<dim;j++){
            mu+=out_dir.get_data(i,j)*trial.get_data(j);
        }
        for(j=0;j<dim;j++){
            trial.subtract_val(j,mu*out_dir.get_data(i,j));
        }
    }
    trial.normalize();
    out_dir.add_row(trial);
    int k;
    for(i=0;i<dim;i++){
        mu=0.0;
        for(j=0;j<dim;j++){
            mu+=out_dir.get_data(i,j)*out_dir.get_data(i,j);
        }
        if(fabs(mu-1.0)>1.0e-5){
            printf("WARNING normsq of %d %e\n",i,mu);
            exit(1);
        }
        for(j=i+1;j<dim;j++){
            mu=0.0;
            for(k=0;k<dim;k++){
                mu+=out_dir.get_data(i,k)*out_dir.get_data(j,k);
            }
            if(fabs(mu)>1.0e-5){
                printf("WARNING %d dot %d %e\n",i,j,mu);
                exit(1);
            }
        }
    }

    array_1d<double> dot_product;


    for(i=0;i<pts.get_rows();i++){
        dot_product.set(i,0.0);
        for(j=0;j<dim;j++){
            dot_product.add_val(i,(pts.get_data(i,j)-pts.get_data(min_dex,j))*out_dir.get_data(0,j));
        }
    }


    FILE *out_file;
    out_file=fopen("quad_test_output.txt", "w");
    for(i=0;i<pts.get_rows();i++){
        mu=fn.get_data(min_dex)+aa.get_data(0)*dot_product.get_data(i)*dot_product.get_data(i);
        fprintf(out_file,"%e %e %e\n",fn.get_data(i),mu,dot_product.get_data(i));
    }
    fclose(out_file);

    printf("%d\n",sigma.get_dim());

    out_file = fopen("quad_bases.txt", "w");
    for(i=0;i<out_dir.get_rows();i++){
        for(j=0;j<out_dir.get_cols();j++){
            fprintf(out_file,"%e ",out_dir.get_data(i,j));
        }
        fprintf(out_file,"\n");
    }
    fclose(out_file);

    out_file = fopen("quad_wgts.txt", "w");
    for(i=0;i<sigma.get_dim();i++){
        fprintf(out_file,"%e\n",sigma.get_data(i));
    }
    fclose(out_file);

}
