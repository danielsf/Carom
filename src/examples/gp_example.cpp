#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "jellyBean.h"
#include "exampleLikelihoods.h"
#include "dalex_driver.h"

double distance_sq(const array_1d<double> &pt1,
                   const array_1d<double> &pt2,
                   const array_1d<double> &ell){


    double mu=0.0;
    int i;
    for(i=0;i<pt1.get_dim();i++){
        mu+=power((pt1.get_data(i)-pt2.get_data(i))/ell.get_data(i),2);
    }
    return sqrt(mu);
}


class GaussianProcess{

    public:

        GaussianProcess(array_2d<double> &pts, array_1d<double> &fn){
            _pts.set_name("gp_pts");
            _fn.set_name("gp_fn");
            _cov_inv.set_name("gp_cov_inv");
            _ell.set_name("gp_ell");
            _cq.set_name("gp_cq");
            _mean_bases.set_name("gp_mean_bases");
            _mean_coeffs.set_name("gp_mean_coeffs");
            _fn_mean.set_name("gp_fn_mean");
            int i,j;
            _pts.set_dim(pts.get_rows(), pts.get_cols());
            _fn.set_dim(pts.get_rows());
            _min_dex=-1;
            for(i=0;i<pts.get_rows();i++){
                _fn.set(i,fn.get_data(i));
                for(j=0;j<pts.get_cols();j++){
                    _pts.set(i,j,pts.get_data(i,j));
                }
                if(_min_dex<0 || fn.get_data(i)<fn.get_data(_min_dex)){
                    _min_dex=i;
                }
            }

            for(i=0;i<pts.get_rows();i++){
                _fn_mean.set(i,_mean(pts(i)));
            }

        }

        void build_cov_inv(const array_1d<double> &ell, double nugget){

            _cov_inv.set_dim(_pts.get_rows(),_pts.get_rows());
            _cov_vector.reset_preserving_room();

            int i,j;
            double ddsq;
            array_2d<double> cov;
            cov.set_dim(_pts.get_rows(),_pts.get_rows());

            for(i=0;i<_pts.get_cols();i++){
                _ell.set(i,raiseup(10.0,ell.get_data(i)));
            }

            _nugget=nugget;

            for(i=0;i<_pts.get_rows();i++){
                cov.set(i,i,1.0+_nugget);
                for(j=i+1;j<_pts.get_rows();j++){
                    ddsq = distance_sq(_pts(i), _pts(j), _ell);
                    cov.set(i,j,exp(-0.5*ddsq));
                    cov.set(j,i,exp(-0.5*ddsq));
                }
            }
            invert_lapack(cov, _cov_inv, 0);

            for(i=0;i<_pts.get_rows();i++){
                _cov_vector.set(i,0.0);
                for(j=0;j<_pts.get_rows();j++){
                    _cov_vector.add_val(i,_cov_inv.get_data(i,j)*(_fn.get_data(j)-_mean(_pts(j))));
                }
            }
        }

        double operator()(const array_1d<double> &pt){
            _cq.reset_preserving_room();
            _cq.set_dim(_pts.get_rows());
            int i;
            double ddsq;
            for(i=0;i<_pts.get_rows();i++){
                ddsq = distance_sq(pt,_pts(i),_ell);
                _cq.set(i,exp(-0.5*ddsq));
            }
            double ans=_mean(pt);
            int j;
            for(i=0;i<_pts.get_rows();i++){
                ans += _cq.get_data(i)*_cov_vector.get_data(i);
            }
            if(fabs(ans-95.3)<1.0e-4){
                for(i=0;i<_pts.get_rows();i++){
                    printf("%e %e\n",_cq.get_data(i),distance_sq(pt,_pts(i),_ell));
                }
                exit(1);
            }
            return ans;
        }

        double _mean(const array_1d<double> &pt){
            if(_mean_bases.get_rows()==0){
                _load_mean_model();
            }
            double ans=0.0;
            double mu;
            int i,j;
            for(i=0;i<pt.get_dim();i++){
                mu=0.0;
                for(j=0;j<pt.get_dim();j++){
                    mu+=(pt.get_data(j)-_pts.get_data(_min_dex,j))*_mean_bases.get_data(i,j);
                }
                ans += mu*mu*_mean_coeffs.get_data(i);
            }
            ans += _fn.get_data(_min_dex);
            return ans;
        }

    private:
        array_2d<double> _pts;
        array_1d<double> _fn;
        array_2d<double> _cov_inv;
        array_1d<double> _ell;
        array_1d<double> _cq;
        array_1d<double> _fn_mean;
        array_1d<double> _cov_vector;
        double _nugget;
        int _min_dex;
        array_2d<double> _mean_bases;
        array_1d<double> _mean_coeffs;

        void _load_mean_model(){
            FILE *in_file;
            in_file = fopen("output/test_180214/quad_bases.txt", "r");
            int i,j,k;
            double mu;
            _mean_bases.reset();
            _mean_bases.set_dim(_pts.get_cols(),_pts.get_cols());
            for(i=0;i<_pts.get_cols();i++){
                for(j=0;j<_pts.get_cols();j++){
                    k=fscanf(in_file,"%le",&mu);
                    if(k!=1){
                        printf("WARNING did not read in basis %d %d\n",i,j);
                        exit(1);
                    }
                    _mean_bases.set(i,j,mu);
                }
            }
            fclose(in_file);
            _mean_coeffs.reset();
            in_file=fopen("output/test_180214/quad_a_wgts.txt", "r");
            for(i=0;i<_pts.get_cols();i++){
                k=fscanf(in_file,"%le",&mu);
                if(k!=1){
                    printf("WARNING did not read in basis coeff %d\n",i);
                    exit(1);
                }
                _mean_coeffs.set(i,mu);
            }
        }



};

class gp_optimizer : public function_wrapper{

    public:
        gp_optimizer(){
            gp=NULL;
            _called=0;
            _pts.set_name("gp_opt_pts");
            _fn.set_name("gp_opt_fn");
        }

        ~gp_optimizer(){}

        virtual int get_called(){
            return _called;
        }

        void set_gp(GaussianProcess *gp_in){
            gp=gp_in;
        }

        void set_data(array_2d<double> &p, array_1d<double> &f){
            _pts.reset();
            _fn.reset();
            int i,j;
            for(i=0;i<p.get_rows();i++){
                _pts.add_row(p(i));
                _fn.add(f.get_data(i));
            }
            _min_rms=2.0*exception_value;
            _best_err=2.0*exception_value;
        }

        virtual double operator()(const array_1d<double> &ell){
            gp->build_cov_inv(ell, raiseup(10.0,ell.get_data(ell.get_dim()-1)));
            double err =0.0;
            double mu;
            int i;
            double delta;
            int mis_char=0;
            for(i=0;i<_pts.get_rows();i++){
                mu=gp[0](_pts(i));
                delta = power(mu-_fn.get_data(i),2);
                if(_fn.get_data(i)<116.0 && mu>116.0){
                    err+=4.0*delta;
                    mis_char++;
                }
                else if(_fn.get_data(i)>116.0 && mu<116.0){
                    err += 4.0*delta;
                    mis_char++;
                }
                else{
                    err+=delta/power(1.0+(_fn.get_data(i)-95.0)/5.0,2);
                }
            }
            double rms=sqrt(err/_pts.get_rows());
            if(err<_best_err){
                _best_err=err;
                _best_mis_char=mis_char;

            }
            printf("err %e best %e - %d\n",err,_best_err,_best_mis_char);
            return err;
        }

    private:
        GaussianProcess *gp;
        int _called;
        array_2d<double> _pts;
        array_1d<double> _fn;
        double _min_rms;
        int _best_mis_char;
        double _best_err;

};


int main(int iargc, char *argv[]){
    gaussianJellyBean12 chifn;
    int i,j;

    char data_name[letters];
    sprintf(data_name,"output/test_180214/output_midpt.txt_end_pts.txt");

    int dim=12;

    array_2d<double> pts;
    array_1d<double> row;
    array_1d<double> fn,min,max;

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
        fscanf(in_file,"%d",&j);
    }
    fclose(in_file);

    Ran chaos(99);

    array_2d<double> test_pts;
    array_1d<double> test_fn;
    char word[letters];
    double roll;
    int n_total_rows=0;
    int as_test;
    in_file = fopen("output/test_180214/output.txt", "r");
    for(i=0;i<dim+3;i++){fscanf(in_file,"%s",word);}
    printf("last word %s\n",word);
    while(fscanf(in_file,"%le",&mu)>0){
        n_total_rows++;
        row.set(0,mu);
        for(i=1;i<dim;i++){
            fscanf(in_file,"%le",&mu);
            row.set(i,mu);
        }
        fscanf(in_file,"%le",&mu);
        fscanf(in_file,"%d",&j);
        as_test = 0;
        if(mu<200.0){
            roll = chaos.doub();
            if(roll<0.01){
                as_test=1;
            }
        }
        if(as_test==1){
            test_pts.add_row(row);
            test_fn.add(mu);
        }
        /*else{
            if(mu<200.0){
                roll = chaos.doub();
                if(roll<0.001){
                    pts.add_row(row);
                    fn.add(mu);
                }
            }
        }*/
    }
    fclose(in_file);
    printf("total rows %d test rows %d\n",n_total_rows, test_pts.get_rows());
    printf("gp pts %d\n",pts.get_rows());


    GaussianProcess gp(pts, fn);

    gp_optimizer gp_opt;
    gp_opt.set_gp(&gp);
    gp_opt.set_data(test_pts, test_fn);

    array_2d<double> seed;
    seed.set_dim(dim+2,dim+1);
    for(i=0;i<dim+2;i++){
        for(j=0;j<dim;j++){
            seed.set(i,j,3.0*chaos.doub()+0.5);
        }
        seed.set(i,dim,-1.0*chaos.doub());
    }

    simplex_minimizer ffmin;
    ffmin.set_dice(&chaos);
    ffmin.use_gradient();
    ffmin.set_chisquared(&gp_opt);
    ffmin.set_abort_max_factor(100);
    array_1d<double> ell;

    ffmin.find_minimum(seed,ell);

    double nugget = 0.01;
    double mean;

    printf("building cov inv\n");
    gp.build_cov_inv(ell,raiseup(10.0,ell.get_data(ell.get_dim()-1)));
    printf("built covinv\n");
    FILE *out_file;
    out_file=fopen("junk.txt", "w");
    fprintf(out_file,"# truth fit delta\n");
    fprintf(out_file,"# ");
    for(i=0;i<ell.get_dim();i++){
        fprintf(out_file,"%e ",ell.get_data(i));
    }
    fprintf(out_file,"\n");
    for(i=0;i<test_pts.get_rows();i++){
        mu=gp(test_pts(i));
        mean = gp._mean(test_pts(i));
        fprintf(out_file,"%e %e %e %e\n",test_fn.get_data(i),mu,mean,fabs(test_fn.get_data(i)-mu));
    }
    fclose(out_file);

}
