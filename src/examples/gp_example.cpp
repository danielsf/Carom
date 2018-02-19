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

        ~GaussianProcess(){}

        GaussianProcess(){
            _pts.set_name("gp_pts");
            _fn.set_name("gp_fn");
            _cov_inv.set_name("gp_cov_inv");
            _ell.set_name("gp_ell");
            _cq.set_name("gp_cq");
            _mean_bases.set_name("gp_mean_bases");
            _mean_coeffs.set_name("gp_mean_coeffs");
            _fn_mean.set_name("gp_fn_mean");
            _proj_1.set_name("gp_proj_1");
            _proj_arr.set_name("gp_proj_2");
       }

        void build(array_2d<double> &pts, array_1d<double> &fn){
            int i,j;
            _pts.reset();
            _fn.reset();
            _proj_arr.reset();
            _cov_inv.reset();
            _ell.reset();
            _cq.reset();
            _fn_mean.reset();
            _cov_vector.reset();
            _mean_bases.reset();
            _mean_coeffs.reset();
            _proj_1.reset();
            _proj_arr.reset();

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

        void build_cov_inv(const array_1d<double> &ell){

            if(ell.get_dim()!=_pts.get_cols()+2){
                printf("WARNING ell %d pts_dim %d\n",ell.get_dim(),_pts.get_cols());
                exit(1);
            }

            _cov_inv.set_dim(_pts.get_rows(),_pts.get_rows());
            _cov_vector.reset_preserving_room();
            _proj_arr.reset();
            _failed_cut=0;
            _passed_cut=0;

            double nugget=raiseup(10.0,ell.get_data(ell.get_dim()-2));
            array_2d<double> dist_arr;
            dist_arr.set_name("dist_arr");
            array_1d<double> dist_to_sort,dist_sorted;
            array_1d<int> dist_dex;
            dist_to_sort.set_name("dist_to_sort");
            dist_sorted.set_name("dist_sorted");
            dist_dex.set_name("dist_dex");
            dist_arr.set_dim(_pts.get_rows(),_pts.get_rows());

            int i,j;
            double ddsq,dd;
            array_2d<double> cov;
            cov.set_dim(_pts.get_rows(),_pts.get_rows());

            for(i=0;i<_pts.get_cols();i++){
                _ell.set(i,raiseup(10.0,ell.get_data(i)));
            }

            _nugget=nugget;

            for(i=0;i<_pts.get_rows();i++){
                cov.set(i,i,1.0+_nugget);
                for(j=i+1;j<_pts.get_rows();j++){
                    ddsq = cov_distance_sq(i,j);
                    dd=sqrt(ddsq);
                    dist_arr.set(i,j,dd);
                    dist_arr.set(j,i,dd);
                    dist_to_sort.add(dd);
                    dist_dex.add(dist_to_sort.get_dim()-1);
                    cov.set(i,j,exp(-0.5*dd));
                    cov.set(j,i,exp(-0.5*dd));
                }
            }

            sort(dist_to_sort, dist_sorted, dist_dex);
            double frac = raiseup(10.0,ell.get_data(ell.get_dim()-1));
            int cutoff_dex;
            if(frac<0.0){
                _cutoff=0.0;
            }
            else if(frac>1.0){
                _cutoff=dist_sorted.get_data(dist_sorted.get_dim()-1);
            }
            else{
                cutoff_dex = int(dist_sorted.get_dim()*frac);
                if(cutoff_dex>=dist_sorted.get_dim()){
                    cutoff_dex=dist_sorted.get_dim()-1;
                }
                _cutoff=dist_sorted.get_data(cutoff_dex);
            }
            for(i=0;i<_pts.get_rows();i++){
                for(j=i+1;j<_pts.get_rows();j++){
                    if(dist_arr.get_data(i,j)>_cutoff){
                        cov.set(i,j,0.0);
                        cov.set(j,i,0.0);
                        _failed_cut++;
                    }
                    else{
                        _passed_cut++;
                    }
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

        void _build_proj(){
           int i,j,k;
           _proj_arr.set_dim(_pts.get_rows(),_pts.get_cols());
           for(i=0;i<_pts.get_rows();i++){
               for(j=0;j<_pts.get_cols();j++){
                   _proj_arr.set(i,j,0.0);
                   for(k=0;k<_pts.get_cols();k++){
                       _proj_arr.add_val(i,j,(_pts.get_data(i,k)-_pts.get_data(_min_dex,k))*
                                              _mean_bases.get_data(j,k));
                   }
               }
           }
        }

        double cov_distance_sq(int d1, int d2){
           if(_proj_arr.get_rows()==0){
               _build_proj();
           }
           int i;
           double ans=0.0;
           for(i=0;i<_pts.get_cols();i++){
               ans+=power((_proj_arr.get_data(d1,i)-_proj_arr.get_data(d2,i))/_ell.get_data(i),2);
           }
           return ans;
        }

        double cov_distance_sq(const array_1d<double> &pt, int dex){
            int i,j;
            if(_proj_arr.get_rows()==0){
               _build_proj();
            }
            /*for(i=0;i<_pts.get_cols();i++){
                _proj_1.set(i,0.0);
                for(j=0;j<_pts.get_cols();j++){
                    _proj_1.add_val(i,(pt.get_data(j)-_pts.get_data(_min_dex,j))*_mean_bases.get_data(i,j));
                }
            }*/
            double ans=0.0;
            for(i=0;i<_pts.get_cols();i++){
                ans+=power((_proj_1.get_data(i)-_proj_arr.get_data(dex,i))/_ell.get_data(i),2);
            }
            return ans;
        }

        double operator()(const array_1d<double> &pt){
            _cq.reset_preserving_room();
            _cq.set_dim(_pts.get_rows());
            int i;
            double ddsq;
            double dd;
            double ans=_mean(pt);
            for(i=0;i<_pts.get_rows();i++){
                ddsq = cov_distance_sq(pt,i);
                dd=sqrt(ddsq);
                if(dd<_cutoff){
                    _cq.set(i,exp(-0.5*dd));
                }
                else{
                    _cq.set(i,0.0);
                }

            }
            int j;
            for(i=0;i<_pts.get_rows();i++){
                ans += _cq.get_data(i)*_cov_vector.get_data(i);
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
                _proj_1.set(i,mu);
                ans += mu*mu*_mean_coeffs.get_data(i);
            }
            ans += _fn.get_data(_min_dex);
            return ans;
        }


        int _failed_cut,_passed_cut;

    private:
        array_2d<double> _pts;
        array_1d<double> _fn;
        array_2d<double> _cov_inv;
        array_1d<double> _ell;
        array_1d<double> _cq;
        array_1d<double> _fn_mean;
        array_1d<double> _cov_vector;
        double _nugget,_cutoff;
        int _min_dex;
        array_2d<double> _mean_bases;
        array_1d<double> _mean_coeffs;
        array_1d<double> _proj_1;
        array_2d<double> _proj_arr;

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

        void reset_called(){
            _called=0;
            _best_err=2.0*exception_value;
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
            _called++;
            gp->build_cov_inv(ell);
            double err =0.0;
            double err_mean=0.0;
            double mu;
            double mu_mean;
            int i;
            double delta,delta_mean;
            int mis_char=0;
            double wgt;
            for(i=0;i<_pts.get_rows();i++){
                mu=gp[0](_pts(i));
                mu_mean=gp[0]._mean(_pts(i));
                delta = power(mu-_fn.get_data(i),2);
                delta_mean=power(mu_mean-_fn.get_data(i),2);
                if(_fn.get_data(i)<116.0 && mu>116.0){
                    mis_char++;
                    wgt=1.0;
                }
                else if(_fn.get_data(i)>116.0 && mu<116.0){
                    mis_char++;
                    wgt=1.0;
                }
                else{
                    wgt=1.0/power(1.0+(_fn.get_data(i)-95.0)/5.0,2);
                }

                err+=delta*wgt;

                err_mean+=delta_mean*wgt;

            }
            double rms=sqrt(err/_pts.get_rows());
            if(err<_best_err){
                _best_err=err;
                _best_mis_char=mis_char;
                printf("err %.3e err_mean %.3e best %e - %d log_cutoff %.2e -- %d %d\n",
                err,err_mean,_best_err,_best_mis_char,ell.get_data(ell.get_dim()-1),
                gp->_failed_cut,gp->_passed_cut);
            }
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

    array_2d<double> gp_pts;
    array_1d<double> row;
    array_1d<double> gp_fn,min,max;

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
        gp_pts.add_row(row);
        fscanf(in_file,"%le",&mu);
        gp_fn.add(mu);
        fscanf(in_file,"%d",&j);
    }
    fclose(in_file);

    Ran chaos(99);

    array_2d<double> test_pts;
    array_1d<double> test_fn;

    array_2d<double> final_pts;
    array_1d<double> final_fn;

    array_2d<double> pool_pts;
    array_1d<double> pool_fn;

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
        else{
            if(mu<200.0){
                roll = chaos.doub();
                if(roll<0.05){
                    final_pts.add_row(row);
                    final_fn.add(mu);
                }
                else{
                    roll = chaos.doub();
                    if(roll<0.01){
                        pool_pts.add_row(row);
                        pool_fn.add(mu);
                    }
                }
            }
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
    printf("test %d pool %d final %d\n",
    test_pts.get_rows(),pool_pts.get_rows(),final_pts.get_rows());

    printf("gp pts %d\n",gp_pts.get_rows());

    array_1d<double> ell;
    array_2d<double> seed;


    double mean;

    double err=0.0;
    double err_mean=0.0;

    FILE *out_file;
    simplex_minimizer *ffmin;
    ffmin = NULL;


    GaussianProcess gp;

    gp_optimizer gp_opt;
    gp_opt.set_gp(&gp);
    gp_opt.set_data(test_pts, test_fn);

    int iteration;
    int use_it=0;
    int added;
    for(iteration=0;iteration<2;iteration++){

        gp_opt.reset_called();
        gp.build(gp_pts, gp_fn);

        seed.set_dim(dim+3,dim+2);
        for(i=0;i<dim+2;i++){
            for(j=0;j<dim;j++){
                seed.set(i,j,3.0*(chaos.doub()-0.5));
            }
            seed.set(i,dim,-1.0*chaos.doub());
            seed.set(i,dim+1,-1.0*chaos.doub());
        }

        if(ffmin != NULL){
            delete ffmin;
        }

        ffmin = new simplex_minimizer;
        ffmin->set_dice(&chaos);
        ffmin->use_gradient();
        ffmin->set_chisquared(&gp_opt);
        ffmin->set_abort_max_factor(10);

        ffmin->find_minimum(seed,ell);

        gp.build_cov_inv(ell);
        added=0;
        for(i=0;i<pool_pts.get_rows();i++){
            mu=gp(pool_pts(i));
            use_it=0;
            if(pool_fn.get_data(i)>116.0 && mu<100.0){
                use_it=1;
            }
            if(pool_fn.get_data(i)<110.0 && fabs(mu-pool_fn.get_data(i))>1.0){
                use_it=1;
            }
            if(use_it==1){
                added++;
                gp_pts.add_row(pool_pts(i));
                gp_fn.add(pool_fn.get_data(i));
            }
        }
        printf("added %d == cutoff %e\n",added,ell.get_data(ell.get_dim()-1));

    }

    err=0.0;
    err_mean=0.0;

    printf("building cov inv\n");
    gp.build_cov_inv(ell);
    printf("built covinv\n");
    out_file=fopen("junk.txt", "w");
    fprintf(out_file,"# truth fit delta\n");
    fprintf(out_file,"# ");
    for(i=0;i<ell.get_dim();i++){
        fprintf(out_file,"%e ",ell.get_data(i));
    }
    fprintf(out_file,"\n");
    for(i=0;i<final_pts.get_rows();i++){
        mu=gp(final_pts(i));
        mean = gp._mean(final_pts(i));
        err+=power(mu-final_fn.get_data(i),2);
        err_mean+=power(mean-final_fn.get_data(i),2);
        fprintf(out_file,"%e %e %e %e -- ",final_fn.get_data(i),mu,mean,fabs(final_fn.get_data(i)-mu));
        for(j=0;j<dim;j++){
            fprintf(out_file,"%e ",final_pts.get_data(i,j));
        }
        fprintf(out_file,"\n");
    }
    fclose(out_file);
    printf("err %e err_mean %e\n",err,err_mean);

}
