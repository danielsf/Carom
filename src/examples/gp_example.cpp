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
            int i,j;
            _pts.set_dim(pts.get_rows(), pts.get_cols());
            _fn.set_dim(pts.get_rows());
            _fn_mean = 0.0;
            for(i=0;i<pts.get_rows();i++){
                _fn.set(i,fn.get_data(i));
                _fn_mean+=fn.get_data(i);
                for(j=0;j<pts.get_cols();j++){
                    _pts.set(i,j,pts.get_data(i,j));
                }
            }
            _fn_mean = _fn_mean/double(pts.get_rows());
            _fn_mean=95.3;
        }
    
        void build_cov_inv(array_1d<double> &ell, double nugget){

            _cov_inv.set_dim(_pts.get_rows(),_pts.get_rows());

            int i,j;
            double ddsq;
            array_2d<double> cov;
            cov.set_dim(_pts.get_rows(),_pts.get_rows());

            for(i=0;i<_pts.get_cols();i++){
                _ell.set(i,ell.get_data(i));
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
            double ans=_fn_mean;
            int j;
            for(i=0;i<_pts.get_rows();i++){
                for(j=0;j<_pts.get_rows();j++){
                    ans += _cq.get_data(i)*_cov_inv.get_data(i,j)*(_fn.get_data(j)-_fn_mean);
                }
            }
            if(fabs(ans-95.3)<1.0e-4){
                for(i=0;i<_pts.get_rows();i++){
                    printf("%e %e\n",_cq.get_data(i),distance_sq(pt,_pts(i),_ell));
                }
                exit(1);
            }
            return ans;
        }

    private:
        array_2d<double> _pts;
        array_1d<double> _fn;
        array_2d<double> _cov_inv;
        array_1d<double> _ell;
        array_1d<double> _cq;
        double _nugget;
        double _fn_mean;

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


    double ell_factor = 0.1;
    double nugget = 0.01;
    array_1d<double> ell;
    for(i=0;i<dim;i++){
        ell.set(i,2.0);
    }
    printf("building cov inv\n");
    gp.build_cov_inv(ell,nugget);
    printf("built covinv\n");
    FILE *out_file;
    out_file=fopen("junk.txt", "w");
    fprintf(out_file,"# truth fit delta\n");
    for(i=0;i<test_pts.get_rows();i++){
        mu=gp(test_pts(i));
        fprintf(out_file,"%e %e %e\n",test_fn.get_data(i),mu,fabs(test_fn.get_data(i)-mu));
    }
    fclose(out_file);

}
