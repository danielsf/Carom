#ifndef GP_H
#define GP_H

#include "wrappers.h"
#include "eigen_wrapper.h"
#include "kd.h"

class gp : public function_wrapper{

    public:
        ~gp(){}
        gp(){
            _fn.set_name("gp_fn");
            _dexes.set_name("gp_dexes");
            _covarin.set_name("gp_covarin");
            _ell_factor=0.1;
            _nn=40;
            _nugget=1.0e-4;
        }
        void build(array_2d<double>&, array_1d<double>&, array_1d<double>&, array_1d<double>&);
        void add_pt(array_1d<double>&, double);
        virtual double operator()(array_1d<double>&);
        void optimize(array_2d<double>&, array_1d<double>&);
        void set_ell_factor(double nn){
            _ell_factor=nn;
        }


    private:

        kd_tree _kd;

        int _nn;
        double _ell,_ell_factor,_fbar,_nugget;
        array_1d<double> _fn;

        array_2d<double> _covarin;
        array_1d<int> _dexes;

        void _set_covarin();
        double _covariogram(array_1d<double>&, array_1d<double>&);


};

#endif
