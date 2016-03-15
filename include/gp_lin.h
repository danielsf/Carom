#ifndef GP_LIN_H
#define GP_LIN_H

#include "wrappers.h"
#include "eigen_wrapper.h"
#include "kd.h"

class gp_lin : public function_wrapper{

    public:

        ~gp_lin(){
            if(_built_here==1){
                if(_kd!=NULL){
                    delete _kd;
                }
                if(_fn!=NULL){
                    delete _fn;
                }
            }
        }

        gp_lin(){
            _kd=NULL;
            _fn=NULL;
            _built_here=1;
            _dexes.set_name("gp_dexes");
            _covarin.set_name("gp_covarin");
            _ell_factor=0.1;
            _nn=40;
            _nugget=1.0e-4;
            _mindex=0;
            _pt_when_mindex=0;
        }

        void build(array_2d<double>&, array_1d<double>&, array_1d<double>&, array_1d<double>&);

        void add_pt(array_1d<double>&, double);

        virtual double operator()(array_1d<double>&);

        void optimize(array_2d<double>&, array_1d<double>&);

        void set_ell_factor(double nn){
            _ell_factor=nn;
        }

        void set_kd_fn(kd_tree *kk, array_1d<double> *ff){
            if(_kd!=NULL && _built_here==1){
                printf("WARNING cannot set _kd; not null\n");
                exit(1);
            }

            if(_fn!=NULL && _built_here==1){
                printf("WARNING cannot set _fn; not null\n");
                exit(1);
            }

            _kd=kk;
            _fn=ff;
            _built_here=0;
        }


    private:

        kd_tree *_kd;
        array_1d<double> *_fn;
        int _built_here;

        int _nn;
        double _ell,_ell_factor,_fbar,_nugget;

        array_2d<double> _covarin;
        array_1d<int> _dexes;

        int _mindex,_pt_when_mindex;

};

#endif
