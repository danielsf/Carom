#ifndef DCHI_SIMPLEX_GP_H
#define DCHI_SIMPLEX_GP_H

#include "gp_lin.h"
#include "dchi_simplex.h"

class dchi_boundary_simplex_gp : public function_wrapper{

    public:
        dchi_boundary_simplex_gp(chisq_wrapper *cc, gp_lin *gg){
            _chifn=cc;
            _interpolator=gg;
            _called=0;
        }
        ~dchi_boundary_simplex_gp(){};

        virtual double operator()(array_1d<double>&);
        virtual int get_called(){
            return _called;
        }

    private:
        chisq_wrapper *_chifn;
        gp_lin *_interpolator;
        int _called;

};

#endif
