#ifndef DCHI_SIMPLEX_GP_H
#define DCHI_SIMPLEX_GP_H

#include "gp_lin.h"
#include "dchi_simplex.h"

class dchi_boundary_simplex_gp : public dchi_boundary_simplex{

    public:
        dchi_boundary_simplex_gp(chisq_wrapper*, gp_lin*, array_1d<int>&);
        ~dchi_boundary_simplex_gp(){};

        virtual double operator()(array_1d<double>&);

    private:
        gp_lin *_interpolator;

};

#endif
