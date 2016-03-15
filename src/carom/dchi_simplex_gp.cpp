#include "dchi_simplex_gp.h"


double dchi_boundary_simplex_gp::operator()(array_1d<double> &pt){
    _called++;
    double mu=_chifn[0](pt);
    if(mu<_chifn->target()){
        return mu;
    }

    double predicted=_interpolator[0](pt);
    return mu-predicted;

}
