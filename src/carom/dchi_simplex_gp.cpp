#include "dchi_simplex_gp.h"

dchi_boundary_simplex_gp::dchi_boundary_simplex_gp(chisq_wrapper *cc,
                            gp_lin *gg, array_1d<int> &ii) :
                          dchi_boundary_simplex(cc, ii){

    _interpolator=gg;

}

double dchi_boundary_simplex_gp::operator()(array_1d<double> &pt){
    double distance=associate_distance(pt);
    double mu;
    int i_found;
    mu=_interpolator[0](pt);

    double dmu=fabs(_chisq->target()-mu);
    double delta=_chisq->target()-_chisq->chimin();

    double exp_term;
    if(_chisq->target()<mu){
        exp_term=exp(-0.1*dmu/delta);
    }
    else{
        exp_term=1.0;
    }

    _called++;
    return dmu-exp_term*delta*distance*2.0;
}
