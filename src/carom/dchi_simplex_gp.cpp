#include "dchi_simplex_gp.h"

dchi_boundary_simplex_gp::dchi_boundary_simplex_gp(chisq_wrapper *cc,
                            gp_lin *gg, array_1d<int> &ii) :
                          dchi_boundary_simplex(cc, ii){

    _interpolator=gg;
    _last_real_call=0;

}

double dchi_boundary_simplex_gp::operator()(array_1d<double> &pt){
    double distance=associate_distance(pt);
    double mu;
    int i_found;

    mu=_interpolator[0](pt);
    if(_called>0 && _called-_last_real_call>50){
        mu=_chifn[0](pt);
        _last_real_call=_called;
    }

    double dmu;

    dmu=fabs(_chifn->target()-mu);

    double delta=_chifn->target()-_chifn->chimin();

    double exp_term;
    if(_chifn->target()<mu){
        exp_term=exp(-0.1*fabs(dmu)/delta);
    }
    else{
        exp_term=1.0;
    }

    _called++;
    return dmu-exp_term*delta*distance*2.0;
}
