#include "gp.h"

void gp::build(array_2d<double> &pt_in, array_1d<double> &fn_in,
               array_1d<double> &min_in, array_1d<double> &max_in){

    int i;
    for(i=0;i<fn_in.get_dim();i++){
        _fn.set(i,fn_in.get_data(i));
    }

    _kd.build_tree(pt_in,min_in,max_in);

}

void gp::add_pt(array_1d<double> &pt, double ff){
    _kd.add(pt);
    _fn.add(ff);
}

double gp::operator()(array_1d<double> &pt){

    array_1d<int> local_dex;
    local_dex.set_name("gp_operator_local_dex");
    array_1d<double> local_distance;
    local_distance.set_name("gp_operator_local_distance");

    _kd.nn_srch(pt, _nn, local_dex, local_distance);

    int i,j;

    array_1d<double> distance,distance_sorted;
    array_1d<int> distance_dex;
    int ct=0;
    for(i=0;i<_nn;i++){
        for(j=i+1;j<_nn;j++){
            distance.add(_kd.distance(local_dex.get_data(i), local_dex.get_data(j)));
            distance_dex.add(ct);
            ct++;
        }
    }

    sort_and_check(distance,distance_sorted,distance_dex);
    _ell=_ell_factor*distance_sorted.get_data(ct/2);

    //printf("\n    nearest dd %e -- %e\n",
    //local_distance.get_data(0),_fn.get_data(local_dex.get_data(0)));
    //printf("    mid fn %e\n",_fn.get_data(local_dex.get_data(_nn/2)));

    int use_it=0;
    if(_dexes.get_dim()!=local_dex.get_dim()){
        use_it=1;
    }
    else{
        for(i=0;i<local_dex.get_dim();i++){
            if(_dexes.contains(local_dex.get_data(i))==0){
                use_it=1;
            }
        }
    }

    if(use_it==1){
        _dexes.reset_preserving_room();
        for(i=0;i<local_dex.get_dim();i++){
            _dexes.set(i,local_dex.get_data(i));
        }
        _set_covarin();
    }

    array_1d<double> gq;
    gq.set_name("gp_operator_gq");
    for(i=0;i<_dexes.get_dim();i++){
        gq.set(i, _covariogram(pt, _kd.get_pt(_dexes.get_data(i))[0]));
    }

    double ans=_fbar;
    for(i=0;i<_dexes.get_dim();i++){
        for(j=0;j<_dexes.get_dim();j++){
            ans+=_covarin.get_data(i,j)*gq.get_data(i)*(_fn.get_data(_dexes.get_data(j))-_fbar);
        }
    }

    return ans;

}

void gp::_set_covarin(){
    _fbar=0.0;
    int i,j,ct;

    for(i=0;i<_dexes.get_dim();i++){
        _fbar+=_fn.get_data(_dexes.get_data(i));
    }
    _fbar=_fbar/double(_dexes.get_dim());
 
    array_2d<double> covar;
    covar.set_cols(_dexes.get_dim());
    covar.set_name("gp_set_covarin_covar");
    double minterm,maxterm;
    minterm=2.0*exception_value;
    maxterm=-2.0*exception_value;
    for(i=0;i<_dexes.get_dim();i++){
        for(j=i;j<_dexes.get_dim();j++){
            covar.set(i,j,_covariogram(_kd.get_pt(_dexes.get_data(i))[0],_kd.get_pt(_dexes.get_data(j))[0]));
            if(i==j){
                covar.add_val(i,j,_nugget);
            }
        }
    }

    array_1d<double> terms,sorted_terms;
    array_1d<int> term_dexes;

    invert_lapack(covar, _covarin, 0);

    ct=0;
    for(i=0;i<_covarin.get_rows();i++){
        for(j=0;j<_covarin.get_cols();j++){
            terms.add(fabs(_covarin.get_data(i,j)));
            term_dexes.add(ct);
            ct++;
            if(fabs(_covarin.get_data(i,j))<minterm){
                minterm=fabs(_covarin.get_data(i,j));
            }
            if(fabs(_covarin.get_data(i,j))>maxterm){
                maxterm=fabs(_covarin.get_data(i,j));
            }
        }
    }

    sort_and_check(terms,sorted_terms,term_dexes);

    //printf("    fb %e\n",_fbar);
    //printf("    minterm %e\n",minterm);
    //printf("    maxterm %e\n",maxterm);
    //printf("    medterm %e\n",sorted_terms.get_data(term_dexes.get_dim()/2));

}


double gp::_covariogram(array_1d<double> &p1, array_1d<double> &p2){
    double dd=_kd.distance(p1,p2);
    return exp(-0.5*dd/_ell);
}


void gp::optimize(array_2d<double> &pts, array_1d<double> &ff){

    double log_ell,ell_best;
    double cost,cost_best,ln10;
    int i;
    array_1d<double> cost_arr,cost_sorted;
    array_1d<int> cost_dex;
    cost_best=2.0*exception_value;
    ell_best=-1.0;
    ln10=log(10.0);
    for(log_ell=-1.0;log_ell<3.0;log_ell+=0.5){
        _ell_factor=exp(ln10*log_ell);
        for(i=0;i<pts.get_rows();i++){
            cost=fabs(this[0](pts(i)[0])-ff.get_data(i))/ff.get_data(i);
            cost_arr.set(i,cost);
            cost_dex.set(i,i);
        }
        sort_and_check(cost_arr,cost_sorted,cost_dex);
        cost=cost_sorted.get_data(pts.get_rows()/2);
        if(cost<cost_best || ell_best<0.0){
             //printf("   best cost %e best ell %e\n",cost,_ell);
             cost_best=cost;
             ell_best=_ell_factor;
        }
        printf("    cost %.2e min %.2e quart %.2e ell %.2e\n",
        cost,cost_sorted.get_data(0),
        cost_sorted.get_data(pts.get_rows()/4),_ell_factor);
    }

    printf("    best cost %e ell %e\n",cost_best,ell_best);
    _ell_factor=ell_best;

}
