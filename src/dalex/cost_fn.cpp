#include "cost_fn.h"

cost_fn::cost_fn(chisq_wrapper *cc, array_1d<int> &aa){
    _envelope=1.0;
    build(cc,aa,1);
}

cost_fn::cost_fn(chisq_wrapper *cc, array_1d<int> &aa, int min_or_med){
    _envelope=1.0;
    build(cc,aa,min_or_med);
}

void cost_fn::build(chisq_wrapper *cc, array_1d<int> &aa, int min_or_med){

    printf("building cost_fn with %d associates\n",aa.get_dim());
    _called=0;

    _pt_cache.set_name("dchi_interior_pt_cache");
    _fn_cache.set_name("dchi_interior_fn_cache");
    _associates.set_name("dchi_interior_fn_associates");
    _relative_norm.set_name("dchi_interior_relative_norm");

    _associates.reset_preserving_room();
    _fn_cache.reset_preserving_room();
    _pt_cache.reset_preserving_room();

    _chifn=cc;

    int i;

    for(i=0;i<aa.get_dim();i++){
        _associates.add(aa.get_data(i));
    }

    array_1d<double> norm;
    array_1d<double> min,max;
    norm.set_name("cost_fn_build_norm");
    min.set_name("cost_fn_build_min");
    max.set_name("cost_fn_build_max");
    double norm_max;
    int j;
    for(i=0;i<_associates.get_dim();i++){
        for(j=0;j<_chifn->get_dim();j++){
            if(i==0 || _chifn->get_pt(_associates.get_data(i),j)<min.get_data(j)){
                min.set(j,_chifn->get_pt(_associates.get_data(i),j));
            }
            if(i==0 || _chifn->get_pt(_associates.get_data(i),j)>max.get_data(j)){
                max.set(j,_chifn->get_pt(_associates.get_data(i),j));
            }
        }
    }

    if(min.get_dim()>0){
        for(i=0;i<_chifn->get_dim();i++){
            norm.set(i,max.get_data(i)-min.get_data(i));
        }
    }
    else{
        for(i=0;i<_chifn->get_dim();i++){
            norm.set(i,1.0);
        }
    }

    for(i=0;i<_chifn->get_dim();i++){
        if(norm.get_data(i)<1.0e-20){
            norm.set(i,1.0);
        }
    }

    array_1d<double> norm_sorted;
    array_1d<int> norm_dex;
    norm_sorted.set_name("cost_fn_build_norm_sorted");
    norm_dex.set_name("cost_fn_build_norm_dex");
    for(i=0;i<norm.get_dim();i++){
        norm_dex.add(i);
    }
    sort(norm, norm_sorted, norm_dex);
    if(min_or_med==1){
        _scalar_norm=norm_sorted.get_data(norm_dex.get_dim()/2);
    }
    else{
        _scalar_norm=norm_sorted.get_data(0);
    }

    _relative_norm.set_dim(_chifn->get_dim());

    for(i=0;i<_chifn->get_dim()/4;i++){
        _relative_norm.set(norm_dex.get_data(i), 0.71);
    }
    for(;i<_chifn->get_dim()/2;i++){
        _relative_norm.set(norm_dex.get_data(i), 0.82);
    }
    for(;i<(_chifn->get_dim()*3)/4;i++){
        _relative_norm.set(norm_dex.get_data(i), 0.91);
    }
    for(;i<_chifn->get_dim();i++){
        _relative_norm.set(norm_dex.get_data(i), 1.0);
    }

    for(i=0;i<_chifn->get_dim();i++){
        if(_relative_norm.get_data(i)<0.1){
            printf("WARNING relative norm %d %e\n",i,_relative_norm.get_data(i));
            exit(1);
        }
    }

}


double cost_fn::nn_distance(const array_1d<double> &pt){
    double dd;
    int i,j,k;
    double dd_avg;
    double ct=0.0;
    dd_avg=0.0;

    for(i=0;i<_associates.get_dim();i++){
        dd=0.0;
        for(j=0;j<_chifn->get_dim();j++){
            dd+=power((pt.get_data(j)-_chifn->get_pt(_associates.get_data(i),j))/(_scalar_norm*_relative_norm.get_data(j)),2);
        }
        if(dd>1.0e-20){
            dd_avg += 1.0/sqrt(dd);
            ct+=1.0;
        }
    }
    if(ct>0.0){
        dd_avg = dd_avg/ct;
        return 1.0/dd_avg;
    }
    return 0.0;
}


int cost_fn::get_called(){
   return _called;
}

double cost_fn::operator()(const array_1d<double> &pt){

    if(_chifn==NULL){
        printf("WARNING cannot call cost_fn operator; _chifn is NULL\n");
        exit(1);
    }

    _called++;

    double mu;
    int i_found;
    _chifn->evaluate(pt,&mu,&i_found);

    if(_associates.get_dim()==0){
        return mu;
    }

    double delta=_chifn->target()-_chifn->chimin();

    double distance;

    double exp_term;

    if(mu>_chifn->target()){
        exp_term=exp((_chifn->target()-mu)/_envelope);
    }
    else{
        exp_term=1.0;
    }

    if(exp_term>1.0e-7){
        distance=nn_distance(pt);
    }
    else{
        distance=0.0;
    }

    double val;
    val = mu-1.0*delta*distance*exp_term;
    _pt_cache.add(i_found);
    _fn_cache.add(val);
    return val;
}
