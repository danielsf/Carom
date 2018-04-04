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

    _bases.set_name("dchi_interior_bases");
    _projected_associates.set_name("dchi_interior_projected_associates");
    _pt_projected.set_name("dchi_interior_pt_projected");
    _pt_cache.set_name("dchi_interior_pt_cache");
    _fn_cache.set_name("dchi_interior_fn_cache");
    _associates.set_name("dchi_interior_fn_associates");
    _relative_norm.set_name("dchi_interior_relative_norm");

    _projected_associates.reset();
    _associates.reset_preserving_room();
    _fn_cache.reset_preserving_room();
    _pt_cache.reset_preserving_room();

    _chifn=cc;

    int i;

    for(i=0;i<aa.get_dim();i++){
        _associates.add(aa.get_data(i));
    }

    _min_or_med=min_or_med;

    if(_bases.get_rows()>0 && _chifn!=NULL){
        _set_d_params();
    }
}


void cost_fn::_set_d_params(){
    printf("in set d params\n");
    _project_associates();
    _set_scalar_norm();
    printf("projected\n");
}


void cost_fn::_set_scalar_norm(){
    if(_chifn==NULL){
        printf("CANNOT set scalar norm; chifn is NULL\n");
        exit(1);
    }

    _relative_norm.reset_preserving_room();

    array_1d<double> norm_sorted;
    array_1d<int> norm_dex;
    norm_sorted.set_name("set_scalar_norm_sorted");
    norm_dex.set_name("set_scalar_norm_dex");

    double harmonic_mean=0.0;
    double geometric_mean=0.0;
    double harmonic_ct=0.0;
    double min,max,mu;
    int idim,ipt,j;
    for(idim=0;idim<_chifn->get_dim();idim++){
        min=2.0*exception_value;
        max=-2.0*exception_value;
        for(ipt=0;ipt<_chifn->get_pts();ipt++){
            if(_chifn->get_fn(ipt)<_chifn->target()){
                mu=0.0;
                for(j=0;j<_chifn->get_dim();j++){
                   mu+=_chifn->get_pt(ipt,j)*_bases.get_data(idim,j);
                }
                if(mu<min){
                    min=mu;
                }
                if(mu>max){
                    max=mu;
                }
            }
        }

        if(max-min>1.0e-20){
            harmonic_ct+=1.0;
            harmonic_mean+=1.0/(max-min);
            geometric_mean+=log(max-min);
            _relative_norm.set(idim, max-min);
        }
        else{
            _relative_norm.set(idim,1.0);
        }
        norm_dex.set(idim,idim);
    }

    sort(_relative_norm, norm_sorted, norm_dex);

    if(harmonic_ct>0.1){
        harmonic_mean = harmonic_ct/harmonic_mean;
        geometric_mean=geometric_mean/harmonic_ct;
        geometric_mean=exp(geometric_mean);
    }
    double median = norm_sorted.get_data(norm_dex.get_dim()/2);
    _scalar_norm = 5.0;
    printf("    set scalar norm to %e\n",_scalar_norm);
    printf("    median %e geometric %e harmonic %e\n",
    median,geometric_mean,harmonic_mean);

    for(idim=0;idim<_relative_norm.get_dim();idim++){
        if(_relative_norm.get_data(idim)>_scalar_norm){
            _relative_norm.set(idim,_scalar_norm);
        }
    }

}

void cost_fn::_project_associates(){
    int ipt,idim,j;
    _projected_associates.reset_preserving_room();
    _projected_associates.set_dim(_associates.get_dim(),
                                  _chifn->get_dim());

    for(ipt=0;ipt<_associates.get_dim();ipt++){
        for(idim=0;idim<_chifn->get_dim();idim++){
            _projected_associates.set(ipt,idim,0.0);
            for(j=0;j<_chifn->get_dim();j++){
                _projected_associates.add_val(ipt,idim,
                        _chifn->get_pt(_associates.get_data(ipt),
                                       j)*_bases.get_data(idim,j));
            }
        }
    }
}

double cost_fn::nn_distance(const array_1d<double> &pt){
    double dd;
    int i,j,k;
    double dd_avg;
    double ct=0.0;
    dd_avg=0.0;

    if(_projected_associates.get_rows()==0){
        _set_d_params();
    }

    for(i=0;i<_associates.get_dim();i++){
        dd=0.0;
        for(j=0;j<_chifn->get_dim();j++){
            dd+=power((pt.get_data(j)-_projected_associates.get_data(i,j))/_relative_norm.get_data(j),2);
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

    _pt_projected.reset_preserving_room();
    _pt_projected.set_dim(_chifn->get_dim());
    int i,j;
    for(i=0;i<_chifn->get_dim();i++){
        _pt_projected.set(i,0.0);
        for(j=0;j<_chifn->get_dim();j++){
            _pt_projected.add_val(i,pt.get_data(j)*_bases.get_data(i,j));
        }
    }

    if(exp_term>1.0e-7){
        distance=nn_distance(_pt_projected);
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
