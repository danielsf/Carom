#include "dchi_simplex.h"


dchi_simplex_base::dchi_simplex_base(chisq_wrapper *chi_in,
                                     array_1d<int> &associates_in){

    _called=0;
    _chifn=chi_in;
    _min_0=_chifn->chimin();
    _associates.set_name("dchi_simplex_associates");
    int i,j;
    for(i=0;i<associates_in.get_dim();i++){
        _associates.add(associates_in.get_data(i));
    }

    array_1d<double> min,max;
    min.set_name("dchi_base_constructor_min");
    max.set_name("dchi_base_constructor_max");
    _norm.set_name("dchi_base_norm");

    double mu;
    if(_associates.get_dim()>0){
        for(i=0;i<_associates.get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                mu=_chifn->get_pt(i,j);
                if(i==0 || mu<min.get_data(j)){
                    min.set(j,mu);
                }

                if(i==0 || mu>max.get_data(j)){
                    max.set(j,mu);
                }
            }
        }

        for(i=0;i<_chifn->get_dim();i++){
            if(max.get_data(i)-min.get_data(i)<1.0e-20){
                _norm.set(i,1.0);
            }
            else{
                _norm.set(i,0.5*(max.get_data(i)-min.get_data(i)));
            }
        }
    }

}

int dchi_simplex_base::get_called(){
    return _called;
}

double dchi_simplex_base::associate_distance(array_1d<double> &pt){

     if(_associates.get_dim()<=0){
         return 0.0;
     }

     double dd,ddmin;
     int i,j;
     ddmin=0.0;
     for(i=0;i<_associates.get_dim();i++){
         dd=0.0;
         for(j=0;j<_chifn->get_dim();j++){
            dd+=power((pt.get_data(j)-_chifn->get_pt(_associates.get_data(i),j))/_norm.get_data(j),2);
         }
         if(i==0 || dd<ddmin){
             ddmin=dd;
         }
     }
     return sqrt(ddmin);
}

double dchi_simplex_base::operator()(array_1d<double> &pt){
    printf("WARNING calling dchi_simplex_base operator; shouldn't do that\n");
    exit(1);
}

dchi_boundary_simplex::dchi_boundary_simplex(chisq_wrapper *cc, array_1d<int> &aa) :
                       dchi_simplex_base(cc, aa){}

double dchi_boundary_simplex::operator()(array_1d<double> &pt){
    double distance=associate_distance(pt);
    double mu;
    int i_found;
    _chifn->evaluate(pt,&mu,&i_found);

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

dchi_multimodal_simplex::dchi_multimodal_simplex(chisq_wrapper *cc, array_1d<int> &aa) :
                         dchi_simplex_base(cc, aa){}

void dchi_multimodal_simplex::set_norm(array_1d<double> &nn){

    if(nn.get_dim()!=_norm.get_dim()){
        printf("WARNING cannot set norm in dchi_multimodal_simplex\n");
        printf("input dim %d need %d\n",nn.get_dim(),_norm.get_dim());
        exit(1);
    }

    int i;
    for(i=0;i<_norm.get_dim();i++){
        _norm.set(i,nn.get_data(i));
    }
}

double dchi_multimodal_simplex::operator()(array_1d<double> &pt){

    _called++;

    double mu;
    int i_found;
    _chifn->evaluate(pt,&mu,&i_found);

    if(_associates.get_dim()==0){
        return mu;
    }

    double delta=_chifn->target()-_chifn->chimin();
    double dmu=fabs(mu-_chifn->target());

    double exp_term;

    if(_chifn->target()<mu){
        exp_term=exp(-0.1*dmu/delta);
    }
    else{
        exp_term=1.0;
    }

    double distance;
    if(exp_term>1.0e-5){
        distance=associate_distance(pt);
    }
    else{
        return mu;
    }

    return mu-exp_term*distance*delta*2.0;
}

dchi_interior_simplex::dchi_interior_simplex(chisq_wrapper *cc, array_1d<int> &aa){

    _called=0;

    _mask.set_name("dchi_interior_mask");
    _median_associate.set_name("dchi_interior_median");

    _just_median=0;
    _chifn=cc;
    _envelope=1.0;

    int i;

    for(i=0;i<aa.get_dim();i++){
        _associates.add(aa.get_data(i));
    }

    array_1d<double> norm;
    for(i=0;i<_chifn->get_dim();i++){
        norm.set(i,_chifn->get_characteristic_length(i));
    }
    array_1d<double> min,max;
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
            if(max.get_data(i)-min.get_data(i)<norm.get_data(i)){
                norm.set(i,max.get_data(i)-min.get_data(i));
            }
        }
    }

    for(i=0;i<_chifn->get_dim();i++){
        _median_associate.set(i,0.5*(max.get_data(i)+min.get_data(i)));
        if(i==0 || norm.get_data(i)<_norm){
            _norm=norm.get_data(i);
        }
    }

}

double dchi_interior_simplex::nn_distance(array_1d<double> &pt){
    double dd;
    int i,j;

    if(_just_median==1){
        dd=0.0;
        for(i=0;i<_chifn->get_dim();i++){
            dd+=power((pt.get_data(i)-_median_associate.get_data(i))/_norm,2);
        }
        return sqrt(dd);
    }

    if(_associates.get_dim()==0){
        return 0.0;
    }
    double dd_min=2.0*exception_value;

    for(i=0;i<_associates.get_dim();i++){
        if(_mask.get_dim()==0 || _mask.get_data(i)==1){
            dd=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                dd+=power((pt.get_data(j)-_chifn->get_pt(_associates.get_data(i),j))/_norm,2);
            }
            if(dd<dd_min){
                dd_min=dd;
            }
        }
    }

    return sqrt(dd_min);
}


int dchi_interior_simplex::get_called(){
   return _called;
}

double dchi_interior_simplex::operator()(array_1d<double> &pt){

    _called++;

    double mu;
    int i_found;
    _chifn->evaluate(pt,&mu,&i_found);

    if(_associates.get_dim()==0){
        return mu;
    }

    double delta=_chifn->target()-_chifn->chimin();

    double distance=nn_distance(pt);

    double exp_term;
    if(mu<_chifn->target()){
        exp_term=1.0;
    }
    else{
        exp_term=exp((_chifn->target()-mu)/_envelope);
    }

    return mu-1.0*distance*delta*exp_term;
}
