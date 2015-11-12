#include "dchi_simplex.h"


dchi_simplex_base::dchi_simplex_base(chisq_wrapper *chi_in,
                                     array_1d<int> &associates_in){

    _called=0;
    _chisq=chi_in;
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
            for(j=0;j<_chisq->get_dim();j++){
                mu=_chisq->get_pt(i,j);
                if(i==0 || mu<min.get_data(j)){
                    min.set(j,mu);
                }

                if(i==0 || mu>max.get_data(j)){
                    max.set(j,mu);
                }
            }
        }

        for(i=0;i<_chisq->get_dim();i++){
            _norm.set(i,max.get_data(i)-min.get_data(i));
        }
    }

}

int dchi_simplex_base::get_called(){
    return _called;
}

double dchi_simplex_base::associate_distance(array_1d<double> &pt){

     if(_associates.get_dim()<=0){
         printf("Cannot call dchi_simplex_base::associate_distance\n");
         printf("_associates.get_dim() %d\n",_associates.get_dim());
         exit(1);
     }

     double dd,ddmin;
     int i,j;
     ddmin=0.0;
     for(i=0;i<_associates.get_dim();i++){
         dd=0.0;
         for(j=0;j<_chisq->get_dim();j++){
            dd+=power((pt.get_data(j)-_chisq->get_pt(_associates.get_data(i),j))/_norm.get_data(j),2);
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
    _chisq->evaluate(pt,&mu,&i_found);

    double dmu=fabs(_chisq->target()-mu);
    double delta=_chisq->target()-_chisq->chimin();

    _called++;
    return dmu-exp(-0.1*dmu/delta)*delta*distance*2.0;
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
    _chisq->evaluate(pt,&mu,&i_found);

    if(_associates.get_dim()==0){
        return mu;
    }

    double delta=_chisq->target()-_chisq->chimin();
    double dmu=fabs(mu-_chisq->chimin());

    double exp_term=exp(-dmu/delta);
    double distance;
    if(exp_term>1.0e-5){
        distance=associate_distance(pt);
    }
    else{
        return mu;
    }

    return mu-exp_term*distance*delta*2.0;
}
