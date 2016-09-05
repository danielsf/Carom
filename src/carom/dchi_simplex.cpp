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


dchi_interior_simplex::dchi_interior_simplex(chisq_wrapper *cc, array_1d<int> &aa){

    _called=0;

    _mask.set_name("dchi_interior_mask");
    _median_associate.set_name("dchi_interior_median");
    _norm.set_name("dchi_interior_norm");
    _bases.set_name("dchi_interior_bases");
    _hyper_center.set_name("dchi_interior_hyper_center");
    _hyper_norm.set_name("dchi_interior_hyper_norm");

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
        if(i==0 || norm.get_data(i)<_scalar_norm){
            _scalar_norm=norm.get_data(i);
        }
    }

    set_bases();
    _set_hyper_ellipse();
    calibrate_model();

}

void dchi_interior_simplex::calibrate_model(){

    array_1d<double> ddsq;
    ddsq.set_name("dchi_calibrate_ddsq");

    array_1d<double> min,max;
    min.set_name("dchi_calibrate_min");
    max.set_name("dchi_calibrate_max");

    double mu;
    int i,j,k;
    for(i=0;i<_associates.get_dim();i++){
        for(j=0;j<_chifn->get_dim();j++){

            mu=0.0;

            for(k=0;k<_chifn->get_dim();k++){
                mu+=_chifn->get_pt(_associates.get_data(i),k)*_bases.get_data(j,k);
            }

            if(i==0 || mu<min.get_data(j)){
                min.set(j,mu);
            }

            if(i==0 || mu>max.get_data(j)){
                max.set(j,mu);
            }
        }
    }

    for(i=0;i<_chifn->get_dim();i++){
        _norm.set(i,max.get_data(i)-min.get_data(i));
    }

    double dd,ddsqsum=0.0;;
    int ip;
    for(i=0;i<_associates.get_dim();i++){
        dd=0.0;
        ip=_associates.get_data(i);
        for(j=0;j<_chifn->get_dim();j++){
            mu=0.0;
            for(k=0;k<_chifn->get_dim();k++){
                mu+=(_chifn->get_pt(_chifn->mindex(),k)-_chifn->get_pt(ip,k))*_bases.get_data(j,k);
            }
            dd+=power(mu/_norm.get_data(j),2);
        }
        ddsq.set(i,dd);
        ddsqsum+=dd;
    }

    _alpha=0.0;
    for(i=0;i<_associates.get_dim();i++){
        ip=_associates.get_data(i);
        _alpha+=ddsq.get_data(i)*(_chifn->get_fn(ip)-_chifn->chimin());
    }

    _alpha=_alpha/ddsqsum;
}


double dchi_interior_simplex::apply_model(array_1d<double> &pt){

    double ddsq,mu;
    int i,j;
    ddsq=0.0;
    for(i=0;i<_chifn->get_dim();i++){
        mu=0.0;
        for(j=0;j<_chifn->get_dim();j++){
            mu+=(_chifn->get_pt(_chifn->mindex(),j)-pt.get_data(j))*_bases.get_data(i,j);
        }
        ddsq+=power(mu/_norm.get_data(i),2);
    }

    return _chifn->chimin()+_alpha*ddsq;
}


double dchi_interior_simplex::nn_distance(array_1d<double> &pt){
    double dd;
    int i,j;

    if(_just_median==1){
        printf("cannot use median distance\n");
        exit(1);
        dd=0.0;
        for(i=0;i<_chifn->get_dim();i++){
            dd+=power((pt.get_data(i)-_median_associate.get_data(i))/_scalar_norm,2);
        }
        return sqrt(dd);
    }

    if(_associates.get_dim()==0){
        return 0.0;
    }
    double dd_min=2.0*exception_value;

    double delta=_chifn->target()-_chifn->chimin();

    for(i=0;i<_associates.get_dim();i++){
        if(_mask.get_dim()==0 || _mask.get_data(i)==1){
            dd=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                dd+=power((pt.get_data(j)-_chifn->get_pt(_associates.get_data(i),j))/_scalar_norm,2);
            }
            if(dd<dd_min){
                dd_min=dd;
            }
        }
    }

    return delta*sqrt(dd_min);
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

    if(mu>_chifn->target()){
        exp_term=exp((_chifn->target()-mu)/_envelope);
    }
    else{
        exp_term=1.0;
    }

    return mu-1.0*distance*exp_term;
}


void dchi_interior_simplex::_principal_set_bases(){

    printf("\nprincipal bases\n\n");

    _bases.reset_preserving_room();

    array_1d<double> dir,dir_best;
    dir.set_name("dchi_princ_bases_dir");
    dir_best.set_name("dchi_princ_bases_dir_best");

    int ip,i,j,ia;
    double dd,dd_best,component;

    while(_bases.get_rows()!=_chifn->get_dim()){
        for(ip=0;ip<_associates.get_dim();ip++){
            ia=_associates.get_data(ip);
            for(i=0;i<_chifn->get_dim();i++){
                dir.set(i,_chifn->get_pt(ia,i)-_chifn->get_pt(_chifn->mindex(),i));
            }
            for(i=0;i<_bases.get_rows();i++){
                component=0.0;
                for(j=0;j<_chifn->get_dim();j++){
                    component+=dir.get_data(j)*_bases.get_data(i,j);
                }
                for(j=0;j<_chifn->get_dim();j++){
                    dir.subtract_val(j,component*_bases.get_data(i,j));
                }
            }

            dd=dir.normalize();
            if(ip==0 || dd>dd_best){
                dd_best=dd;
                for(i=0;i<_chifn->get_dim();i++){
                    dir_best.set(i,dir.get_data(i));
                }
            }
        }

        _bases.add_row(dir_best);
    }


    for(i=0;i<_chifn->get_dim();i++){
        for(j=i;j<_chifn->get_dim();j++){
            component=0.0;
            for(ia=0;ia<_chifn->get_dim();ia++){
                component+=_bases.get_data(i,ia)*_bases.get_data(j,ia);
            }


            if(i==j){
                if(fabs(component-1.0)>0.001){
                    printf("WARNING basis %d norm %e\n",i,component);
                    exit(1);
                }
            }
            else{
                if(fabs(component)>0.001){
                    printf("WARNING dot product between bases %d %d is %e\n",
                    i,j,component);
                    exit(1);
                }
            }
        }
    }

}


void dchi_interior_simplex::_random_set_bases(){

    array_1d<double> vv;
    vv.set_name("dchi_set_bases_vv");
    int i,j;
    double component;

    _bases.reset_preserving_room();
    while(_bases.get_rows()<_chifn->get_dim()){
        for(i=0;i<_chifn->get_dim();i++){
            vv.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        for(i=0;i<_bases.get_rows();i++){
            component=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                component+=vv.get_data(j)*_bases.get_data(i,j);
            }
            for(j=0;j<_chifn->get_dim();j++){
                vv.subtract_val(j,component*_bases.get_data(i,j));
            }
        }
        component=vv.normalize();
        if(component>1.0e-10){
            _bases.add_row(vv);
        }
    }

}


void dchi_interior_simplex::_set_hyper_ellipse(){

    _hyper_center.reset_preserving_room();
    _hyper_norm.reset_preserving_room();

    double component;
    array_1d<double> cc,cc_sorted;
    array_1d<int> cc_dex;
    cc.set_name("hyper_ellipse_cc");
    cc_sorted.set_name("hyper_ellipse_cc_sorted");
    cc_dex.set_name("hyper_ellipse_cc_dex");

    double x1,x2;
    int ix,i,j;
    for(ix=0;ix<_bases.get_rows();ix++){
        cc.reset_preserving_room();
        cc_sorted.reset_preserving_room();
        cc_dex.reset_preserving_room();

        for(i=0;i<_associates.get_dim();i++){
            component=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                component+=_chifn->get_pt(_associates.get_data(i),j)*_bases.get_data(ix,j);
            }
            cc.set(i,component);
            cc_dex.set(i,i);
        }
        sort(cc,cc_sorted,cc_dex);

        i=cc_dex.get_dim();

        x1=cc_sorted.get_data(1/6);
        x2=cc_sorted.get_data((5*i)/6);
        _hyper_center.set(ix,0.5*(x1+x2));
        _hyper_norm.set(ix,0.5*(x2-x1));
    }

}


double dchi_interior_simplex::_hyper_ellipse_distance(array_1d<double> &pt){
    double component;
    int i,j;
    double dd=0.0;
    for(i=0;i<_chifn->get_dim();i++){
        component=0.0;
        for(j=0;j<_chifn->get_dim();j++){
            component+=pt.get_data(j)*_bases.get_data(i,j);
        }
        dd+=power((component-_hyper_center.get_data(i))/_hyper_norm.get_data(i),2);
    }
    return sqrt(dd);
}
