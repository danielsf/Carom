#include "cost_fn.h"

cost_fn::cost_fn(chisq_wrapper *cc, array_1d<int> &aa){

    _called=0;

    _median_associate.set_name("dchi_interior_median");
    _bases.set_name("dchi_interior_bases");
    _norm.set_name("dchi_interior_norm");
    _pt_cache.set_name("dchi_interior_pt_cache");
    _fn_cache.set_name("dchi_interior_fn_cache");

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
    }

    array_1d<double> norm_sorted;
    array_1d<int> norm_dex;
    for(i=0;i<norm.get_dim();i++){
        norm_dex.add(i);
    }
    sort(norm, norm_sorted, norm_dex);
    _scalar_norm=norm_sorted.get_data(norm_dex.get_dim()/2);

    _set_bases();
    _set_norm();

}


double cost_fn::nn_distance(const array_1d<double> &pt){
    double dd;
    int i,j,k;

    int n_dd=3;

    if(_associates.get_dim()<n_dd){
        return 0.0;
    }

    array_1d<double> dd_arr,dd_arr_sorted;
    array_1d<int> dd_dex;
    dd_arr.set_name("cost_fn_dd_arr");
    dd_arr_sorted.set_name("cost_fn_dd_arr_sorted");
    dd_dex.set_name("cost_fn_dd_dex");
    for(i=0;i<_associates.get_dim();i++){
        dd=0.0;
        for(j=0;j<_chifn->get_dim();j++){
            dd+=power((pt.get_data(j)-_chifn->get_pt(_associates.get_data(i),j))/_scalar_norm,2);
        }
        if(dd_arr.get_dim()<n_dd){
            dd_arr.add(dd);
            dd_dex.add(_associates.get_data(i));
            if(dd_arr.get_dim()==n_dd){
                sort(dd_arr, dd_arr_sorted, dd_dex);
            }
        }
        else{
            if(dd<dd_arr_sorted.get_data(n_dd-1)){
                for(j=0;j<n_dd-1 && dd>dd_arr_sorted.get_data(j);j++);
                for(k=n_dd-1;k>j;k--){
                    dd_arr_sorted.set(k,dd_arr_sorted.get_data(k-1));
                    dd_dex.set(k,dd_dex.get_data(k-1));
                }
                dd_arr_sorted.set(j,dd);
                dd_dex.set(j,_associates.get_data(i));
            }
        }
    }

    array_1d<double> avg_pt;
    avg_pt.set_name("cost_fn_avg_pt");
    for(i=0;i<_chifn->get_dim();i++){
        avg_pt.set(i,0.0);
    }
    for(i=0;i<n_dd;i++){
        for(j=0;j<_chifn->get_dim();j++){
            avg_pt.add_val(j,_chifn->get_pt(dd_dex.get_data(i),j));
        }
    }
    for(i=0;i<_chifn->get_dim();i++){
        avg_pt.divide_val(i,double(n_dd));
    }
    double ans=0.0;
    for(i=0;i<_chifn->get_dim();i++){
        ans+=power((pt.get_data(i)-avg_pt.get_data(i))/_scalar_norm,2);
    }
    ans=sqrt(ans);

    return ans;
}


int cost_fn::get_called(){
   return _called;
}

double cost_fn::operator()(const array_1d<double> &pt){

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


    double val;
    val = mu-1.0*delta*distance*exp_term;
    _pt_cache.add(i_found);
    _fn_cache.add(val);
    return val;
}


void cost_fn::_principal_set_bases(){

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
                    printf("associates %d\n",_associates.get_dim());
                    exit(1);
                }
            }
        }
    }

}


void cost_fn::_random_set_bases(){

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


void cost_fn::_set_norm(){

    _norm.reset_preserving_room();

    double component;
    array_1d<double> cc,cc_sorted;
    array_1d<int> cc_dex;
    cc.set_name("hyper_ellipse_cc");
    cc_sorted.set_name("hyper_ellipse_cc_sorted");
    cc_dex.set_name("hyper_ellipse_cc_dex");

    int ix1,ix2;
    ix1=0;
    ix2=(5*_associates.get_dim())/6;

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

        x1=cc_sorted.get_data(ix1);
        x2=cc_sorted.get_data(ix2);
        _norm.set(ix,0.5*(x2-x1));
    }

}
