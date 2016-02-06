#include "gp_lin.h"

void gp_lin::build(array_2d<double> &pt_in, array_1d<double> &fn_in,
               array_1d<double> &min_in, array_1d<double> &max_in){

    if(_kd!=NULL){
        printf("WARNING cannot build _kd; is not NULL\n");
        exit(1);
    }

    if(_fn!=NULL){
        printf("WARNING cannot build _fn; is not NULL\n");
        exit(1);
    }

    _kd=new kd_tree;
    _fn=new array_1d<double>;

    _fn->set_name("gp_lin_fn");

    int i;
    for(i=0;i<fn_in.get_dim();i++){
        _fn->set(i,fn_in.get_data(i));
    }

    _kd->build_tree(pt_in,min_in,max_in);

    _built_here=1;

}

void gp_lin::add_pt(array_1d<double> &pt, double ff){
    if(_built_here==0){
        printf("cannot call gp_lin::add_pt; kd tree was not built here\n");
        exit(1);
    }

    if(_kd==NULL || _fn==NULL){
        printf("cannot call gp_lin::add_pt; null pointers\n");
        exit(1);
    }
    _kd->add(pt);
    _fn->add(ff);
}

double gp_lin::operator()(array_1d<double> &pt){

    array_1d<int> local_dex;
    local_dex.set_name("gp_operator_local_dex");
    array_1d<double> local_distance;
    local_distance.set_name("gp_operator_local_distance");

    _kd->nn_srch(pt, _nn, local_dex, local_distance);

    array_1d<double> dir,mu,trial;
    dir.set_name("operator_dir");
    mu.set_name("operator_mu");
    trial.set_name("operator_trial");

    int i,ix1,ix2;
    int p1,p2;
    double d1,d2;
    double dd,dmu,slope,aa,denom;


    array_1d<double> distance,distance_sorted;
    array_1d<int> distance_dex;
    distance.set_name("operator_distance");
    distance_sorted.set_name("operator_distance_sorted");
    distance_dex.set_name("operator_distance_dex");

    array_1d<double> mu_sorted;
    array_1d<int> mu_dex;

    for(ix1=0;ix1<local_dex.get_dim();ix1++){
        p1=local_dex.get_data(ix1);
        for(ix2=ix1+1;ix2<local_dex.get_dim();ix2++){
            p2=local_dex.get_data(ix2);
            for(i=0;i<pt.get_dim();i++){
                dir.set(i,_kd->get_pt(p2,i)-_kd->get_pt(p1,i));
            }
            dd=dir.normalize();
            dmu=_fn->get_data(p2)-_fn->get_data(p1);
            slope=dmu/dd;
            aa=0.0;
            denom=0.0;
            for(i=0;i<pt.get_dim();i++){
                denom+=dir.get_data(i)*dir.get_data(i);
                aa+=dir.get_data(i)*(pt.get_data(i)-_kd->get_pt(p1,i));
            }
            aa=aa/denom;
            mu.add(_fn->get_data(p1)+slope*aa);
            mu_dex.add(mu.get_dim()-1);
            dd=0.0;
            for(i=0;i<pt.get_dim();i++){
                trial.set(i,_kd->get_pt(p1,i)+aa*dir.get_data(i));
            }
            dd=_kd->distance(pt,trial);
            /*dd+=_kd->distance(p1,p2);
            d1=_kd->distance(p1,pt);
            d2=_kd->distance(p2,pt);
            if(d1<d2){
                dd+=d1;
            }
            else{
                dd+=d2;
            }*/
            distance.add(dd);
            distance_dex.add(distance.get_dim()-1);

        }
    }

    sort_and_check(distance,distance_sorted,distance_dex);
    _ell=_ell_factor*(distance_sorted.get_data(distance.get_dim()/2)-distance_sorted.get_data(0));

    sort_and_check(mu,mu_sorted,mu_dex);

    double total_wgt=0.0;
    array_1d<double> wgt;
    double ww;
    for(i=0;i<distance.get_dim();i++){
        ww=exp(-0.5*(distance.get_data(i)-distance_sorted.get_data(0))/_ell);
        total_wgt+=ww;
        wgt.add(ww);
    }


    double sum=0.0;
    for(i=0;i<mu.get_dim()-1 && sum<total_wgt*0.5;i++){
        sum+=wgt.get_data(mu_dex.get_data(i));
    }

    return mu_sorted.get_data(i);

}

void gp_lin::optimize(array_2d<double> &pts, array_1d<double> &ff){

    double log_ell,ell_best;
    double cost,cost_best,ln10;
    int i;
    double worst_dmu,worst_ff,worst_mu,dmu,mu;
    array_1d<double> cost_arr,cost_sorted;
    array_1d<int> cost_dex;
    cost_best=2.0*exception_value;
    ell_best=-1.0;
    ln10=log(10.0);
    worst_dmu=-2.0*exception_value;
    for(log_ell=-1.0;log_ell<3.0;log_ell+=0.25){
        cost=0.0;
        _ell_factor=exp(ln10*log_ell);
        for(i=0;i<pts.get_rows();i++){
            mu=this[0](pts(i)[0]);
            if(fabs(ff.get_data(i))<fabs(mu)){
                dmu=fabs((mu-ff.get_data(i))/ff.get_data(i));
            }
            else{
                dmu=fabs((mu-ff.get_data(i))/mu);
            }
            cost_arr.set(i,dmu);
            cost_dex.set(i,i);
            if(dmu>worst_dmu){
                worst_dmu=dmu;
                worst_mu=mu;
                worst_ff=ff.get_data(i);
            }
             //printf("%e %e\n\n",ff.get_data(i),mu);
        }
        //printf("        worst %e %e %e\n",worst_ff,worst_mu,worst_dmu);
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
        //exit(1);
    }

    printf("    best cost %e ell %e\n",cost_best,ell_best);
    _ell_factor=ell_best;

}

