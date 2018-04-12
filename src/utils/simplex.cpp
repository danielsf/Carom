#include "simplex.h"

simplex_minimizer::simplex_minimizer(){
    _alpha=1.0;
    _beta=0.5;
    _gamma=2.1;

    initialize();
}

simplex_minimizer::simplex_minimizer(double aa, double bb, double gg){
    _alpha=aa;
    _beta=bb;
    _gamma=gg;

    initialize();
}

simplex_minimizer::~simplex_minimizer(){}

void simplex_minimizer::initialize(){
    _last_found_tol=1.0e-4;
    _cost=NULL;
    _chisquared=NULL;
    _dice=NULL;
    _limit=-1;
    _n_gradients=0;
    _bases.reset();
    _bases.set_name("simplex_bases");

    _min_temp=-3.0;

    _use_gradient=0;

    _abort_max_factor=10;

    _freeze_temp=-1;

    _transform.set_name("simplex_transform");
    _origin.set_name("simplex_origin");
    _min_pt.set_name("simplex_min_pt");
    _ff.set_name("simplex_ff");
    _pstar.set_name("simplex_pstar");
    _pstarstar.set_name("simplex_pstarstar");
    _pts.set_name("simplex_pts");
    _last_improved_ff.set_name("simplex_last_improved_ff");
    _last_improved_pts.set_name("simplex_last_improved_pts");

    _is_a_model=0;
}

void simplex_minimizer::set_abort_max_factor(int ii){
    _abort_max_factor=ii;
}

void simplex_minimizer::set_chisquared(function_wrapper *ff){
    _chisquared=ff;
}

void simplex_minimizer::set_cost(function_wrapper *cc){
    _cost=cc;
}

void simplex_minimizer::set_dice(Ran *dd){
    _dice=dd;
}

void simplex_minimizer::is_a_model(){
    _is_a_model=1;
}

void simplex_minimizer::use_gradient(){
    _use_gradient=1;
}

void simplex_minimizer::do_not_use_gradient(){
    _use_gradient=0;
}

void simplex_minimizer::freeze_temp(){
    _freeze_temp=1;
}

void simplex_minimizer::unfreeze_temp(){
    _freeze_temp=0;
}

void simplex_minimizer::set_minmax(array_1d<double> &min, array_1d<double> &max){
    if(min.get_dim()!=max.get_dim()){
        printf("WARNING simplex_minimizer::set_minmax min %d max %d\n",
        min.get_dim(),max.get_dim());

        exit(1);
    }

    int i;
    _origin.set_dim(min.get_dim());
    _transform.set_dim(min.get_dim());
    for(i=0;i<min.get_dim();i++){
        _origin.set(i,min.get_data(i));
        _transform.set(i,max.get_data(i)-min.get_data(i));
    }

}

void simplex_minimizer::paranoia(){

    int need_to_thaw=0;
    if(_freeze_called==0){
        _freeze_called=1;
        need_to_thaw=1;
    }

    int i;
    double mu;
    for(i=0;i<_pts.get_rows();i++){
        mu=evaluate(_pts(i));
        _ff.set(i,mu);
    }

    find_il();

    if(need_to_thaw==1){
        _freeze_called=0;
    }

}

void simplex_minimizer::find_il(){
    if(_ff.get_dim()==0 ||
       _pts.get_rows()==0 ||
       _pts.get_cols()==0 ||
       _pts.get_rows()!=_ff.get_dim()){

       printf("WARNING cannot call find_il %d %d %d\n",
       _ff.get_dim(),_pts.get_rows(),_pts.get_cols());

       exit(1);
   }

   int i;
   for(i=0;i<_ff.get_dim();i++){
       if(i==0 || _ff.get_data(i)<_ff.get_data(_il))_il=i;
       if(i==0 || _ff.get_data(i)>_ff.get_data(_ih))_ih=i;
   }

   //_min_ff=_ff.get_data(_il);
   //printf("    find_il set min %e %d\n",_min_ff,_il);
   /*double mu;
   mu=evaluate(_pts(_il)[0]);
   if(fabs(mu-_min_ff)>1.0e-2){
       printf("WARNING find_il should have found %e\n",mu);
       exit(1);
   }*/
}

double simplex_minimizer::evaluate(const array_1d<double> &pt){
    /*
    pt will be in the transformed dimensions of the simplex
    */

    if(_chisquared==NULL){
        printf("WARNING cannot call simplex_minimizer::evaluate\n");
        printf("_chisquared is NULL\n");
        exit(1);
    }

    int i,j;
    array_1d<double> vv;
    vv.set_name("simplex_minimizer_vv");
    for(i=0;i<pt.get_dim();i++){
        vv.set(i,_origin.get_data(i)+pt.get_data(i)*_transform.get_data(i));
    }

    double fval,raw;
    fval=_chisquared[0](vv);

    if(_min_pt.get_dim()==0){
        for(i=0;i<pt.get_dim();i++){
            _min_pt.set(i,vv.get_data(i));
        }
    }

    if(fval<_true_min_ff){
        _true_min_ff=fval;
        for(i=0;i<pt.get_dim();i++){
            _min_pt.set(i,vv.get_data(i));
        }
    }

    raw=fval;

    double cval;
    cval=0.0;
    if(_cost!=NULL){
        cval=evaluate_cost(vv);
        fval+=cval;
    }

   /* printf("    %e %e %d %d -- %d\n",
    fval,raw,_called_evaluate,_chisquared->get_called(),
    _chisquared->get_called()-_called_evaluate);*/

    if(fval<_min_ff){
        if(_ff_last_found-fval>fabs(_ff_last_found*_last_found_tol)){
            _last_found=_called_evaluate;
            _ff_last_found=fval;
        }
        _min_ff=fval;
        //printf("    setting min %e\n",_min_ff);
        //printf("min %e true %e cost %e raw %e\n",
        //_min_ff,_true_min_ff,cval,raw);

        if(_ff.get_dim()==pt.get_dim()+1){
            for(i=0;i<_pts.get_rows();i++){
                _last_improved_ff.set(i,_ff.get_data(i));
                for(j=0;j<_pts.get_cols();j++){
                    _last_improved_pts.set(i,j,_pts.get_data(i,j));
                }
            }
        }
    }

    if(_called_evaluate>_last_cooled_off+500 && !(_temp<_min_temp) && _freeze_temp==0){
        cool_off();
        if(_freeze_called==0)i=1;
        else i=0;

        _freeze_called=1;
        _freeze_temp=1;
        fval=evaluate(pt);
        _freeze_temp=0;
        if(i==1)_freeze_called=0;
    }


    if(_freeze_called==0)_called_evaluate++;

    return fval;
}

void simplex_minimizer::cool_off(){
    if(_freeze_temp==1 || _temp<_min_temp) return;

    _temp-=1.0;

    int i,need_thaw_temp,need_thaw_called;

    _min_ff=2.0*exception_value;

    if(_freeze_called==0)need_thaw_called=1;
    else need_thaw_called=0;

    if(_freeze_temp==0)need_thaw_temp=1;
    else need_thaw_temp=0;

    _freeze_called=1;
    _freeze_temp=1;

    for(i=0;i<_pts.get_rows();i++){
        _ff.set(i,evaluate(_pts(i)));
    }
    _fstar=evaluate(_pstar);
    _fstarstar=evaluate(_pstarstar);

    expand();

    _last_found=_called_evaluate;
    _last_cooled_off=_called_evaluate;
    if(need_thaw_called==1)_freeze_called=0;
    if(need_thaw_temp==1)_freeze_temp=0;
    /*printf("cooled off temp %e min %e true %e\n",_temp,_min_ff,_true_min_ff);
    for(i=0;i<_pts.get_rows();i++){
        printf("    %e\n",_ff.get_data(i));
    }*/

}

double simplex_minimizer::evaluate_cost(const array_1d<double> &vv){
    /*
    vv will already be transformed into natural coordinates
    */

    if(_cost==NULL){
        printf("WARNING cannot call simplex_minimizer::evaluate_cost\n");
        printf("_cost is NULL\n");
        exit(1);
    }

    if(_temp<_min_temp)return 0.0;

    double cval;
    cval=_cost[0](vv);

    if(_freeze_temp==0)_called_cost++;

    return exp(_temp)*cval;
}

void simplex_minimizer::is_it_safe(char *word){

    if(_chisquared==NULL){
        printf("WARNING in simplex::%s -- _chisquared is NULL\n",word);
        exit(1);
    }

    if(_cost!=NULL){
        if(_dice==NULL){
            printf("WARNING in simplex::%s -- _cost is not NULL, but _chisquared is\n",word);
            exit(1);
        }
    }
}

double simplex_minimizer::get_minimum(){
    return _true_min_ff;
}

void simplex_minimizer::find_minimum(array_2d<double> &seed, array_1d<double> &min_pt){
    array_1d<double> errors,seed_pt;
    errors.set_name("mapping_errors");
    seed_pt.set_name("mapping_seed");
    double mu,mu_min;
    int i_min;
    int i,j,k;
    double max, min;
    for(i=0;i<seed.get_rows();i++){
        mu=_chisquared[0](seed(i));
        if(i==0 || mu<mu_min){
            mu_min=mu;
            i_min=i;
        }
    }
    for(i=0;i<seed.get_cols();i++){
        seed_pt.set(i,seed.get_data(i_min,i));
    }

    for(i=0;i<seed.get_cols();i++){
        for(j=0;j<seed.get_rows();j++){
            if(_bases.get_rows()>0){
                mu=0.0;
                for(k=0;k<seed.get_cols();k++){
                    mu+=seed.get_data(j,k)*_bases.get_data(i,k);
                }
            }
            else{
                mu=seed.get_data(j,i);
            }
            if(j==0 || mu<min){
                min=mu;
            }
            if(j==0 ||mu>max){
                max=mu;
            }
        }
        errors.set(i,0.5*(max-min));
    }

    find_minimum(seed_pt,errors,min_pt);

}

void simplex_minimizer::find_minimum(array_1d<double> &seed_pt,
                                     array_1d<double> &errors,
                                     array_1d<double> &min_pt){

    is_it_safe("find_minimum");

    MnUserParameters input_params;
    char word[100];
    _minuit_fn.set_dim(seed_pt.get_dim());
    _minuit_fn.set_chisquared(_chisquared);
    int i,j;
    double mu;
    mu=_chisquared[0](seed_pt);
    printf("    minimizer starts with %e %d\n",mu,_bases.get_rows());

    if(_bases.get_rows()>0){
        _minuit_fn.set_bases(_bases);
        for(i=0;i<seed_pt.get_dim();i++){
            mu=0.0;
            for(j=0;j<seed_pt.get_dim();j++){
                mu+=seed_pt.get_data(j)*_bases.get_data(i,j);
            }
            sprintf(word,"pp%d",i);
            input_params.Add(word,mu,errors.get_data(i));
        }
    }
    else{
        for(i=0;i<seed_pt.get_dim();i++){
            sprintf(word,"pp%d",i);
            input_params.Add(word,seed_pt.get_data(i),errors.get_data(i));
        }

    }

    MnMigrad migrad_runner(_minuit_fn, input_params);
    FunctionMinimum minuit_result = migrad_runner();
    std::vector<double> pp_out = minuit_result.UserParameters().Params();
    if(_bases.get_rows()==0){
        for(i=0;i<seed_pt.get_dim();i++){
            min_pt.set(i,pp_out[i]);
            _min_pt.set(i,pp_out[i]);
        }
    }
    else{
        for(i=0;i<seed_pt.get_dim();i++){
            min_pt.set(i,0.0);
            _min_pt.set(i,0.0);
            for(j=0;j<seed_pt.get_dim();j++){
                min_pt.add_val(i,pp_out[j]*_bases.get_data(j,i));
                _min_pt.add_val(i,pp_out[j]*_bases.get_data(j,i));
            }
        }
    }
    printf("    leaving find minimum -- %d -- %e\n",min_pt.get_dim(),minuit_result.Fval());
}


double simplex_minimizer::get_dx(int dex){

    double max,min;
    int i;
    for(i=0;i<_pts.get_rows();i++){
        if(i==0 || _pts.get_data(i,dex)<min){
            min=_pts.get_data(i,dex);
        }
        if(i==0 || _pts.get_data(i,dex)>max){
            max=_pts.get_data(i,dex);
        }
    }

    if(max-min<1.0e-20){
        return 1.0;
    }

    return max-min;

}


void simplex_minimizer::calculate_gradient(const array_1d<double> &pp, array_1d<double> &grad){

    array_1d<double> trial;
    trial.set_name("simplex_calculate_gradient_trial");
    int i;
    for(i=0;i<pp.get_dim();i++){
        trial.set(i,pp.get_data(i));
    }

    double dx=0.01;
    double y1,y2,x1,x2,xnorm;
    for(i=0;i<pp.get_dim();i++){
        xnorm=get_dx(i);
        trial.add_val(i,dx*xnorm);
        y1=evaluate(trial);
        x1=trial.get_data(i);
        trial.subtract_val(i,2.0*dx*xnorm);
        y2=evaluate(trial);
        x2=trial.get_data(i);
        grad.set(i,(y1-y2)/(x1-x2));
        trial.set(i,pp.get_data(i));
    }

}

void simplex_minimizer::gradient_cloud(){
    if(_dice==NULL){
        printf("WARNING cannot use simplex gradient; _dice is NULL\n");
        exit(1);
    }
    find_il();

    double true_min_0=_true_min_ff;

    int need_thaw_temp,need_thaw_called;

    if(_freeze_temp==0)need_thaw_temp=1;
    else need_thaw_temp=0;

    if(_freeze_called==0)need_thaw_called=1;
    else need_thaw_called=0;

    _freeze_temp=1;
    _freeze_called=1;

    array_1d<double> gradient;
    gradient.set_name("simplex_gradient");

    double step,dd;
    int i,j,k;
    step=-1.0;

    for(i=0;i<_pts.get_rows();i++){
        for(j=i+1;j<_pts.get_rows();j++){
            dd=0.0;
            for(k=0;k<_pts.get_cols();k++){
                dd+=power(_pts.get_data(i,k)-_pts.get_data(j,k),2);
            }
            dd=sqrt(dd);
            if(dd>step){
                step=dd;
            }
        }
    }

    int ix;
    for(ix=0;ix<_pts.get_rows();ix++){

        calculate_gradient(_pts(ix), gradient);
        gradient.normalize();
        for(i=0;i<_pts.get_cols();i++){
            gradient.multiply_val(i,-1.0);
        }
        for(i=0;i<_pts.get_cols();i++){
            _pts.add_val(ix,i,0.5*step*gradient.get_data(i));
        }

    }

    double mu;

    for(i=0;i<_pts.get_rows();i++){
        mu=evaluate(_pts(i));
        _ff.set(i,mu);
    }
    find_il();
    if(need_thaw_temp==1)_freeze_temp=0;
    if(need_thaw_called==1)_freeze_called=0;
    _last_called_gradient=_called_evaluate;

    printf("    delta from gradient_cloud %e\n",_true_min_ff-true_min_0);

}


void simplex_minimizer::gradient_minimizer(){
    if(_dice==NULL){
        printf("WARNING cannot use simplex gradient; _dice is NULL\n");
        exit(1);
    }
    _n_gradients++;
    find_il();

    int need_thaw_temp,need_thaw_called;

    if(_freeze_temp==0)need_thaw_temp=1;
    else need_thaw_temp=0;

    if(_freeze_called==0)need_thaw_called=1;
    else need_thaw_called=0;

    _freeze_temp=1;
    _freeze_called=1;

    array_1d<double> gradient;
    gradient.set_name("simplex_gradient");

    double step,dd;
    int i,j,k;
    step=-1.0;

    for(i=0;i<_pts.get_rows();i++){
        for(j=i+1;j<_pts.get_rows();j++){
            dd=0.0;
            for(k=0;k<_pts.get_cols();k++){
                dd+=power(_pts.get_data(i,k)-_pts.get_data(j,k),2);
            }
            dd=sqrt(dd);
            if(dd>step){
                step=dd;
            }
        }
    }

    calculate_gradient(_pts(_il), gradient);
    gradient.normalize();
    for(i=0;i<_pts.get_cols();i++){
        gradient.multiply_val(i,-1.0);
    }


    array_1d<double> perturbation;
    double mu;
    for(i=0;i<_pts.get_rows();i++){
        if(i!=_il){
            for(j=0;j<_pts.get_cols();j++){
                perturbation.set(j,normal_deviate(_dice,0.0,1.0));
            }
            perturbation.normalize();
            for(j=0;j<_pts.get_cols();j++){
                _pts.add_val(i,j,step*0.1*gradient.get_data(j)+0.05*step*perturbation.get_data(j));
            }
        }
        else{
            for(j=0;j<_pts.get_cols();j++){
                _pts.add_val(i,j,step*0.5*gradient.get_data(j));
            }
        }
        mu=evaluate(_pts(i));
        _ff.set(i,mu);
    }
    find_il();
    if(need_thaw_temp==1)_freeze_temp=0;
    if(need_thaw_called==1)_freeze_called=0;
    _last_called_gradient=_called_evaluate;

}


void simplex_minimizer::expand(){
    if(_dice==NULL){
        printf("WARNING simplex cannot expand; it does not have a _dice value\n");
        exit(1);
    }
    if(_pts.get_rows()==0) return;

    //printf("expanding\n");

    int need_thaw_temp,need_thaw_called;

    if(_freeze_temp==0)need_thaw_temp=1;
    else need_thaw_temp=0;

    if(_freeze_called==0)need_thaw_called=1;
    else need_thaw_called=0;

    _freeze_temp=1;
    _freeze_called=1;

    int i,j;
    array_1d<double> span_min,span_max;
    span_min.set_name("simplex_expand_span_min");
    span_max.set_name("simplex_expand_span_max");
    for(i=0;i<_pts.get_rows();i++){
        for(j=0;j<_pts.get_cols();j++){
            if(i==0 || _pts.get_data(i,j)<span_min.get_data(j)){
                span_min.set(j,_pts.get_data(i,j));
            }
            if(i==0 || _pts.get_data(i,j)>span_max.get_data(j)){
                span_max.set(j,_pts.get_data(i,j));
            }
        }
    }

    double mu;
    find_il();
    for(i=0;i<_pts.get_rows();i++){
        if(i!=_il){
            for(j=0;j<_pts.get_cols();j++){
                _pts.set(i,j,_pts.get_data(_il,j)+4.0*(_dice->doub()-0.5)*(span_max.get_data(j)-span_min.get_data(j)));
            }

        }
        mu=evaluate(_pts(i));
        if(i==_il && fabs(mu-_ff.get_data(_il))>1.0e-2){
            printf("WARNING at _il %e %e -- %e\n",mu,_ff.get_data(_il),_temp);
            //exit(1);
        }
        _ff.set(i,mu);
    }
    find_il();
    if(need_thaw_temp==1)_freeze_temp=0;
    if(need_thaw_called==1)_freeze_called=0;
}


void simplex_minimizer::get_pt(int ii, array_1d<double> &out){
    int i;
    for(i=0;i<_pts.get_cols();i++){
        out.set(i,_pts.get_data(ii,i)*_transform.get_data(i)+_origin.get_data(i));
    }
}

void simplex_minimizer::get_minpt(array_1d<double> &out){
    int i;
    for(i=0;i<_min_pt.get_dim();i++){
        out.set(i,_min_pt.get_data(i));
    }
}
