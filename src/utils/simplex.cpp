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
    _cost=NULL;
    _chisquared=NULL;
    _dice=NULL;

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

void simplex_minimizer::use_gradient(){
    _use_gradient=1;
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

double simplex_minimizer::evaluate(array_1d<double> &pt){
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
        if(_min_ff-fval>fabs(_min_ff*1.0e-4)){
            _last_found=_called_evaluate;
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
        _ff.set(i,evaluate(_pts(i)[0]));
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

double simplex_minimizer::evaluate_cost(array_1d<double> &vv){
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
    is_it_safe("find_minimum");

    if(seed.get_rows()!=seed.get_cols()+1){
        printf("WARNING you gave simplex minimizer %d points of %d dim\n",
        seed.get_rows(),seed.get_cols());
        printf("Need dim+1 points\n");
        exit(1);
    }


    if(_cost!=NULL){
        _temp=0.0;
    }
    else{
        _temp=-1000.0;
    }

    int ibefore=_chisquared->get_called();

    if(_freeze_temp<0)_freeze_temp=0;

    _freeze_called=0;
    _called_cost=0;
    _last_found=0;
    _called_evaluate=0;
    _last_called_gradient=0;
    _last_cooled_off=0;


    _min_ff=2.0*exception_value;
    _true_min_ff=2.0*exception_value;
    _min_pt.reset();
    _pstar.reset();
    _pstarstar.reset();
    _last_improved_ff.reset();
    _pts.reset();
    _last_improved_pts.reset();

    int i;
    if(_origin.get_dim()==0){
        for(i=0;i<seed.get_cols();i++){
            _origin.set(i,0.0);
            _transform.set(i,1.0);
        }
    }

    int j;
    _pts.set_cols(seed.get_cols());
    _last_improved_pts.set_cols(seed.get_cols());
    for(i=0;i<seed.get_rows();i++){
        for(j=0;j<seed.get_cols();j++){
            _pts.set(i,j,(seed.get_data(i,j)-_origin.get_data(j))/_transform.get_data(j));
        }
    }

    array_1d<double> initial_min,initial_max;
    initial_min.set_name("simplex_initial_min");
    initial_max.set_name("simplex_initial_max");
    for(i=0;i<_pts.get_rows();i++){
        for(j=0;j<_pts.get_cols();j++){
            if(i==0 || _pts.get_data(i,j)<initial_min.get_data(j)){
                initial_min.set(j,_pts.get_data(i,j));
            }
            if(i==0 || _pts.get_data(i,j)>initial_max.get_data(j)){
                initial_max.set(j,_pts.get_data(i,j));
            }
        }
    }

    double mu;
    for(i=0;i<seed.get_rows();i++){
        mu=evaluate(_pts(i)[0]);
        _ff.set(i,mu);
    }

    find_il();

    int abort_max=_abort_max_factor*seed.get_cols();
    if(abort_max<100){
        abort_max=100;
    }
    int dim=seed.get_cols();
    double spread,gradient_threshold;

    array_1d<double> pbar,buffer;
    pbar.set_name("simplex_pbar");
    buffer.set_name("simplex_buffer");

    spread=_ff.get_data(_ih)-_ff.get_data(_il);

    while(_called_evaluate-_last_found<abort_max){
       for(i=0;i<dim;i++){
           pbar.set(i,0.0);
           for(j=0;j<dim+1;j++){
               if(j!=_ih){
                   pbar.add_val(i,_pts.get_data(j,i));
               }
           }
           pbar.divide_val(i,double(dim));
       }

       for(i=0;i<dim;i++){
           _pstar.set(i,(1.0+_alpha)*pbar.get_data(i)-_alpha*_pts.get_data(_ih,i));
       }
       _fstar=evaluate(_pstar);

       if(_fstar<_ff.get_data(_ih) && _fstar>_ff.get_data(_il)){
           _ff.set(_ih,_fstar);
           for(i=0;i<dim;i++){
               _pts.set(_ih,i,_pstar.get_data(i));
           }
       }
       else if(_fstar<_ff.get_data(_il)){
           for(i=0;i<dim;i++){
               _pstarstar.set(i,_gamma*_pstar.get_data(i)+(1.0-_gamma)*pbar.get_data(i));
           }
           _fstarstar=evaluate(_pstarstar);

           if(_fstarstar<_ff.get_data(_il)){
               for(i=0;i<dim;i++){
                   _pts.set(_ih,i,_pstarstar.get_data(i));
               }
               _ff.set(_ih,_fstarstar);
           }
           else{
               for(i=0;i<dim;i++){
                   _pts.set(_ih,i,_pstar.get_data(i));
               }
               _ff.set(_ih,_fstar);
           }
       }

       find_il();

       j=1;
       for(i=0;i<dim+1;i++){
           if(_fstar<_ff.get_data(i) && i!=_ih){
               j=0;
           }
       }

       if(j==1){
           for(i=0;i<dim;i++){
               _pstarstar.set(i,_beta*_pts.get_data(_ih,i)+(1.0-_beta)*pbar.get_data(i));
           }
           _fstarstar=evaluate(_pstarstar);

           if(_fstarstar<_ff.get_data(_ih)){
               for(i=0;i<dim;i++)_pts.set(_ih,i,_pstarstar.get_data(i));
               _ff.set(_ih,_fstarstar);
           }
           else{
               find_il();
               for(i=0;i<dim+1;i++){
                   if(i!=_il){
                       for(j=0;j<dim;j++){
                           buffer.set(j,0.5*(_pts.get_data(i,j)+_pts.get_data(_il,j)));
                       }
                       mu=evaluate(buffer);
                       _ff.set(i,mu);
                       for(j=0;j<dim;j++){
                           _pts.set(i,j,buffer.get_data(j));
                       }
                   }
               }
           }
       }

       find_il();
       spread=_ff.get_data(_ih)-_ff.get_data(_il);

       if(_temp<_min_temp)gradient_threshold=1.0e-1;
       else gradient_threshold=0.1*_min_ff;

       if(spread<gradient_threshold &&
           _use_gradient==1 &&
           _called_evaluate>abort_max/2+_last_called_gradient){
           gradient_minimizer();
       }

       if(_called_evaluate-_last_found>=abort_max && !(_temp<_min_temp)){
           expand();
           _last_found=_called_evaluate;
       }
       //printf("min %.3e treu %.3e spread %.3e since last %d %.3e\n",_min_ff,_true_min_ff,spread,_last_found-_called_evaluate,_temp);
       //printf("spread %e %e %e\n\n",spread,_temp,_min_ff);
    }

    for(i=0;i<dim;i++){
        min_pt.set(i,_min_pt.get_data(i));
    }
    printf("    leaving simplex %d %d %d\n",_called_evaluate,_last_found,abort_max);
    printf("    temp %e\n",_temp);
    printf("    actually called %d\n",_chisquared->get_called()-ibefore);
    printf("    _true_min_ff %e\n",_true_min_ff);

    _freeze_temp=-1;
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


void simplex_minimizer::calculate_gradient(array_1d<double> &pp, array_1d<double> &grad){

    array_1d<double> trial;
    trial.set_name("simplex_calculate_gradient_trial");
    int i;
    for(i=0;i<pp.get_dim();i++){
        trial.set(i,pp.get_data(i));
    }

    double dx=0.01;
    double y1,y2,xnorm;
    y1=evaluate(trial);
    for(i=0;i<pp.get_dim();i++){
        xnorm=get_dx(i);
        trial.add_val(i,dx*xnorm);
        y2=evaluate(trial);
        grad.set(i,(y1-y2)/(pp.get_data(i)-trial.get_data(i)));
        trial.set(i,pp.get_data(i));
    }

}


void simplex_minimizer::gradient_minimizer(){
    if(_dice==NULL){
        printf("WARNING cannot use simplex gradient; _dice is NULL\n");
        exit(1);
    }
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

    calculate_gradient(_pts(_il)[0], gradient);
    gradient.normalize();
    for(i=0;i<_pts.get_cols();i++){
        gradient.multiply_val(i,-1.0);
    }


    array_1d<double> perturbation;
    double mu;
    for(i=0;i<_pts.get_rows();i++){
        for(j=0;j<_pts.get_cols();j++){
            perturbation.set(j,normal_deviate(_dice,0.0,1.0));
        }
        perturbation.normalize();
        for(j=0;j<_pts.get_cols();j++){
            _pts.add_val(i,j,step*0.1*gradient.get_data(j)+0.05*step*perturbation.get_data(j));
        }
        mu=evaluate(_pts(i)[0]);
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
        mu=evaluate(_pts(i)[0]);
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

void simplex_minimizer::get_minpt(array_1d<double> &out){
    int i;
    for(i=0;i<_min_pt.get_dim();i++){
        out.set(i,_min_pt.get_data(i));
    }
}
