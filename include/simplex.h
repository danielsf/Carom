#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "containers.h"
#include "goto_tools.h"
#include "wrappers.h"

#include <vector>
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/FunctionMinimum.h"

using namespace ROOT::Minuit2;

class _minuit_wrapper : public FCNBase{
    private:
        function_wrapper *_chisquared;
        int _dim;
        mutable array_1d<double> _pp;
        array_2d<double> _bases;

    public:

    _minuit_wrapper(){
        _chisquared=NULL;
        _dim=-1;
        _pp.set_name("_minuit_wrapper_pp");
        _bases.set_name("_minuit_wrapper_bases");
    }

    ~_minuit_wrapper(){}

    void set_bases(array_2d<double> &b_in){
        int i,j;
        _bases.reset();
        _bases.set_dim(b_in.get_rows(), b_in.get_cols());
        for(i=0;i<b_in.get_rows();i++){
            for(j=0;j<b_in.get_cols();j++){
                _bases.set(i,j,b_in.get_data(i,j));
            }
        }
    }

    void set_chisquared(function_wrapper *cc){
        _chisquared=cc;
    }

    void set_dim(int ii){
        _dim=ii;
        _pp.reset_preserving_room();
        _pp.set_dim(ii);
    }

    virtual double Up() const {return 1.0;}

    double operator()(const std::vector<double> &pp) const{
        if(_dim<0){
            printf("CANNOT call _minuit_wrapper operator; dim<0\n");
            exit(1);
        }

        if(_chisquared==NULL){
            printf("CANNOT call _minuit_wrapper operator; _chisquared is NULL\n");
            exit(1);
        }

        int i,j;

        if(_bases.get_rows()==0){
            for(i=0;i<_dim;i++){
                _pp.set(i,pp[i]);
            }
        }
        else{
            for(i=0;i<_dim;i++){
                _pp.set(i,0.0);
                for(j=0;j<_dim;j++){
                    _pp.add_val(i,pp[j]*_bases.get_data(j,i));
                }
            }
        }

        return _chisquared[0](_pp);
    }

};


class simplex_minimizer{

public:

    simplex_minimizer();
    simplex_minimizer(double,double,double);
    ~simplex_minimizer();

    void set_chisquared(function_wrapper*);
    void set_cost(function_wrapper*);
    void set_dice(Ran*);
    void set_minmax(array_1d<double>&, array_1d<double>&);
    void set_abort_max_factor(int);
    void use_gradient();
    void do_not_use_gradient();
    void freeze_temp();
    void unfreeze_temp();
    void get_minpt(array_1d<double>&);
    void get_pt(int,array_1d<double>&);
    void is_a_model();

    /*
    seed_pt
    errors
    min_pt (output)
    */
    void find_minimum(array_1d<double>&,
                      array_1d<double>&,
                      array_1d<double>&);

    void find_minimum(array_2d<double>&, array_1d<double>&);
    double get_minimum();

    void set_limit(int ii){
        _limit=ii;
    }

    void set_bases(array_2d<double> &b_in){
         int i,j;
         _bases.reset();
         _bases.set_dim(b_in.get_rows(), b_in.get_cols());
         for(i=0;i<b_in.get_rows();i++){
             for(j=0;j<b_in.get_cols();j++){
                 _bases.set(i,j,b_in.get_data(i,j));
             }
         }
    }

private:

    double evaluate(const array_1d<double>&);
    double evaluate_cost(const array_1d<double>&);
    void cool_off();

    void gradient_minimizer();
    void gradient_cloud();
    void calculate_gradient(const array_1d<double>&, array_1d<double>&);
    double get_dx(int);
    void expand();

    void find_il();
    void paranoia();

    void initialize();

    void is_it_safe(char*);

    double _temp,_min_ff,_true_min_ff,_fstar,_fstarstar,_min_temp;
    double _alpha,_beta,_gamma;
    int _il,_ih,_called_cost,_freeze_temp,_use_gradient;
    int _freeze_called,_last_called_gradient;
    int _last_found,_called_evaluate,_abort_max_factor;
    double _ff_last_found;
    double _last_found_tol;
    int _last_cooled_off;
    int _is_a_model;
    int _limit;
    int _n_gradients;
    array_1d<double> _transform, _origin,_ff,_pstar,_pstarstar,_min_pt;
    array_1d<double> _last_improved_ff;
    array_1d<double> _min,_max;
    array_2d<double> _pts,_last_improved_pts;
    array_2d<double> _bases;

    Ran *_dice;
    function_wrapper *_chisquared, *_cost;
    _minuit_wrapper _minuit_fn;

    /*
    cost will need to be a function_wrapper sub-class
    that has pointers to all of the aps nodes.
    */

};

#endif
