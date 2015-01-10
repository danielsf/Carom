#ifndef WRAPPERS_H
#define WRAPPERS_H

#include <math.h>
#include <time.h>
#include <stdio.h>

#include "goto_tools.h"
#include "containers.h"
#include "kd.h"
#include "chisq.h"

class function_wrapper{
public:
    function_wrapper();
    ~function_wrapper();
    virtual double operator()(array_1d<double>&);
    virtual int get_called();
};

class chisq_wrapper : public function_wrapper{

public:
    chisq_wrapper();
    ~chisq_wrapper();
    
    void initialize(int);
    
    void set_chisquared(chisquared*);

    void set_target(double);
    void set_seed(int);
    void set_deltachi(double);
    void set_characteristic_length(int, double);
    void set_min(array_1d<double>&);
    void set_max(array_1d<double>&);
    void set_ddmin(double);
    
    double target();
    int get_pts();
    int get_dim();
    virtual int get_called();
    
    double random_double();
    int random_int();

    void evaluate(array_1d<double>&, double*, int*);
    double get_fn(int);
    double get_pt(int,int);
    
    void nn_srch(array_1d<double>&,int,array_1d<int>&,array_1d<double>&);

private:
    double _chimin,_deltachi,_target,_ddmin;
    int _adaptive_target,_seed,_called,_mindex;
    
    array_1d<double> _characteristic_length,_range_min,_range_max,_fn;
    
    chisquared *_chifn;
    kd_tree *_kptr;
    Ran *_dice;
    
    int is_valid(array_1d<double>&, int*);
    void is_it_safe(char*);
    
    array_1d<double> _valid_dd;
    array_1d<int> _valid_neigh;
};

#endif
