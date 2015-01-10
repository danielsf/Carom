#ifndef WRAPPER_H
#define WRAPPER_H

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

class chisq_wrapper : function_wrapper{

public:
    chisq_wrapper();
    ~chisq_wrapper();
    
    void set_chisquared(*chisquared);
    
    virtual double operator()(array_1d<double>&);

    void set_target(double);
    void set_seed(int);
    void set_deltachi(double);
    void set_characteristic_length(int, double);
    void set_min(array_1d<double>&);
    void set_max(array_1d<double>&);
    void set_ddmin(double);

    double target();
    
    double random_double();
    int random_int();

    void evaluate(array_1d<double>&, double*, int*);

private:
    double _chimin,_deltachi,_target,_ddmin;
    int _adaptive_target,_seed;
    
    array_1d<double> _characteristic_length,_range_min,_range_max;
    
    chisquared *chifn;
    kd_tree *kptr;
    Ran *dice;
    
    int is_valid(array_1d<double>&);
}

#endif
