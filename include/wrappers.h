#ifndef WRAPPERS_H
#define WRAPPERS_H

#include <math.h>
#include <time.h>
#include <stdio.h>

#include "containers.h"

class function_wrapper{
public:
    function_wrapper();
    ~function_wrapper();
    virtual double operator()(array_1d<double>&);
    virtual int get_called();
};

#endif
