#ifndef CAROM_H
#define CAROM_H

#include <stdio.h>
#include <time.h>
#include "wrappers.h"

class carom{

public:

    carom();
    ~carom();
    
    void initialize(int);
    
    void set_seed(int);
    void set_min(array_1d<double>&);
    void set_max(array_1d<double>&);
    void set_characteristic_length(int,double);
    void set_deltachi(double);
    void set_target(double);
    void set_write_every(int);
    void set_outname(char*);
    void set_timingname(char*);
    
    void set_chisquared(chisquared*);

    void write_pts();
    
private:

    chisq_wrapper chifn;
    int _write_every,_last_written;
    double _time_started;
    
    char _outname[letters],_timingname[letters];

};

#endif
