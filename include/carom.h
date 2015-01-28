#ifndef CAROM_H
#define CAROM_H

#include <stdio.h>
#include <time.h>
#include "eigen_wrapper.h"
#include "wrappers.h"
#include "simplex.h"
#include "node.h"

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

    void search();
    void write_pts();
    
    int get_called();
    int active_nodes();
    
private:

    chisq_wrapper _chifn;
    arrayOfNodes _nodes;
    int _write_every,_last_written;
    int _ct_simplex,_ct_node,_calls_to_simplex;
    double _time_started;
    
    char _outname[letters],_timingname[letters];

    void simplex_search();
    void assess_node(int);

};

class gp_cost : public function_wrapper{

public:

    gp_cost();
    ~gp_cost();
    void set_chifn(chisq_wrapper*);
    virtual double operator()(array_1d<double>&);
    virtual int get_called();

private:
    chisq_wrapper *_chifn;
    
    array_2d<double> _covarin,_covar;
    array_1d<int> _neigh,_neigh_buff;
    array_1d<double> _dd,_qq;
    double _ell,_fbar;
    int _called;

    void is_it_safe(char*);

};

#endif
