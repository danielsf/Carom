#ifndef MCMC_H
#define MCMC_H

#include "mcmc/chain.h"
#include "chisq.h"
#include "simplex.h"
#include "eigen_wrapper.h"

class mcmc{

public:

    ~mcmc();
    mcmc();
    mcmc(int, int, function_wrapper*);
    void initialize(int, int, function_wrapper*);
    void set_burnin(int);
    void set_name_root(char*);
    void set_min(int,double);
    void set_max(int,double);
    void guess_bases(array_1d<double>&,double,int);
    void guess_bases(double,int);
    double acceptance_rate();
    double update_bases();
    void write_burnin();
    void write_chains();
    void write_timing(int);
    void write_timing(char*);

    void sample(int);

private:

    arrayOfChains _chains;
    function_wrapper *_chisq;
    Ran *_dice;

    double _factor,_time_started;

    array_2d<double> _bases;
    array_1d<double> _sigma,_guess_max,_guess_min;

    char _name_root[letters];

    int _burn_in;

    void initialize();
    void validate_bases();
    void find_fisher_eigen(array_1d<double>&, array_2d<double>&);
    double find_minimum_point(array_1d<double>&);
    void find_fisher_matrix(array_1d<double>&, array_2d<double>&);
    void mcmc_bisection(double, array_1d<double>&, double, array_1d<double>&, array_1d<double>&);

};

#endif
