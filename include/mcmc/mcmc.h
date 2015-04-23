#ifndef MCMC_H
#define MCMC_H

#include "mcmc/chain.h"
#include "chisq.h"
#include "simplex.h"
#include "eigen_wrapper.h"

class mcmc{

public:

    ~mcmc();
    mcmc(int, chisquared*);
    void set_burn_in(int);
    void set_covar_lim(double);
    void set_name_root(char*);
    void set_min(int,double);
    void set_max(int,double);

private:

    arrayOfChains _chains;
    chisquared *_chisq;
    Ran *dice;

    array_2d<double> _bases;
    array_1d<double> _sigma,_centerpt,_guess_max,_guess_min;

    char _name_root[letters];

    int _burn_in,_last_set;
    double _max_covar,_covar_lim;

    void initialize();
    void guess_bases(array_2d<double>&);
    void find_fisher_matrix(array_2d<double>&, array_1d<double>&);
    void bisection(array_1d<double>&, double, array_1d<double>&, array_1d<double>&);

};

#endif
