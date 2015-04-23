#ifndef MCMC_H
#define MCMC_H

#include "mcmc/chain.h"
#include "chisq.h"
#include "simplex.h"
#include "eigen_wrapper.h"

class mcmc{

public:

    ~mcmc();
    mcmc(int, int, chisquared*);
    void set_burn_in(int);
    void set_covar_lim(double);
    void set_name_root(char*);
    void set_min(int,double);
    void set_max(int,double);
    void guess_bases(int);
    double acceptance_rate();

private:

    arrayOfChains _chains;
    chisquared *_chisq;
    Ran *_dice;

    array_2d<double> _bases;
    array_1d<double> _sigma,_guess_max,_guess_min;

    char _name_root[letters];

    int _burn_in,_last_set;

    void initialize();
    void find_fisher_eigen(array_2d<double>&, array_1d<double>&, double*);
    void find_fisher_matrix(array_2d<double>&, array_1d<double>&, double*);
    void bisection(array_1d<double>&, double, array_1d<double>&, array_1d<double>&);

};

#endif
