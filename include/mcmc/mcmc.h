#ifndef MCMC_H
#define MCMC_H

#include "mcmc/chain.h"
#include "chisq.h"

class mcmc{

public:

    ~mcmc()
    mcmc(int,int);
    void set_chisq(chisquared&);

private:

    arrayOfChains _chains;
    chisquared *chisq;

    array_2d<double> _covar,_bases;
    array_1d<double> _sigma;

    double thinby(chain, int, array_1d<double>&, array_1d<double>&, array_2d<double>&);

};

#endif
