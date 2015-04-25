#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mcmc/mcmc.h"
#include "wmap_likelihood_function.h"


int main(int iargc, char *argv[]){

//d=8 -> delta_chisq=15.5
//d=5 -> delta_chisq=11

int i,j;
int seed=99;
int nsamples=50000;

if(iargc>1)seed=atoi(argv[1]);

if(iargc>2){
    nsamples=atoi(argv[2]);
}

if(seed<0){
    seed=int(time(NULL));
    if(seed>10000)seed=seed%10000;
}


printf("seed %d\n",seed);

//declare the chisquared function APS will be searching
//ellipses_integrable chisq(dim,ncenters);

wmap_likelihood chisq;

//set the maximum and minimum values in parameter space
array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");

min.set(0,0.01);
max.set(0,0.04);

min.set(1,0.01);
max.set(1,0.3);

min.set(2,0.4);
max.set(2,1.0);

min.set(3,0.005);
max.set(3,0.15);

min.set(4,0.7);
max.set(4,1.3);

min.set(5,2.0);
max.set(5,4.0);

for(i=0;i<6;i++){
    chisq.set_max(i,max.get_data(i));
    chisq.set_min(i,min.get_data(i));
}

mcmc mcmc_test(4,seed,&chisq);
for(i=0;i<6;i++){
    mcmc_test.set_min(i,min.get_data(i));
    mcmc_test.set_max(i,max.get_data(i));
}

mcmc_test.guess_bases(12.6,1);
mcmc_test.set_burnin(2000,1000);
mcmc_test.set_name_root("chains/mcmc_test_150425");
mcmc_test.sample(20000);

}
