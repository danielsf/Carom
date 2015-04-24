#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mcmc/mcmc.h"
#include "chisq.h"


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

ellipses_integrable chisq(6,2);

//set the maximum and minimum values in parameter space
array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");

for(i=0;i<6;i++){
    min.set(i,-20.0);
    max.set(i,20.0);
}


mcmc mcmc_test(4,seed,&chisq);
for(i=0;i<6;i++){
    mcmc_test.set_min(i,min.get_data(i));
    mcmc_test.set_max(i,max.get_data(i));
}

mcmc_test.guess_bases(1);
mcmc_test.set_burnin(2000,1000);
mcmc_test.set_name_root("chains/test_chain");
mcmc_test.sample(5000);

}
