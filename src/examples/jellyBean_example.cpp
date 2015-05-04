#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "carom.h"

int main(int iargc, char *argv[]){

//d=8 -> delta_chisq=15.5
//d=5 -> delta_chisq=11

int i,j;
int seed=99;
int dim=5;
int nsamples=100000;
double delta_chisq=11.0;

if(iargc>1)seed=atoi(argv[1]);
if(iargc>2)dim=atoi(argv[2]);
if(iargc>3)delta_chisq=atof(argv[3]);

if(iargc>4){
    nsamples=atoi(argv[4]);
}

if(seed<0){
    seed=int(time(NULL));
    if(seed>10000)seed=seed%10000;
}

char timingname[letters],outname[letters];

//what is the name of the file where APS will store its timing information
sprintf(timingname,"output/jellyBean_d%d_s%d_timing.sav",dim,seed);

//what is the name of the file where APS will output the points it sampled
sprintf(outname,"output/jellyBean_d%d_s%d_output.sav",dim,seed);

printf("seed %d\n",seed);

//declare the chisquared function APS will be searching
//ellipses_integrable chisq(dim,ncenters);

jellyBean chisq(dim,1.0,20.0);

//declare APS
//the '20' below is the number of nearest neighbors to use when seeding the
//Gaussian process
//
//the '11.0' is the \Delta\chi^2 corresponding to a 95% confidence limit
//on a 5-dimensional parameter space

carom carom_test;
carom_test.set_target(delta_chisq);
carom_test.set_seed(seed);

//pass chisq to the aps object
carom_test.set_chisquared(&chisq);

//how often will APS stop and write its output
carom_test.set_write_every(3000);

//set the maximum and minimum values in parameter space
array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");
max.set_dim(dim);
min.set_dim(dim);

for(i=0;i<dim;i++){
    min.set(i,-30.0);
    max.set(i,30.0);
}


carom_test.set_timingname(timingname);
carom_test.set_outname(outname);
carom_test.set_min(min);
carom_test.set_max(max);

//initialize aps with 1000 random samples
printf("time to initialize\n");
chisq.reset_timer();
carom_test.initialize(1000);
int active_nodes=1;
carom_test.search(nsamples);

}
