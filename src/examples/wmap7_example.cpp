#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "carom.h"
#include "wmap_likelihood_function.h"


int main(int iargc, char *argv[]){

//d=8 -> delta_chisq=15.5
//d=5 -> delta_chisq=11

int i,j;
int seed=99;
int dim,ncenters;
int nsamples=50000;
int dof=1199;

if(iargc>1)seed=atoi(argv[1]);

if(iargc>2){
    nsamples=atoi(argv[2]);
}

if(seed<0){
    seed=int(time(NULL));
    if(seed>10000)seed=seed%10000;
}

char timingname[letters],outname[letters];

//what is the name of the file where APS will store its timing information
sprintf(timingname,"output/wmap7_s%d_timing.sav",seed);

//what is the name of the file where APS will output the points it sampled
sprintf(outname,"output/wmap7_s%d_output.sav",seed);

printf("seed %d\n",seed);

//declare the chisquared function APS will be searching
//ellipses_integrable chisq(dim,ncenters);

wmap_likelihood chisq;

carom carom_test;
//carom_test.set_deltachi(12.6);
carom_test.set_target(1280.7);
carom_test.set_seed(seed);
carom_test.set_dof(dof);
carom_test.set_confidence_limit(0.95);

//pass chisq to the aps object
carom_test.set_chisquared(&chisq);

//how often will APS stop and write its output
carom_test.set_write_every(100);

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

carom_test.set_timingname(timingname);
carom_test.set_outname(outname);
carom_test.set_min(min);
carom_test.set_max(max);

//initialize aps with 1000 random samples
printf("time to initialize\n");
chisq.reset_timer();
carom_test.initialize(100);
int active_nodes=1;
carom_test.search(nsamples);

}
