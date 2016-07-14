#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "jellyBean.h"
#include "exampleLikelihoods.h"
#include "maps.h"

int main(int iargc, char *argv[]){

//d=8 -> delta_chisq=15.5
//d=5 -> delta_chisq=11

int chisq_dex=0;

int i,j;
int init=1000;
int seed=99;
int dim=4;
int nsamples=-1;
double width=1.0;
double delta_chisq=-1.0;
double abs_target=-1.0;
double confidence_limit=0.95;

char timingname[letters],outname[letters];

//what is the name of the file where APS will store its timing information
sprintf(timingname,"output/jellyBean_d%d_s%d_timing.sav",dim,seed);

//what is the name of the file where APS will output the points it sampled
sprintf(outname,"output/jellyBean_d%d_s%d_output.sav",dim,seed);

for(i=1;i<iargc;i++){
    if(argv[i][0]=='-'){
        switch(argv[i][1]){
            case 'h':
                printf("s = seed\nc = delta_chi\na = abs target\nn = nsamples\n");
                printf("d = dim\no = outputname\nt = timingname\n");
                printf("p = confidence limit\ni= n_init\n");
                exit(1);
                break;
            case 'i':
                i++;
                init=atoi(argv[i]);
                break;
            case 'x':
                i++;
                chisq_dex=atoi(argv[i]);
                break;
            case 'p':
                i++;
                confidence_limit=atof(argv[i]);
                break;
            case 's':
                i++;
                seed=atoi(argv[i]);
                break;
            case 'c':
                i++;
                delta_chisq=atof(argv[i]);
                break;
            case 'a':
                i++;
                abs_target=atof(argv[i]);
                break;
            case 'n':
                i++;
                nsamples=atoi(argv[i]);
                break;
            case 'd':
                i++;
                dim=atoi(argv[i]);
                break;
            case 'w':
                i++;
                width=atof(argv[i]);
                break;
            case 'o':
                i++;
                for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                    outname[j]=argv[i][j];
                }
                outname[j]=0;
                break;
            case 't':
                i++;
                for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                    timingname[j]=argv[i][j];
                }
                timingname[j]=0;
                break;

        }
    }
}


if(seed<0){
    seed=int(time(NULL));
    if(seed>10000)seed=seed%10000;
}


printf("seed %d\n",seed);

//declare the chisquared function APS will be searching
//ellipses_integrable chisq(dim,ncenters);

//jellyBeanData chisq(dim,1,width,100,0.4,0.4,0.02,20.0);

jellyBeanData *chisq;

if(chisq_dex==0){
    if(dim==4){
        chisq=new gaussianJellyBean4;
        printf("not going to let you do 4d non-integrable\n");
        exit(1);
    }
    else if(dim==12){
        chisq=new gaussianJellyBean12;
    }
    else if(dim==24){
        chisq=new gaussianJellyBean24;
    }
    else{
        printf("WARNING do not know what to do with dim %d\n",dim);
        exit(1);
    }
}
else if(chisq_dex==1){
    if(dim==4){
        chisq=new integrableJellyBean;
    }
    else if(dim==12){
        printf("creating new class\n");
        chisq=new integrableJellyBean12;
    }
    else{
        printf("WARNING do not have integrable jellyBean for dim %d\n",dim);
    }
}
else{
    printf("WARNING do not know what to do with chisq_dex %d\n",chisq_dex);
    exit(1);
}

printf("done constructing chisq\n");

//chisq->print_mins();

//declare APS
//the '20' below is the number of nearest neighbors to use when seeding the
//Gaussian process
//
//the '11.0' is the \Delta\chi^2 corresponding to a 95% confidence limit
//on a 5-dimensional parameter space

maps carom_test;
if(delta_chisq>0.0 && abs_target>0.0){
    printf("WARNING cannot set delta chisq and target\n");
    exit(1);
}
else if(delta_chisq<0.0 && abs_target<0.0){
     printf("WARNING did not set either delta chisq or target\n");
     exit(1);
}

if(delta_chisq>0.0){
    carom_test.set_deltachi(delta_chisq);
}
else{
    carom_test.set_target(abs_target);
}
carom_test.set_seed(seed);
carom_test.set_confidence_limit(confidence_limit);

//pass chisq to the aps object
carom_test.set_chisquared(chisq);

//how often will APS stop and write its output
carom_test.set_write_every(3000);

//set the maximum and minimum values in parameter space
array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");
max.set_dim(chisq->get_dim());
min.set_dim(chisq->get_dim());

for(i=0;i<chisq->get_dim();i++){
    min.set(i,-40.0);
    max.set(i,40.0);
}

/*min.set(1,0.0);
max.set(1,80.0);
min.set(2,-15.0);
max.set(2,65.0);
min.set(3,0.0);
max.set(3,80.0);*/

carom_test.set_timingname(timingname);
carom_test.set_outname(outname);
carom_test.set_min(min);
carom_test.set_max(max);

//initialize aps with 1000 random samples
printf("time to initialize\n");
chisq->reset_timer();
carom_test.initialize(init);
int active_nodes=1;
printf("ready to search\n");
carom_test.search(nsamples);
//chisq->print_mins();

array_1d<double> v0,v1;
chisq->get_basis(0,v0);
chisq->get_basis(1,v1);
printf("\nfirst two bases\n");
for(i=0;i<chisq->get_dim();i++){
    printf("%e %e\n",v0.get_data(i),v1.get_data(i));
}

printf("\nwidths\n");
for(i=0;i<chisq->get_dim();i++){
    printf("%e\n",chisq->get_width(0,i));
}

}
