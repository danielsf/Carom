#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mcmc/mcmc.h"
#include "exampleLikelihoods.h"


int main(int iargc, char *argv[]){

//d=8 -> delta_chisq=15.5
//d=5 -> delta_chisq=11

int i,j;
int doGuess=0;
double chiLimit=12.6;
int seed=99;
int nsamples=50000;
int burnin=3000;
int dim=4;
int nChains=4;
char nameroot[letters];
char logname[letters];

logname[0]=0;

for(i=1;i<iargc;i++){
    if(argv[i][0]=='-'){
        switch(argv[i][1]){
            case 'h':
                printf("s = seed\nb = burnin\nn = nsamples\n");
                printf("d = dim\no = nameroot\nc = nChains\n");
                printf("g = doGuess (limit)\n");
                exit(1);
            case 's':
                i++;
                seed=atoi(argv[i]);
                break;
            case 'b':
                i++;
                burnin=atoi(argv[i]);
                break;
            case 'n':
                i++;
                nsamples=atoi(argv[i]);
                break;
            case 'd':
                i++;
                dim=atoi(argv[i]);
                break;
            case 'c':
                i++;
                nChains=atoi(argv[i]);
                break;
            case 'g':
                doGuess=1;
                i++;
                chiLimit=atof(argv[i]);
                break;
            case 'o':
                i++;
                for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                    nameroot[j]=argv[i][j];
                }
                nameroot[j]=0;
                break;
            case 'l':
                i++;
                for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                    logname[j]=argv[i][j];
                }
                logname[j]=0;

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

jellyBeanData *chisq;

if(dim==4){
    chisq=new gaussianJellyBean4;
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

chisq->enable_logging();

//set the maximum and minimum values in parameter space
array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");

for(i=0;i<dim;i++){
    min.set(i,-40.0);
    max.set(i,40.0);
}


mcmc mcmc_test(nChains,seed,chisq);
for(i=0;i<dim;i++){
    mcmc_test.set_min(i,min.get_data(i));
    mcmc_test.set_max(i,max.get_data(i));
}

array_1d<double> guess_pt;
if(doGuess==1){
    /*guess_pt.set(0,-1.280762331898815987e+01);
    guess_pt.set(1,3.843115617405354012e+00);
    guess_pt.set(2,-2.519775847962862159e+01);
    guess_pt.set(3,3.355540163929887409e+00);
    guess_pt.set(4,-1.922908387893200732e+00);
    guess_pt.set(5,2.537194900407488962e+00);
    guess_pt.set(6,-2.516959009732638197e+00);
    guess_pt.set(7,1.305366119640559397e+01);
    guess_pt.set(8,-1.022685240009726604e+00);
    guess_pt.set(9,9.789352071177765069e+00);
    guess_pt.set(10,-9.509963639307684957e+00);
    guess_pt.set(11,1.741218781414513117e+01);
    guess_pt.set(12,-2.913483504016468473e+01);
    guess_pt.set(13,-2.119019883927364489e+01);
    guess_pt.set(14,-6.021183229708010343e-01);
    guess_pt.set(15,5.574040528747243428e+00);
    guess_pt.set(16,6.070909881360219806e+00);
    guess_pt.set(17,2.289511182100715914e+00);
    guess_pt.set(18,-4.321251923659932714e+00);
    guess_pt.set(19,4.167919802338357349e-01);
    guess_pt.set(20,-7.665486974048813629e+00);
    guess_pt.set(21,-4.703327291811577382e+00);
    guess_pt.set(22,1.311021940994296564e+01);
    guess_pt.set(23,6.575075801015756838e+00);

    printf("chi guess %e\n",chisq[0](guess_pt));
    mcmc_test.guess_bases(guess_pt,chiLimit,1);*/
    mcmc_test.guess_bases(chiLimit,1);
}
mcmc_test.set_burnin(burnin);
mcmc_test.set_name_root(nameroot);
mcmc_test.sample(nsamples);
mcmc_test.write_timing(0);

FILE *log_test;
if(logname[0]!=0){
    log_test = fopen(logname, "w");
    fprintf(log_test,"# ");
    for(i=0;i<dim;i++){
        fprintf(log_test,"p%d ",i);
    }
    fprintf(log_test,"chisq mu sig ling\n");
    for(i=0;i<chisq->pt_log.get_rows();i++){
        for(j=0;j<dim;j++){
            fprintf(log_test,"%.12e ",chisq->pt_log.get_data(i,j));
        }
        fprintf(log_test,"%.12e 0 0 0\n",chisq->fn_log.get_data(i));
    }
}

}
