#include "dalex_driver.h"
#include "planck.h"

int main(int iargc, char *argv[]){

    int planck_dim = 33;
    int seed=11923;
    int i,j;
    char timing_name[letters];
    char out_name[letters];
    char prior_name[letters];
    char warm_name[letters];
    out_name[0] = 0;
    timing_name[0] = 0;
    prior_name[0] = 0;
    warm_name[0] = 0;

    char log_name[letters];
    log_name[0]=0;

    for(i=1;i<iargc;i++){
        if(argv[i][0] == '-'){
            switch(argv[i][1]){
                case 'p':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        prior_name[j] = argv[i][j];
                    }
                    prior_name[j] = 0;
                    break;
                case 't':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        timing_name[j] = argv[i][j];
                    }
                    timing_name[j] = 0;
                    break;
                case 'w':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        warm_name[j] = argv[i][j];
                    }
                    warm_name[j] = 0;
                    break;
                case 'o':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        out_name[j]=argv[i][j];
                    }
                    out_name[j]=0;
                    break;
                case 'l':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        log_name[j]=argv[i][j];
                    }
                    log_name[j]=0;
                    break;
                case 's':
                    i++;
                    seed = atoi(argv[i]);
                    break;
            }
        }
    }

    if(out_name[0] == 0){
        printf("out name not specified\n");
        exit(1);
    }
    if(timing_name[0] == 0){
         printf("timing name not specified\n");
         exit(1);
    }
    if(prior_name[0] == 0){
        printf("prior name not specified\n");
        exit(1);
    }

    array_1d<double> min,max;
    double mu1,mu2;

    FILE *prior_file;
    prior_file = fopen(prior_name, "r");    
    for(i=0;i<planck_dim;i++){
        fscanf(prior_file,"%le %le",&mu1,&mu2);
        min.set(i,mu1);
        max.set(i,mu2);
    }
    fclose(prior_file);

    DalexPlanckLikelihood chifn;
    for(i=0;i<planck_dim;i++){
        chifn.set_min(i,min.get_data(i));
        chifn.set_max(i,max.get_data(i));
    }    

    printf("chifn dim %d\n",chifn.get_dim());

    FILE *log_file;
    if(log_name[0]!=0){
        log_file=fopen(log_name,"w");
        fprintf(log_file,"timing: %s\n",timing_name);
        fprintf(log_file,"output: %s\n",out_name);
        fclose(log_file);
    }

    dalex_driver dalex_test;
    dalex_test.set_log_file_name(log_name);
    dalex_test.set_deltachi(47.41); // set delta chi^2 defining chi^2_lim
    dalex_test.set_seed(seed); // seed the random number generator

    dalex_test.set_min(min); // set the minimum range in parameter space
    dalex_test.set_max(max); // set the maximum range in parameter space
    // the above limits just bound the initial random samples and the
    // seeds of the first optimization search; if Dalex finds interesting
    // regions of parameter space outside of these bounds, it will go there

    dalex_test.set_chisquared(&chifn); // pass in the chi^2 function
    dalex_test.set_timingname(timing_name);
    dalex_test.set_outname(out_name);
    dalex_test.initialize(2*planck_dim); // randomly sample 2D points
    dalex_test.set_write_every(10000); // how often to write output
    if(warm_name[0]==0){
        dalex_test.search(10000000); // sample 1,700,000 points
    }
    else{
        dalex_test.warm_start(warm_name,10000000);
    }
}
