#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "jellyBean.h"
#include "exampleLikelihoods.h"
#include "dalex_driver.h"

int main(int iargc, char *argv[]){
    gaussianJellyBean12 chifn;
    int i,j;

    char data_name[letters];
    sprintf(data_name,"output/test_180214/output.txt_end_pts.txt");

    int dim=12;

    array_2d<double> pts;
    array_1d<double> row;
    array_1d<double> fn,min,max;

    for(i=0;i<dim;i++){
        min.set(i,2.0*exception_value);
        max.set(i,-2.0*exception_value);
    }

    double mu;

    FILE *in_file;
    in_file=fopen(data_name, "r");
    while(fscanf(in_file,"%le",&mu)>0){
        row.set(0,mu);
        if(mu<min.get_data(0)){
            min.set(0,mu);
        }
        if(mu>max.get_data(0)){
            max.set(0,mu);
        }
        for(i=1;i<dim;i++){
            fscanf(in_file,"%le",&mu);
            row.set(i,mu);
            if(mu<min.get_data(i)){
                min.set(i,mu);
            }
            if(mu>max.get_data(i)){
                max.set(i,mu);
            }
        }
        pts.add_row(row);
        fscanf(in_file,"%le",&mu);
        fn.add(mu);
        fscanf(in_file,"%d",&j);
    }
    fclose(in_file);
}
