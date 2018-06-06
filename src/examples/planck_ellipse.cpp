#include "ellipse.h"

int main(int iargc, char *argv[]){

    char planck_name[500];
    sprintf(planck_name,"output/planck/planck_out_high_q2.txt");
    FILE *in_file;
    in_file=fopen(planck_name,"r");
    Ran dice(8812);
    int dim=33;
    char word[letters];
    int i;
    for(i=0;i<dim+3;i++){
        fscanf(in_file,"%s",word);
    }
    double mu;
    array_2d<double> data;
    data.set_name("data");
    double chisq_min=2903.33;
    double delta_chisq=47.41;
    array_1d<double> row;
    row.set_name("row");
    double roll;
    while(fscanf(in_file,"%le",&mu)>0){
        row.set(0,mu);
        for(i=1;i<dim;i++){
            fscanf(in_file,"%le",&mu);
            row.set(i,mu);
        }
        fscanf(in_file,"%le",&mu);
        if(mu<chisq_min+delta_chisq){
            roll=dice.doub();
            if(roll<0.2){
                data.add_row(row);
            }
        }
        fscanf(in_file,"%le",&mu);
    }
    fclose(in_file);

    ellipse planck_ellipse;
    printf("building ellipse\n");
    planck_ellipse.build(data);

    FILE *out_file;
    out_file = fopen("output/planck/ellipse_pts.txt","w");
    double sgn;
    int j;
    for(i=0;i<dim;i++){
        for(sgn=-1.0;sgn<1.1;sgn+=1.0){
            for(j=0;j<dim;j++){
                mu=planck_ellipse.center(j)+sgn*planck_ellipse.radii(i)*planck_ellipse.bases(i,j);
                fprintf(out_file,"%e ",mu);
            }
            fprintf(out_file,"\n");
        }
    }
    fclose(out_file);

    out_file = fopen("output/planck/ellipse_bases.txt","w");
    for(i=0;i<dim;i++){
        fprintf(out_file,"%e ",planck_ellipse.radii(i));
    }
    fprintf(out_file,"\n");
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            fprintf(out_file,"%e ",planck_ellipse.bases(i,j));
        }
        fprintf(out_file,"\n");
    }
    fclose(out_file);

}
