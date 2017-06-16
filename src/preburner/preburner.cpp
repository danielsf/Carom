#include "ellipse.h"

int main(int iargc, char *argv[]){

    int i,j;
    char in_file[letters];
    char out_file[letters];
    double dchi;
    int dim;

    dim=-1;
    dchi=-1.0;
    in_file[0]=0;
    out_file[0]=0;

    for(i=1;i<iargc;i++){
        if(argv[i][0] == '-'){
            switch(argv[i][1]){
                case 'i':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        in_file[j]=argv[i][j];
                    }
                    in_file[j]=0;
                    break;
                case 'o':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        out_file[j]=argv[i][j];
                    }
                    out_file[j]=0;
                    break;
                case 'd':
                    i++;
                    dim=atoi(argv[i]);
                    break;
                case 'c':
                    i++;
                    dchi=atof(argv[i]);
                    break;
            }
        }
    }

    if(in_file[0]==0){
        printf("must specify in_file\n");
        exit(1);
    }

    if(out_file[0]==0){
        printf("must specify out_file\n");
        exit(1);
    }

    if(dim<0){
        printf("must specify dim\n");
        exit(1);
    }

    if(dchi<0.0){
        printf("must specify dchi\n");
        exit(1);
    }

    FILE *input;
    char word[letters];
    int n_cols=0;
    input = fopen(in_file, "r");
    word[0]=0;
    while(compare_char("log",word)==0){
        fscanf(input,"%s", word);
        if(compare_char("#",word)==0){
            n_cols++;
        }
    }
    printf("n_cols %d\n",n_cols);

    array_2d<double> pts;
    pts.set_name("preburner_pts");
    double chi_min;
    double xx;
    int n_rows = 0;

    chi_min=exception_value;
    pts.set_cols(dim);
    while(fscanf(input,"%le", &xx)>0){
        for(i=1;i<dim;i++){
            fscanf(input,"%le",&xx);
        }
        fscanf(input,"%le",&xx);
        if(xx<chi_min){
            chi_min=xx;
        }
        for(i=0;i<n_cols-dim-1;i++){
            fscanf(input,"%le",&xx);
        }
        n_rows++;
    }
    fclose(input);
    printf("chi_min %e n_rows %d\n",chi_min,n_rows);
    array_1d<double> one_pt;
    one_pt.set_name("preburner_one_pt");

    input = fopen(in_file, "r");
    for(i=0;i<n_cols+1;i++){
        fscanf(input,"%s",word);
    }
    for(i=0;i<n_rows;i++){
        for(j=0;j<dim;j++){
            fscanf(input,"%le", &xx);
            one_pt.set(j,xx);
        }
        fscanf(input,"%le",&xx);
        if(xx<chi_min+dchi){
            pts.add_row(one_pt);
        }
        for(j=0;j<n_cols-dim-1;j++){
            fscanf(input,"%le",&xx);
        }
    }
    fclose(input);

    ellipse good_ellipse;
    good_ellipse.use_geo_center();
    good_ellipse.build(pts);

    FILE *output;
    output=fopen(out_file, "w");
    printf("writing to %s\n",out_file);
    for(i=0;i<dim;i++){
        fprintf(output,"%e ",good_ellipse.center(i));
    }
    fprintf(output,"\n");
    for(i=0;i<dim;i++){
        fprintf(output,"%e ",good_ellipse.radii(i));
    }
    fprintf(output,"\n");
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            fprintf(output,"%e ",good_ellipse.bases(i,j));
        }
        fprintf(output,"\n");
    }
}
