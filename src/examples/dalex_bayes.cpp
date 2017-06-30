#include "containers.h"
#include "goto_tools.h"
#include "kd.h"

int main(int iargc, char *argv[]){

    int i,j,dim;
    char in_name[letters];
    char out_name[letters];
    array_1d<int> xdexes,ydexes;
    xdexes.set_name("xdexes");
    ydexes.set_name("ydexes");
    in_name[0]=0;
    out_name[0]=0;
    dim=-1;
    for(i=1;i<iargc;i++){
        if(argv[i][0]=='-'){
            switch(argv[i][1]){
                case 'i':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        in_name[j]=argv[i][j];
                    }
                    in_name[j]=0;
                    break;
                case 'o':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        out_name[j]=argv[i][j];
                    }
                    out_name[j]=0;
                    break;
                case 'd':
                    i++;
                    dim=atoi(argv[i]);
                    break;
                case 'x':
                    i++;
                    xdexes.add(atoi(argv[i]));
                    i++;
                    ydexes.add(atoi(argv[i]));
                    break;
            }
        }
    }

    if(dim<0){
        printf("need to specify dim\n");
        exit(1);
    }
    if(in_name[0]==0){
        printf("need to specify in_name\n");
        exit(1);
    }
    if(out_name[0]==0){
        printf("need to specify out_name\n");
        exit(1);
    }
    if(xdexes.get_dim()==0){
        printf("need to specify dimensions to plot\n");
        exit(1);
    }
    if(xdexes.get_dim()!=ydexes.get_dim()){
        printf("somehow got %d xes but %d ys\n",
        xdexes.get_dim(),ydexes.get_dim());
        exit(1);
    }

    int n_cols=0;
    char word[letters];
    word[0]=0;
    FILE *in_file;
    in_file = fopen(in_name, "r");
    while(compare_char("log", word)==0){
        fscanf(in_file,"%s",word);
        if(compare_char("#",word)==0){
            n_cols++;
        }
    }

    array_2d<double> dalex_pts;
    array_1d<double> dalex_chisq;
    array_1d<double> pt;
    double xx;
    dalex_pts.set_name("dalex_pts");
    dalex_chisq.set_name("dalex_chisq");
    pt.set_name("pt");

    printf("n_cols %d\n",n_cols);
    while(fscanf(in_file,"%le",&xx)>0){
        pt.set(0,xx);
        for(i=1;i<dim;i++){
            fscanf(in_file,"%le",&xx);
            pt.set(i,xx);
        }
        dalex_pts.add_row(pt);
        fscanf(in_file,"%le",&xx);
        dalex_chisq.add(xx);
        for(i=dim+1;i<n_cols;i++){
            fscanf(in_file,"%le",&xx);
        }
    }

    fclose(in_file);

    printf("dalex_pts %d %d; dalex_chisq %d\n",
    dalex_pts.get_rows(),dalex_pts.get_cols(),dalex_chisq.get_dim());

    if(dalex_pts.get_rows()!=dalex_chisq.get_dim()){
        printf("somehow got %d pts but %d chisq\n",
        dalex_pts.get_rows(),dalex_chisq.get_dim());

        exit(1);
    }

}
