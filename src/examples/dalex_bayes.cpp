#include "containers.h"
#include "goto_tools.h"
#include "kd.h"

int main(int iargc, char *argv[]){

    double pixel_factor=0.1;
    int i,j,dim;
    char in_name[letters];
    char out_name[letters];
    double delta_chisq=-1.0;
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
                case 'c':
                    i++;
                    delta_chisq=atof(argv[i]);
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
    if(delta_chisq<0.0){
        printf("need to specify delta_chisq\n");
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

    array_1d<double> dalex_chisq;
    double xx;
    dalex_chisq.set_name("dalex_chisq");

    printf("n_cols %d\n",n_cols);
    while(fscanf(in_file,"%le",&xx)>0){
        for(i=1;i<dim;i++){
            fscanf(in_file,"%le",&xx);
        }
        fscanf(in_file,"%le",&xx);
        dalex_chisq.add(xx);
        for(i=dim+1;i<n_cols;i++){
            fscanf(in_file,"%le",&xx);
        }
    }

    fclose(in_file);

    double chisq_min=exception_value;
    for(i=0;i<dalex_chisq.get_dim();i++){
        if(dalex_chisq.get_data(i)<chisq_min){
            chisq_min=dalex_chisq.get_data(i);
        }
    }

    array_1d<double> pt;
    pt.set_name("pt");
    array_2d<double> good_pts;
    good_pts.set_name("good_pts");

    in_file = fopen(in_name, "r");
    while(compare_char("log", word)==0){
        fscanf(in_file,"%s",word);
        if(compare_char("#",word)==0){
            n_cols++;
        }
    }

    while(fscanf(in_file,"%le",&xx)>0){
        pt.set(0,xx);
        for(i=1;i<dim;i++){
            fscanf(in_file,"%le",&xx);
            pt.set(i,xx);
        }
        fscanf(in_file,"%le",&xx);
        if(xx<=chisq_min+delta_chisq){
            good_pts.add_row(pt);
        }
        for(i=dim+1;i<n_cols;i++){
            fscanf(in_file,"%le",&xx);
        }
    }
    fclose(in_file);

    array_1d<double> dx;
    dx.set_name("dx");
    array_1d<double> xmin,xmax;
    xmin.set_name("xmin");
    xmax.set_name("xmax");
    for(i=0;i<good_pts.get_rows();i++){
        for(j=0;j<dim;j++){
            if(j>=xmin.get_dim() || good_pts.get_data(i,j)<xmin.get_data(j)){
                xmin.set(j,good_pts.get_data(i,j));
            }
            if(j>=xmax.get_dim() || good_pts.get_data(i,j)>xmax.get_data(j)){
                xmax.set(j,good_pts.get_data(i,j));
            }
        }
    }

    for(i=0;i<dim;i++){
        dx.set(i,pixel_factor*(xmax.get_data(i)-xmin.get_data(i)));
    }

    kd_tree *forest;
    forest = new kd_tree[dim*dim];

    int ix, iy;
    array_2d<double> pts2d;
    pts2d.set_name("pts2d");
    array_1d<double> max2d,min2d;
    max2d.set_name("max2d");
    min2d.set_name("min2d");
    for(ix=0;ix<dim;ix++){
        for(iy=ix+1;iy<dim;iy++){
            pts2d.reset_preserving_room();
            max2d.reset_preserving_room();
            min2d.reset_preserving_room();
            pts2d.set_dim(good_pts.get_rows(),2);
            for(i=0;i<good_pts.get_rows();i++){
                pts2d.set(i,0,good_pts.get_data(i,ix));
                pts2d.set(i,1,good_pts.get_data(i,iy));
                if(i==0 || good_pts.get_data(i,ix)<min2d.get_data(0)){
                    min2d.set(0,good_pts.get_data(i,ix));
                }
                if(i==0 || good_pts.get_data(i,ix)>max2d.get_data(0)){
                    max2d.set(0,good_pts.get_data(i,ix));
                }
                if(i==0 || good_pts.get_data(i,iy)<min2d.get_data(1)){
                    min2d.set(1,good_pts.get_data(i,iy));
                }
                if(i==0 || good_pts.get_data(i,iy)>max2d.get_data(1)){
                    max2d.set(1,good_pts.get_data(i,iy));
                }
            }
            forest[ix*dim+iy].build_tree(pts2d,min2d,max2d);
        }
    }

    delete [] forest;
}
