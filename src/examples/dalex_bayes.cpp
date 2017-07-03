#include "containers.h"
#include "goto_tools.h"
#include "kd.h"
#include "exampleLikelihoods.h"

void select_start_pts(int dim, double dchi, int n_chains,
                      char *in_file, array_2d<double> &pts_out){

    int i,j;

    FILE *input;
    input = fopen(in_file,"r");
    int n_cols=0;
    char word[100];

    word[0]=0;
    while(compare_char("log",word)==0){
        fscanf(input,"%s", word);
        if(compare_char("#",word)==0){
            n_cols++;
        }
    }
    printf("n_cols %d\n",n_cols);

    array_2d<double> good_pts;
    good_pts.set_name("good_pts");
    array_1d<double> min_pt;
    min_pt.set_name("min_pt");
    array_1d<double> one_pt;
    one_pt.set_name("one_pt");
    double chi_min;
    double xx;
    int n_rows = 0;

    chi_min=exception_value;
    good_pts.set_cols(dim);
    while(fscanf(input,"%le", &xx)>0){
        one_pt.set(0,xx);
        for(i=1;i<dim;i++){
            fscanf(input,"%le",&xx);
            one_pt.set(i,xx);

        }
        fscanf(input,"%le",&xx);
        if(xx<chi_min){
            chi_min=xx;
            for(j=0;j<dim;j++){
                min_pt.set(j,one_pt.get_data(j));
            }
        }
        for(i=0;i<n_cols-dim-1;i++){
            fscanf(input,"%le",&xx);
        }
        n_rows++;
    }
    fclose(input);
    printf("chi_min %e n_rows %d\n",chi_min,n_rows);

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
            good_pts.add_row(one_pt);
        }
        for(j=0;j<n_cols-dim-1;j++){
            fscanf(input,"%le",&xx);
        }
    }
    fclose(input);

    pts_out.reset_preserving_room();

    /*for(i=0;i<dim/2;i++){
        pts_out.add_row(min_pt);
    }
    return;*/

    int k;
    double dd;
    double dd_min;
    double dd_min_max;
    array_1d<int> chosen;
    chosen.set_name("chosen");
    int to_use;
    while(pts_out.get_rows()<n_chains){
        dd_min_max=-2.0*exception_value;
        for(i=0;i<good_pts.get_rows();i++){
            if(chosen.contains(i)==1){
                continue;
            }
            dd_min=0.0;
            for(j=0;j<dim;j++){
                dd_min+=power(min_pt.get_data(j)-good_pts.get_data(i,j),2);
            }
            for(j=0;j<pts_out.get_rows();j++){
                dd=0.0;
                for(k=0;k<dim;k++){
                    dd+=power(good_pts.get_data(i,k)-pts_out.get_data(j,k),2);
                }
                if(dd<dd_min){
                    dd_min=dd;
                }
            }
            if(dd_min>dd_min_max){
                dd_min_max=dd_min;
                to_use=i;
            }

        }
        chosen.add(to_use);
        pts_out.add_row(good_pts(to_use));
    }

}


void pixellate(const array_1d<double> &pt,
               const array_1d<double> &dx,
               const array_1d<double> &min,
               array_1d<int> &px){

    int i,j;
    for(i=0;i<pt.get_dim();i++){
        for(j=0;min.get_data(i)+j*dx.get_data(i)<pt.get_data(i);j++);
        if(min.get_data(i)+j*dx.get_data(i)-pt.get_data(i)>0.5*dx.get_data(i)){
            j--;
        }
        if(fabs(min.get_data(i)+j*dx.get_data(i)-pt.get_data(i))>0.5*dx.get_data(i)){
            printf("pixellate failed\n");
            exit(1);
        }
        px.set(i,j);
    }

}

int main(int iargc, char *argv[]){

    double pixel_factor=0.05;
    int i,j,k,dim;
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

    word[0]=0;
    in_file = fopen(in_name, "r");
    while(compare_char("log", word)==0){
        fscanf(in_file,"%s",word);
    }
    printf("last word %s\n",word);

    int n_pts=0;
    while(fscanf(in_file,"%le",&xx)>0){
        n_pts++;
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
    printf("final xx %e\n",xx);

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

    printf("set min max %d %d %e\n",dim,xmax.get_dim(),chisq_min);
    printf("good pts %d out of %d %e\n",good_pts.get_rows(),n_pts,delta_chisq);

    for(i=0;i<dim;i++){
        dx.set(i,pixel_factor*(xmax.get_data(i)-xmin.get_data(i)));
    }

    int ix, iy;

    array_2d<int> pixel_list;
    pixel_list.set_name("pixel_list");
    array_1d<int> pixel;
    pixel.set_name("pixel");
    int is_valid,is_same;
    for(i=0;i<good_pts.get_rows();i++){
        pixellate(good_pts(i),dx,xmin,pixel);
        is_valid=1;
        for(j=0;j<pixel_list.get_rows();j++){
            is_same=1;
            for(k=0;k<dim;k++){
                if(pixel_list.get_data(j,k)!=pixel.get_data(k)){
                    is_same=0;
                    break;
                }
            }
            if(is_same==1){
                is_valid=0;
                break;
            }
        }
        if(is_valid==1){
            pixel_list.add_row(pixel);
        }
    }

    printf("n pixels %d\n",pixel_list.get_rows());
    FILE *out_file;
    out_file=fopen("pixel_test.txt","w");
    for(i=0;i<pixel_list.get_rows();i++){
        for(j=0;j<dim;j++){
            fprintf(out_file,"%le ",xmin.get_data(j)+pixel_list.get_data(i,j)*dx.get_data(j));
        }
        fprintf(out_file,"\n");
    }
    fclose(out_file);
    exit(1);

    gaussianJellyBean12 chifn;

    int n_accepted=0;
    int n_rejected=0;
    int accept_it;
    array_1d<double> chi_old;
    array_2d<double> sample_pts;
    array_1d<int> degen;
    degen.set_name("degen");
    chi_old.set_name("chi_old");
    sample_pts.set_name("sample_pts");
    double chi_new;
    Ran dice(88);
    int n_chains=2*dim;

    int n_samples=100000;
    int i_sample;
    int i_chain;
    pt.reset_preserving_room();
    int i_pixel;
    out_file=fopen(out_name,"w");
    fclose(out_file);

    /*select_start_pts(12,delta_chisq,n_chains,in_name,sample_pts);*/


    double sample_chi_min=2.0*exception_value;

    for(i_chain=0;i_chain<n_chains;i_chain++){
        i_pixel=dice.int32()%pixel_list.get_rows();
        for(i=0;i<dim;i++){
            pt.set(i,xmin.get_data(i)+pixel_list.get_data(i_pixel,i)*dx.get_data(i));
        }
        for(i=0;i<dim;i++){
            pt.add_val(i,dx.get_data(i)*(dice.doub()-0.5));
        }
        sample_pts.add_row(pt);
    }

    for(i=0;i<n_chains;i++){
        xx=chifn(sample_pts(i));
        if(xx<sample_chi_min){
            sample_chi_min=xx;
        }
        chi_old.set(i,xx);
        degen.set(i,1);
    }


    for(i_sample=0;i_sample<n_samples;i_sample++){
        for(i_chain=0;i_chain<n_chains;i_chain++){
            i_pixel=dice.int32()%pixel_list.get_rows();
            for(i=0;i<dim;i++){
                pt.set(i,xmin.get_data(i)+pixel_list.get_data(i_pixel,i)*dx.get_data(i));
            }
            for(i=0;i<dim;i++){
                pt.add_val(i,dx.get_data(i)*(dice.doub()-0.5));
            }
            chi_new=chifn(pt);
            if(chi_new<sample_chi_min){
                sample_chi_min=chi_new;
            }
            accept_it=0;
            if(chi_new<chi_old.get_data(i_chain)){
                accept_it=1;
            }
            else{
                if(dice.doub()<exp(-0.5*(chi_new-chi_old.get_data(i_chain)))){
                    accept_it=1;
                }
            }
            if(accept_it==0){
                n_rejected++;
                degen.add_val(i_chain,1);
            }
            else{
                n_accepted++;
                out_file=fopen(out_name,"a");
                fprintf(out_file,"%d %e ",degen.get_data(i_chain),chi_old.get_data(i_chain));
                for(i=0;i<dim;i++){
                    fprintf(out_file,"%e ",sample_pts.get_data(i_chain,i));
                }
                fprintf(out_file,"\n");
                fclose(out_file);
                degen.set(i_chain,1);
                chi_old.set(i_chain,chi_new);
                for(i=0;i<dim;i++){
                    sample_pts.set(i_chain,i,pt.get_data(i));
                }
            }
        }
        if(i_sample%10000==0){
            printf("n_accept %.3e n_reject %.3e chimin %e\n",
            float(n_accepted),float(n_rejected),sample_chi_min);
        }
    }

    for(i_chain=0;i_chain<n_chains;i_chain++){
        out_file=fopen(out_name,"a");
        fprintf(out_file,"%d %e ",degen.get_data(i_chain),chi_old.get_data(i_chain));
        for(i=0;i<dim;i++){
            fprintf(out_file,"%e ",sample_pts.get_data(i_chain,i));
        }
        fprintf(out_file,"\n");
        fclose(out_file);
    }

}
