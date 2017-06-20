#include "mcmc.h"
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


int main(int iargc, char *argv[]){

nonGaussianLump12 chifn;
mcmc_sampler sampler;
sampler.set_seed(124);
sampler.set_chisq_fn(&chifn);

char preburner_name[letters];
sprintf(preburner_name, "test_ellipse_bases.txt");

int dim=12;

array_2d<double> bases;
array_1d<double> radii,center;
array_2d<double> points;

bases.set_dim(dim,dim);
radii.set_dim(dim);

Ran dice(47);

FILE *input;
input=fopen(preburner_name,"r");
double xx;
int i;
for(i=0;i<dim;i++){
    fscanf(input,"%le",&xx);
    center.set(i,xx);
}
int j;
for(i=0;i<dim;i++){
    fscanf(input,"%le",&xx);
    radii.set(i,xx);
}
for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
        fscanf(input,"%le",&xx);
        bases.set(i,j,xx);
    }
}

sampler.set_bases(bases,radii);

/*select_start_pts(12,21.03,24,
                "output/workspace/jellyBean_d12_s90_output.sav",
                points);*/

select_start_pts(12,21.03,24,
                "output/workspace/lump_d12_s13_output.sav",
                points);


/*for(i=0;i<dim;i++){
    points.add_row(center);
}*/

sampler.set_start_pts(points);
sampler.set_name_root("chains/jb_170619/test_chain");

printf("sampling\n");
sampler.sample(100000);
sampler.write_chains();

}
