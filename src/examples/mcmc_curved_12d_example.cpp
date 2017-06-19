#include "mcmc.h"
#include "exampleLikelihoods.h"

int main(int iargc, char *argv[]){

gaussianJellyBean12 chifn;
mcmc_sampler sampler;
sampler.set_seed(124);
sampler.set_chisq_fn(&chifn);

char preburner_name[letters];
sprintf(preburner_name, "test_jelly_bases.txt");

int dim=12;

array_2d<double> bases;
array_1d<double> radii,center;
array_2d<double> points;

bases.set_dim(dim,dim);
radii.set_dim(dim);
points.set_dim(dim,dim);

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
    radii.set(i,xx*0.1);
}
for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
        fscanf(input,"%le",&xx);
        bases.set(i,j,xx);
    }
}

sampler.set_bases(bases,radii);
for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
        points.set(i,j,center.get_data(j));
    }
}

sampler.set_start_pts(points);
sampler.set_name_root("chains/jb_170619/test_chain");

sampler.sample(5000);
sampler.write_chains();

}
