#include "mcmc/chain.h"

int main(){

char inNameRoot[letters],outNameRoot[letters];
int nChains,dim;

dim=6;
nChains=4;

//sprintf(inNameRoot,"chains/mcmc_test_150424");
//sprintf(outNameRoot,"processedChains/mcmc_test_150424");

sprintf(inNameRoot,"/Users/danielsf/physics/recreate_getdist/ieuchains_1304/wmap7_reformatted");
sprintf(outNameRoot,"processedChains/ieu15k");

arrayOfChains chains(nChains, dim, NULL);

int i;
char inName[letters];
for(i=0;i<nChains;i++){
    sprintf(inName,"%s_%d.txt",inNameRoot,i);
    chains(i)->read_chain(inName);
}

printf("acceptance rate %e\n",chains.acceptance_rate());

array_1d<double> V,W,R;
R.set_name("R");
V.set_name("V");
W.set_name("W");

chains.get_independent_samples(0.1,15000);
chains.calculate_R(R,V,W);
printf("got R,V,W\n");
int ix,iy;
for(ix=0;ix<dim;ix++){
    if(ix!=3){
        for(iy=ix+1;iy<dim;iy++){
            if(iy!=3){
                chains.plot_contours(ix,iy,0.95,outNameRoot);
            }
        }
    }
}

for(i=0;i<dim;i++){
    printf("%d %e %e %e\n",i,R.get_data(i),V.get_data(i),W.get_data(i));
}

int j;

FILE *output;
output=fopen("processedChains/independent_samples.txt","w");
for(i=0;i<chains.get_n_samples();i++){
    for(j=0;j<dim;j++){
        fprintf(output,"%e ",chains.get_sample(i,j));
    }
    fprintf(output,"\n");
}
fclose(output);

}
