#include "mcmc/chain.h"

int main(){

char inNameRoot[letters],outNameRoot[letters];
int nChains,dim;

double chimin,chimax,dchi;

chimin=0.0;
chimax=20.0;
dchi=0.3;

int burnin,limit;

burnin=1000;
limit=-1;

dim=5;
nChains=4;

sprintf(inNameRoot,"chains/jellyBean_chain_0");
sprintf(outNameRoot,"processedChains/jellyBean_all");

arrayOfChains chains(nChains, dim, NULL);

int i;
char inName[letters];
for(i=0;i<nChains;i++){
    sprintf(inName,"%s_%d.txt",inNameRoot,i);
    chains(i)->read_chain(inName);
}

printf("acceptance rate %e\n",chains.acceptance_rate());
printf("points %d\n",chains.get_points());

array_1d<double> V,W,R;
R.set_name("R");
V.set_name("V");
W.set_name("W");

//chains.get_independent_samples(0.1,burnin,limit);
chains.use_all(burnin,limit);
chains.acceptance_statistics(burnin,limit);

chains.calculate_R(R,V,W);
printf("got R,V,W\n");
int ix,iy;
for(ix=0;ix<dim;ix++){
    for(iy=ix+1;iy<dim;iy++){
            chains.plot_contours(ix,iy,0.95,outNameRoot);
            printf("plotted %d %d\n",ix,iy);
        }
}

chains.plot_chisquared_histogram(burnin+limit,chimin,chimax,dchi,outNameRoot);

for(i=0;i<dim;i++){
    printf("%d %e %e %e\n",i,R.get_data(i),V.get_data(i),W.get_data(i));
}


}
