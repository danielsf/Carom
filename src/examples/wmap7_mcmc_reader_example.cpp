#include "mcmc/chain.h"

int main(){

char inNameRoot[letters],outNameRoot[letters];
int nChains,dim;

dim=6;
nChains=4;

sprintf(inNameRoot,"chains/mcmc_test_150424");
sprintf(outNameRoot,"processedChains/mcmc_test_150424");

arrayOfChains chains(nChains, dim, NULL);

int i;
char inName[letters];
for(i=0;i<nChains;i++){
    sprintf(inName,"%s_%d.txt",inNameRoot,i);
    chains(i)->read_chain(inName);
}

array_1d<double> V,W,R;
R.set_name("R");
V.set_name("V");
W.set_name("W");

chains.get_independent_samples(0.1,-1);
chains.calculate_R(R,V,W);
int ix,iy;
for(ix=0;ix<6;ix++){
    if(ix!=3){
        for(iy=i+1;iy<6;iy++){
            if(iy!=3){
                chains.plot_contours(ix,iy,0.95,outNameRoot);
            }
        }
    }
}


}
