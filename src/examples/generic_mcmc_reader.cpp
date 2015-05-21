#include <stdlib.h>
#include "mcmc/chain.h"

int main(int iargc, char *argv[]){

char inNameRoot[letters],outNameRoot[letters];
int nChains,dim;

double chimin,chimax,dchi;

chimin=0.0;
chimax=20.0;
dchi=0.3;

int burnin,limit;
double confidenceLimit=0.95;

burnin=1000;
limit=-1;

dim=22;
nChains=4;

sprintf(inNameRoot,"chains/jellyBean_d22_chain_0");
sprintf(outNameRoot,"processedChains/jellyBean_d22_all");

array_1d<int> xdexes,ydexes;
xdexes.set_name("xdexes");
ydexes.set_name("ydexes");

int i,j;
for(i=1;i<iargc;i++){
    if(argv[i][0]=='-'){
        switch(argv[i][1]){
            case 'h':
                printf("b = burnin\ni = inNameRoot\no = outNameRoot\n");
                printf("d = dim\nc = nChains\nx = xdex, ydex\n");
                printf("l = limit\np = confidence limit\n");
                exit(1);
            case 'b':
                i++;
                burnin=atoi(argv[i]);
                break;
            case 'i':
                i++;
                for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                    inNameRoot[j]=argv[i][j];
                }
                inNameRoot[j]=0;
                break;
            case 'o':
                i++;
                for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                    outNameRoot[j]=argv[i][j];
                }
                outNameRoot[j]=0;
                break;
            case 'd':
                i++;
                dim=atoi(argv[i]);
                break;
            case 'c':
                i++;
                nChains=atoi(argv[i]);
                break;
            case 'l':
                i++;
                limit=atoi(argv[i]);
                break;
            case 'x':
                i++;
                xdexes.add(atoi(argv[i]));
                i++;
                ydexes.add(atoi(argv[i]));
                break; 
            case 'p':
                i++;
                confidenceLimit=atof(argv[i]);
                break;   
        
        }
    }
}

arrayOfChains chains(nChains, dim, NULL);

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

for(i=0;i<dim;i++){
    printf("%d %e %e %e\n",i,R.get_data(i),V.get_data(i),W.get_data(i));
}

int ix,iy;
if(xdexes.get_dim()==0){
    for(ix=0;ix<dim;ix++){
        for(iy=ix+1;iy<dim;iy++){
            xdexes.add(ix);
            ydexes.add(iy);
        }
    }
}

for(i=0;i<xdexes.get_dim();i++){
    ix=xdexes.get_data(i);
    iy=ydexes.get_data(i);
    chains.plot_contours(ix,iy,confidenceLimit,outNameRoot);
    printf("plotted %d %d\n",ix,iy);

}

int trueLimit;
if(limit<0)trueLimit=limit;
else trueLimit=burnin+limit;

chains.plot_chisquared_histogram(trueLimit,chimin,chimax,dchi,outNameRoot);


}
