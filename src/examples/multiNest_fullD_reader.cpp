#include <stdlib.h>
#include "goto_tools.h"

int main(int iargc, char *argv[]){

char inName[letters];
char outName[letters];

inName[0]=0;
outName[0]=0;

int dim=0;
double fraction=0.95;

int i,j;
for(i=1;i<iargc;i++){
    if(argv[i][0]=='-'){
        switch(argv[i][1]){
            case 'i':
                i++;
                transcribe(argv[i],inName);
                break;
            case 'o':
                i++;
                transcribe(argv[i],outName);
                break;
            case 'd':
                i++;
                dim=atoi(argv[i]);
                break;
            case 'p':
                i++;
                fraction=atof(argv[i]);
                break;
            case 'h':
                printf("i = inName\no = outName\nd = dim\n");
                printf("p = confidence limit\n");
                exit(1);
                break;
        }
    }
}

if(dim==0 || inName[0]==0 || outName[0]==0){
    printf("WARNING could not proceed dim %d inname %s outname %s\n",
    dim,inName,outName);
    exit(1);
}

array_1d<double> wgt,wgt_sorted;
array_2d<double> points;
array_1d<int> dexes;

wgt.set_name("wgt");
wgt_sorted.set_name("wgt_sorted");
points.set_name("points");
dexes.set_name("dexes");

double total,sum,mu;
int ct;

points.set_cols(dim);

FILE *input;
total=0.0;
input=fopen(inName,"r");
if(input==NULL){
    printf("\nWARNING could not open %s\n\n",inName);
    exit(1);
}
ct=0;
while(fscanf(input,"%le",&mu)>0){
    total+=mu;
    wgt.add(mu);
    fscanf(input,"%le",&mu);
    for(j=0;j<dim;j++){
        fscanf(input,"%le",&mu);
        points.set(ct,j,mu);
    }
    dexes.add(ct);
    ct++;
}
fclose(input);

sort(wgt,wgt_sorted,dexes);

int dex;
input=fopen(outName,"w");
sum=0.0;
for(i=dexes.get_dim()-1;i>=0 && sum<fraction*total;i--){
    dex=dexes.get_data(i);
    sum+=wgt_sorted.get_data(i);
    for(j=0;j<dim;j++){
        fprintf(input,"%e ",points.get_data(dex,j));
    }
    fprintf(input,"\n");
}
fclose(input);


}
