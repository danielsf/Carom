#include <stdlib.h>
#include "goto_tools.h"

int main(int iargc, char *argv[]){

char inName[letters];
char outName[letters];

inName[0]=0;
outName[0]=0;

int dim=0;

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
            case 'h':
                printf("i = inName\no = outName\nd = dim\n");
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

FILE *input;

}
