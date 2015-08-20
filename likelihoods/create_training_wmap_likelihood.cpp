#include "kd.h"

int main(int iargc, char *argv[]){

    char in_name[letters];
    in_name[0]=0;

    char out_name[letters];
    out_name[0]=0;

    char tree_name[letters];
    tree_name[0]=0;

    int i,j;
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

                case 't':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        tree_name[j]=argv[i][j];
                    }
                    tree_name[j]=0;
                break;

                case 'o':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        out_name[j]=argv[i][j];
                    }
                    out_name[j]=0;
                break;

            }
        }
    }

    if(in_name[0]==0){
        printf("WARNING did not set in_name\n");
        exit(1);
    }

    if(out_name[0]==0){
        printf("WARNING did not set out_name\n");
        exit(1);
    }

    if(tree_name[0]==0){
        printf("WARNING did not set tree_name\n");
        exit(1);
    }



}
