#include <stdlib.h>
#include "kde.h"

int main(int iargc, char *argv[]){

array_1d<double> mins,maxes;
mins.set_name("mins");
maxes.set_name("maxes");

array_1d<int> xdexes,ydexes;
xdexes.set_name("xdexes");
ydexes.set_name("ydexes");
char inName[letters],outNameRoot[letters];
int dim,smoothby=3;
double confidenceLimit=0.95;
double pixelFactor=0.01;

int i,j;
for(i=1;i<iargc;i++){
    if(argv[i][0]=='-'){
        switch(argv[i][1]){
            case 'o':
                i++;
                transcribe(argv[i],outNameRoot);
                break;
            case 'i':
                i++;
                transcribe(argv[i],inName);
                break;
            case 'd':
                i++;
                dim=atoi(argv[i]);
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
            case 's':
                i++;
                smoothby=atoi(argv[i]);
                break;
            case 'f':
                i++;
                pixelFactor=atof(argv[i]);
                break;
            case 'h':
                printf("o = outname\ni = inName\nd = dim\nx = dexes\n");
                printf("p = confidence limit\ns = smoothby\n");
                printf("f = pixelFactor\n");
                exit(1);
                break;
        }
    }
}

FILE *input;
array_1d<double> wgt,vv;
array_2d<double> data;

double mu;

for(i=0;i<dim;i++){
    mins.set(i,2.0*exception_value);
    maxes.set(i,-2.0*exception_value);
}

wgt.set_name("wgt");
vv.set_name("vv");
data.set_name("data");
input=fopen(inName,"r");
while(fscanf(input,"%le",&mu)>0){
    wgt.add(mu);
    fscanf(input,"%le",&mu);
    for(i=0;i<dim;i++){
        fscanf(input,"%le",&mu);
        vv.set(i,mu);
        
        if(mu<mins.get_data(i))mins.set(i,mu);
        if(mu>maxes.get_data(i))maxes.set(i,mu);
    }
    data.add_row(vv);
}

kde density;
density.set_data(&data,wgt);

int ix,iy;
if(xdexes.get_dim()==0){
    for(ix=0;ix<dim;ix++){
        for(iy=ix+1;iy<dim;iy++){
            xdexes.add(ix);
            ydexes.add(iy);
        }
    }
}

double dx,dy;

char outname[2*letters];
for(i=0;i<xdexes.get_dim();i++){
    ix=xdexes.get_data(i);
    iy=ydexes.get_data(i);
    
    dx=pixelFactor*(maxes.get_data(ix)-mins.get_data(ix));
    dy=pixelFactor*(maxes.get_data(iy)-mins.get_data(iy));
    
    sprintf(outname,"%s_%d_%d_contour.txt",outNameRoot,ix,iy);
    density.plot_boundary(ix,dx,iy,dy,confidenceLimit,outname,3);
    printf("plotted %d %e %d %e\n",ix,dx,iy,dy);
}

}
