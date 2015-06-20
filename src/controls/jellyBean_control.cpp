#include "controls/control_integrator.h"
#include "jellyBean.h"

int main(){
    
    int dim=4;
    jellyBeanData chisq(dim,1,100,0.4,0.4,0.05,20.0);
    
    
    
    array_1d<double> min,max,dx;    
    min.set(0,-10.0);
    max.set(0,25.0);
    
    min.set(1,-13.0);
    max.set(1,10.0);
    
    min.set(2,-9.0);
    max.set(2,14.0);
    
    min.set(3,-5.0);
    max.set(3,25.0);
    
    int i;
    for(i=0;i<4;i++){
        dx.set(i,0.01*(max.get_data(i)-min.get_data(i)));
    }
    
    control_integrator integrator(chisq,min,max,dx,"output/scratch/control");
    integrator.run_analysis();

}
