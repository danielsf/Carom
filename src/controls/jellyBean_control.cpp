#include "controls/control_integrator.h"
#include "jellyBean.h"

int main(){
    
    double curvature_radius=20.0;
    double radial_sigma=0.05;
    double angular_width=0.4;
    
    int dim=4;
    jellyBeanData chisq(dim,1,100,0.4,angular_width,radial_sigma,curvature_radius);
    
    
    
    array_1d<double> min,max,dx;    
    min.set(0,-10.0);
    max.set(0,20.0);
    
    min.set(1,-30.0);
    max.set(1,10.0);
    
    min.set(2,-10.0);
    max.set(2,20.0);
    
    min.set(3,-23.0);
    max.set(3,40.0);
    
    int i;
    for(i=0;i<4;i++){
        dx.set(i,0.01*(max.get_data(i)-min.get_data(i)));
    }
    
    control_integrator integrator(chisq,min,max,dx,"output/scratch/control");
    integrator.run_analysis();

}
