#include "controls/control_integrator.h"
#include "jellyBean.h"

int main(){

    double curvature_radius=20.0;
    double radial_sigma=0.02;
    double angular_width=0.4;

    int dim=4;
    jellyBeanData chisq(dim,1,100,0.4,angular_width,radial_sigma,curvature_radius);
    //ellipseData chisq(dim,1,100,0.4);



    array_1d<double> min,max,dx,raw_min,raw_max;

    //for jellyBean
    raw_min.set(0,-7.762234e-01);
    raw_max.set(0,8.082285e+00);
    raw_min.set(1,-9.808949e+00);
    raw_max.set(1, -3.364088e+00);
    raw_min.set(2,-1.413286e-01);
    raw_max.set(2,5.594381e+00);
    raw_min.set(3,3.089276e+00);
    raw_max.set(3,1.252334e+01);

    //for ellipse
    /*raw_min.set(0,7.283857e+00);
    raw_max.set(0,7.369758e+00);
    raw_min.set(1, -2.856360e+01);
    raw_max.set(1,-2.854585e+01);
    raw_min.set(2,6.163477e+00);
    raw_max.set(2,6.202135e+00);
    raw_min.set(3,2.768558e+01);
    raw_max.set(3,2.774378e+01);*/

    //for ellipse nonGaussian
    /*raw_min.set(0,7.275115e+00);
    raw_max.set(0,7.369179e+00);
    raw_min.set(1,-2.857240e+01);
    raw_max.set(1,-2.854504e+01);
    raw_min.set(2,6.162698e+00);
    raw_max.set(2,6.213527e+00);
    raw_min.set(3,2.768490e+01);
    raw_max.set(3,2.775592e+01);*/

    int i;
    double dd;
    for(i=0;i<4;i++){
        dd=raw_max.get_data(i)-raw_min.get_data(i);
        min.set(i,raw_min.get_data(i)-2.0*dd);
        max.set(i,raw_max.get_data(i)+2.0*dd);
    }

    for(i=0;i<4;i++){
        dx.set(i,(max.get_data(i)-min.get_data(i))*0.005);
        printf("dx %e\n",dx.get_data(i));
    }

    dx.multiply_val(0,0.3333);

    control_integrator integrator(chisq,min,max,dx,"controls/jellyBeanData/jellyBean");
    integrator.run_analysis();

    printf("\nrange\n");
    for(i=0;i<4;i++){
        printf("%e %e\n",min.get_data(i),max.get_data(i));
    }

    printf("that took %e\n",double(time(NULL))-start);
}
