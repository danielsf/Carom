#include "controls/control_integrator.h"
#include "jellyBean.h"
#include "exampleLikelihoods.h"

int main(){

    double start=double(time(NULL));
    double curvature_radius=20.0;
    double radial_sigma=0.02;
    double angular_width=0.4;

    int dim=4;
    //jellyBeanData chisq(dim,1,100,0.4,angular_width,radial_sigma,curvature_radius);
    //ellipseData chisq(dim,1,100,0.4);
    integrableJellyBean chisq;


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

    //for gaussianJellyBean
    raw_min.set(0,1.561333e+01);
    raw_max.set(0,2.680832e+01);
    raw_min.set(1,-2.049267e+01);
    raw_max.set(1, -9.829029e+00);
    raw_min.set(2,9.726976e+00);
    raw_max.set(2,1.610506e+01);
    raw_min.set(3,1.212955e+01);
    raw_max.set(3,2.189316e+01);

    //for integrableJellyBean
    raw_min.set(0,-3.0);
    raw_max.set(0,13.0);
    raw_min.set(1,-2.0);
    raw_max.set(1, 15.0);
    raw_min.set(2,-12.0);
    raw_max.set(2,-4.0);
    raw_min.set(3,-16.0);
    raw_max.set(3,-3.0);

    int i;
    for(i=0;i<4;i++){
        raw_min.subtract_val(i,7.0);
        raw_max.add_val(i,7.0);
    }
    double dd;
    for(i=0;i<4;i++){
        //dd=raw_max.get_data(i)-raw_min.get_data(i);
        min.set(i,raw_min.get_data(i));
        max.set(i,raw_max.get_data(i));
    }

    for(i=0;i<4;i++){
        dx.set(i,(max.get_data(i)-min.get_data(i))*0.005);
        printf("dx %e\n",dx.get_data(i));
    }
    //dx.multiply_val(0,0.5);
    //dx.multiply_val(1,0.5);
    //dx.multiply_val(2,2.0);
    //dx.multiply_val(3,2.0);

    control_integrator integrator(chisq,min,max,dx,"controls/draft_160907/gentle_integrable_detailed");
    array_1d<double> cc;
    cc.add(0.95);
    cc.add(0.68);
    integrator.run_analysis(cc);

    printf("\nrange\n");
    for(i=0;i<4;i++){
        printf("%e %e\n",integrator.get_min(i),integrator.get_max(i));
    }

    printf("that took %e\n",double(time(NULL))-start);
}
