#include <stdio.h>
#include <stdlib.h>
#include "wmap_likelihood_function.h"

double target = 1280.7;
wmap_2d_likelihood chisq;

void bisect(array_1d<double> &lowball, double flow,
            array_1d<double> &highball, double fhigh,
            array_1d<double> &best){

    int i;
    double ftest,dfbest,distance;
    array_1d<double> trial;
    trial.set_name("trial");
    distance=10.0;

    if(fabs(flow-target)<fabs(fhigh-target)){
        dfbest=fabs(flow-target);
        best.set(0,lowball.get_data(0));
        best.set(1,lowball.get_data(1));
    }
    else{
        dfbest=fabs(fhigh-target);
        best.set(0,highball.get_data(0));
        best.set(1,highball.get_data(1));
    }

    while(distance>1.0e-6){
        for(i=0;i<2;i++){
            trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
        }
        ftest=chisq(trial);

        if(ftest<target){
            lowball.set(0,trial.get_data(0));
            lowball.set(1,trial.get_data(1));
        }
        else{
            highball.set(0,trial.get_data(0));
            highball.set(1,trial.get_data(1));
        }

        if(fabs(ftest-target)<dfbest){
            best.set(0,trial.get_data(0));
            best.set(1,trial.get_data(1));
            dfbest=fabs(ftest-target);
        }

        distance=power(lowball.get_data(0)-highball.get_data(0),2);
        distance+=power(highball.get_data(1)-highball.get_data(1),2);
        distance=sqrt(distance);

    }


}

int main(){

array_1d<double> center;
center.set_name("center");
center.set(0, 2.194935109505580934e-02);
center.set(1, 3.070865283499835119e+00);

array_1d<double> lowball,highball,dir,best;
lowball.set_name("lowball");
highball.set_name("highball");
dir.set_name("dir");
best.set_name("best");

double flow,fhigh,fbest,theta,fcenter,rr;
fcenter=chisq(center);

FILE *output;

for(theta=0.0;theta<2.0*pi;theta+=0.1*pi){
    lowball.set(0,center.get_data(0));
    lowball.set(1,center.get_data(1));
    highball.set(0,center.get_data(0));
    highball.set(1,center.get_data(1));
    flow=fcenter;

    dir.set(0,cos(theta));
    dir.set(1,sin(theta));
    fhigh=-2.0*exception_value;
    rr=1.0;
    while(fhigh<target){
        highball.add_val(0,rr*dir.get_data(0));
        highball.add_val(1,rr*dir.get_data(1));
        fhigh=chisq(highball);
        rr*=2.0;

    }

    bisect(lowball,flow,highball,fhigh,best);
    fbest=chisq(best);
    output=fopen("output/wmap2d_control_output.sav","a");
    fprintf(output,"%e %e %e\n",best.get_data(0),best.get_data(1),fbest);
    fclose(output);

}

}
