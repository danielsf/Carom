#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "wmap_likelihood_function.h"

wmap_likelihood::wmap_likelihood() : chisquared(6){}

wmap_likelihood::~wmap_likelihood(){}

double wmap_likelihood::operator()(array_1d<double> &v) const{

  //this is the function that calls the likelihood function
  //to evaluate chi squared

  //*v is the point in parameter space you have chosen to evaluate

  //as written, it calls the CAMB likelihood function for the CMB anisotropy
  //spectrum

  int i,start,k;
  double params[14],chisquared,omm,base1,base2,amp1,amp2,base;
  double d1,d2,sncc,cc1,cc2,dcc,ccmaxc;
  double cltt[3000],clte[3000],clee[3000],clbb[3000],el[3000];

  double *dir;

  double before=double(time(NULL));


  FILE *output;

  for(i=0;i<dim;i++){
      if((v.get_data(i)<mins.get_data(i) && mins.get_data(i)<exception_value) ||
      (v.get_data(i)>maxs.get_data(i) && maxs.get_data(i)>-1.0*exception_value)){
          time_spent+=double(time(NULL))-before;

          return 2.0*exception_value;
      }
  }

  called++;

  while((v.get_data(0)+v.get_data(1))/(v.get_data(2)*v.get_data(2))>1.0){
    v.multiply_val(0,0.9);
    v.multiply_val(1,0.9);
    v.multiply_val(2,1.1);
    //in the event that total omega_matter>1

  }

  for(i=0;i<6;i++)params[i]=v.get_data(i);

  for(i=0;i<3000;i++)el[i]=0;
  params[2]=100.0*params[2];
  params[5]=exp(params[5])*1.0e-10;


  camb_wrap_(params,el,cltt,clte,clee,clbb); //a function to call CAMB

  for(start=0;el[start]<1.0;start++);

  //printf("cltt start %e %d\n",cltt[start],start);
  if(cltt[start]>=-10.0){
  wmaplikeness_(&cltt[start],&clte[start],&clee[start],\
  &clbb[start],&chisquared); //a function to call the WMAP likelihood code
  }
  else chisquared=2.0*exception_value;

 //printf("done with likelihood\n");

  if(chisquared<0.01)chisquared=2.0*exception_value;
  			//in case the model was so pathological that the
			//likelihood code crashed and returned
			//chisquared=0  (this has been known to happen)




 /////////////////

  time_spent+=double(time(NULL))-before;

  //printf("got chisquared %e\n",chisquared);
  return chisquared;

}

wmap_2d_likelihood::wmap_2d_likelihood() : chisquared(2){}

wmap_2d_likelihood::~wmap_2d_likelihood(){}

double wmap_2d_likelihood::operator()(array_1d<double> &pt) const{

   double before=double(time(NULL));
   int i;
   for(i=0;i<dim;i++){
        if((pt.get_data(i)<mins.get_data(i) && mins.get_data(i)<exception_value) ||
        (pt.get_data(i)>maxs.get_data(i) && maxs.get_data(i)>-1.0*exception_value)){
            time_spent+=double(time(NULL))-before;

            return 2.0*exception_value;
        }
    }


    array_1d<double> vv;
    vv.set_name("wmap_2d_likelihood_operator_vv");
    //vv.set(0, 2.194935109505580934e-02);
    vv.set(0, pt.get_data(0));
    vv.set(1, 1.142206281115405175e-01);
    vv.set(2, 6.871984247221720743e-01);
    vv.set(3, 1.015650744398735594e-02);
    vv.set(4, 9.590415377620740145e-01);
    //vv.set(5, 3.070865283499835119e+00);
    vv.set(5, pt.get_data(1));

    called++;
    time_spent+=double(time(NULL))-before;
    return _wmap(vv);


}
