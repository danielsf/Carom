#include "containers.h"

#ifndef GOTO_H
#define GOTO_H
#define pi 3.141592654

#define exception_value 1.0e30

enum{iSimplex,iCompass,iRicochet,iNodeSimplex};

void kill(char*);

double raiseup(double,double);


inline double power(double arg,int raised){

  //return arg raised to the integer power 'raised'

  int n;
  double ans;

  if(raised==0)return 1.0;
  else{ans=1.0;
  for(n=0;n<raised;n++){
    ans=ans*arg;
  }

  return ans;
  }

}

void transcribe(char*,char*);

struct Ran{

//this structure will be based on the Xorshift random number generator
// discovered by George Marsaglia and published in
//Journal of Statistical Software, volume 8, no. 14 pp 1-6

//parameters are drawn from the table on page 347 of 
//Numerical Recipes (3rd edition) 
//William H. press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery
//Cambridge University Press, 2007

unsigned long long x;
Ran(unsigned long long seed){

x=seed^88172645463325252LL;
x^=(x<<21);
x^=(x>>35);
x^=(x<<4);
//printf("starting rand with %ld from seed %d\n",x,seed);
}

void thework(){
  x^=(x<<21);
  x^=(x>>35);
  x^=(x<<4);
}

double doub(){
  thework();
  return x*5.42101086242752217e-20;
}

int int32(){
  thework();
  int ans=int(x);
  if(ans<0)ans=-1*ans;
  return ans;
}

};

void polint(double*,double*,int,double,double*,double*);

double interpolate(double*,double*,double,int);

void sort(double*,int*,int);

void check_sort(double*,int*,int);

void sort_and_check(double*,double*,int*,int);

double normal_deviate(Ran*,double,double);

void naive_gaussian_solver(array_1d<double>&,array_1d<double>&,
array_1d<double>&,int);

double compare_arr(array_1d<double>&,array_1d<double>&);

int compare_int_arr(array_1d<int>&, array_1d<int>&);

struct chisquared_distribution{

    double _fix;

    double _pdf_fn(double x, double dof,double *lnpdf){
        double logans;
        logans=-0.5*x+(0.5*dof-1.0)*log(x)-_fix;

        lnpdf[0]=logans;
        return exp(logans);
    }

    double maximum_likelihood_chisquared_value(double dof){
        _fix=0.0;
        double x,lnpdf,pdf,maxpdf,x_max;
        maxpdf=-1.0;
        for(x=0.0;x<dof+1000.0;x+=1.0){
            pdf=_pdf_fn(x,dof,&lnpdf);
            if(pdf>maxpdf){
                maxpdf=pdf;
                x_max=x;
            }
        }

        return x;
    }

    double confidence_limit(double dof, double pct){
        _fix=0.0;

        double x,pdf,pdfold,xold,total,lim,ans,lnpdf;
        double maxf,maxx,lmax=-1.0e10;
        double start,step,stop,target;
        double x1,x2;

        maxf=-1.0e20;

        target=-1.0;

        start=1.0e-10;
        stop=dof+1000.0;
        step=dof*0.01;

        xold=0.0;
        pdfold=0.0;
        total=0.0;

        _fix=0.0;
        for(x=start;x<stop;x+=step){
            pdf=_pdf_fn(x,dof,&lnpdf);
            if(lnpdf>lmax)lmax=lnpdf;
        }
        _fix=lmax;

        lmax=-1.0e10;

        for(x=start;x<stop;x+=step){
            pdf=_pdf_fn(x,dof,&lnpdf);
            if(lnpdf>lmax){
                maxf=pdf;
                maxx=x;
                lmax=lnpdf;
            }
            total+=0.5*(pdf+pdfold)*(x-xold);
 
            xold=x;
            pdfold=pdf;
        }


        lim=pct*total;
        ans=0.0;
        xold=start;

        pdfold=0.0;
        for(x=start;ans<lim;x+=step){
            pdf=_pdf_fn(x,dof,&lnpdf);
            ans+=0.5*(pdfold+pdf)*(x-xold);
            xold=x;
            pdfold=pdf;
 
        }

        return xold;

}


};

#endif
