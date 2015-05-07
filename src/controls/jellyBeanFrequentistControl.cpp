#include "chisq.h"
#include "kd.h"

array_1d<double> curvature_center, radial_direction;
double globalXmin,globalXmax,globalYmin,globalYmax;

void simple_bisection(array_1d<double> &lowball_in, array_1d<double> &highball_in, double target,
                      chisquared *fn, array_1d<double> &output){

    array_1d<double> lowball,highball,trial;
    lowball.set_name("simple_bisection_lowball");
    highball.set_name("simple_bisection_highball");
    trial.set_name("simple_bisection_trial");

    double flow,fhigh,ftrial,fbest;
    
    fbest=2.0*exception_value;
    
    int i;
    for(i=0;i<fn->get_dim();i++){
        lowball.set(i,lowball_in.get_data(i));
        highball.set(i,highball_in.get_data(i));
        output.set(i,lowball.get_data(i));
    }
    
    flow=fn[0](lowball);
    flow=fn[0](highball);
    fbest=flow;
    
    int ct;
    for(ct=0;ct<100;ct++){
        for(i=0;i<fn->get_dim();i++){
            trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
        }
        
        ftrial=fn[0](trial);
        if(ftrial<target){
            for(i=0;i<fn->get_dim();i++){
                lowball.set(i,trial.get_data(i));
            }
        }
        else{
            for(i=0;i<fn->get_dim();i++){
                highball.set(i,trial.get_data(i));
            }
        }
        
        if(fabs(ftrial-target)<fabs(fbest-target)){
            fbest=ftrial;
            for(i=0;i<fn->get_dim();i++){
                output.set(i,trial.get_data(i));
            }
        }
    }
    
}

double marginalized_call(array_1d<double> &pt_in, int iMargin, chisquared *fn){

    int i;
    array_1d<double> pt;
    pt.set_name("marginalized_call_pt");
    
    for(i=0;i<pt_in.get_dim();i++){
        pt.set(i,pt_in.get_data(i));
    }

    double chisq,chisq_best;
    
    chisq_best=2.0*exception_value;
    double xx,dx;
    
    double minval,maxval;
    if(iMargin==0){
        minval=globalXmin;
        maxval=globalXmax;
    }
    else{
        minval=globalYmin;
        maxval=globalYmax;
    }
    
    dx=0.1;
    
    for(xx=minval;xx<maxval;xx+=dx){
        pt.set(iMargin,xx);
        chisq=fn[0](pt);
        if(chisq<chisq_best){
            chisq_best=chisq;
        }
    }

    return chisq_best;

}


void marginalized_bisection(array_1d<double> &lowball_in, array_1d<double> &highball_in, double target,
                      chisquared *fn, array_1d<double> &output, int iMargin){

    array_1d<double> lowball,highball,trial;
    lowball.set_name("simple_bisection_lowball");
    highball.set_name("simple_bisection_highball");
    trial.set_name("simple_bisection_trial");

    double flow,fhigh,ftrial,fbest;
    
    fbest=2.0*exception_value;
    
    int i;
    for(i=0;i<fn->get_dim();i++){
        lowball.set(i,lowball_in.get_data(i));
        highball.set(i,highball_in.get_data(i));
        output.set(i,lowball.get_data(i));
    }
    
    flow=marginalized_call(lowball,iMargin,fn);
    flow=marginalized_call(highball,iMargin,fn);
    fbest=flow;
    
    int ct;
    for(ct=0;ct<100;ct++){

        for(i=0;i<fn->get_dim();i++){
            trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
        }
        
        ftrial=marginalized_call(trial,iMargin,fn);
        if(ftrial<target){
            for(i=0;i<fn->get_dim();i++){
                lowball.set(i,trial.get_data(i));
            }
        }
        else{
            for(i=0;i<fn->get_dim();i++){
                highball.set(i,trial.get_data(i));
            }
        }
        
        if(fabs(ftrial-target)<fabs(fbest-target)){
            fbest=ftrial;
            for(i=0;i<fn->get_dim();i++){
                output.set(i,trial.get_data(i));
            }
        }
    }
    
}


int main(){

int dim=15;
double chi_target=25.0;

array_1d<int> xdexes,ydexes;

xdexes.set_name("xdexes");
ydexes.set_name("ydexes");

xdexes.add(0);
xdexes.add(0);
xdexes.add(1);
xdexes.add(6);

ydexes.add(1);
ydexes.add(5);
ydexes.add(3);
ydexes.add(8);

char outNameRoot[letters];
sprintf(outNameRoot,"controls/jellyBeanD15control/jellyBeanD15freq");

double angularWidth,radiusOfCurvature;

angularWidth=1.0;
radiusOfCurvature=20.0;

jellyBean chisq(dim,angularWidth,radiusOfCurvature);


array_1d<double> true_center,widths;
true_center.set_name("true_center");
curvature_center.set_name("curvature_center");
widths.set_name("widths");
radial_direction.set_name("radial_direction");

int ix,iy;
for(ix=0;ix<dim;ix++){
    true_center.set(ix,chisq.get_real_center(0,ix));
    widths.set(ix,chisq.get_width(0,ix));
}
chisq.get_curvature_center(curvature_center);

for(ix=0;ix<2;ix++){
    radial_direction.set(ix,true_center.get_data(ix)-curvature_center.get_data(ix));
}
radial_direction.normalize();

printf("radial width %e\n",widths.get_data(1));

double total,sum;

array_1d<double> chival,chival_sorted,trial;
array_1d<double> min,max,dist,origin,lastpt;
array_2d<double> xx,xxTree,boundary;
array_1d<int> dexes,neigh;

trial.set_name("trial");
dexes.set_name("dexes");
chival.set_name("chival");
chival_sorted.set_name("chival_sorted");
min.set_name("min");
max.set_name("max");
xx.set_name("xx");
xxTree.set_name("xxTree");
boundary.set_name("boundary");
neigh.set_name("neigh");
dist.set_name("dist");
origin.set_name("origin");
lastpt.set_name("lastpt");

int i,j,dd,lastdex;

dd=100;

double xval,yval,mu,dx,dy,marginval;
double theta,volume,dtheta,rval,rr;
double xmin,xmax,ymin,ymax;
double globalDx,globalDy,xMargin;
int rowct;

min.set_dim(2);
max.set_dim(2);
xx.set_cols(2);

kd_tree *boundaryPointsTree;

boundaryPointsTree=NULL;


if(xdexes.get_dim()==0){
    for(ix=0;ix<dim;ix++){
        for(iy=ix+1;iy<dim;iy++){
            xdexes.add(ix);
            ydexes.add(iy);
        }
    }
}

char filename[3*letters];
FILE *output;

int crossed_middle;
array_1d<double> best_pt,highball,lowball;
double chi_best,sgn,dr;
best_pt.set_name("best_pt");
highball.set_name("highball");
lowball.set_name("lowball");
trial.set_name("trial");

dtheta=0.01;

double chi_worst,minval,maxval;
int ii,ct,iMargin,iVar;
boundary.set_cols(2);
for(ii=0;ii<xdexes.get_dim();ii++){
    ix=xdexes.get_data(ii);
    iy=ydexes.get_data(ii);
    
    if(ix>iy){
        i=ix;
        ix=iy;
        iy=i;
    }
    
    boundary.reset_preserving_room();

      printf("about to delete tree %d %d\n",ix,iy);
       if(boundaryPointsTree!=NULL){
           delete boundaryPointsTree;
           boundaryPointsTree=NULL;
       }
       printf("done with that\n");
      
       if(ix>=2 && iy>=2){
           dtheta=0.01;
           ct=0;
           for(theta=0.0;theta<=2.0*pi;theta+=dtheta){
               mu=-2.0*exception_value;
               for(i=0;i<dim;i++){
                   trial.set(i,true_center.get_data(i));
               }
               
               while(mu>chi_target){
                   trial.add_val(ix,cos(theta));
                   trial.add_val(iy,sin(theta));
                   mu=chisq(trial);
               }
           
               simple_bisection(true_center,trial,chi_target,&chisq,best_pt);
               boundary.set(ct,0,best_pt.get_data(ix));
               boundary.set(ct,1,best_pt.get_data(iy));
               ct++;
           }
       }
       else if(ix<2 && iy<2){
           globalXmin=2.0*exception_value;
           globalXmax=-2.0*exception_value;
           globalYmin=2.0*exception_value;
           globalYmax=-2.0*exception_value;
           
           dtheta=0.01;
           ct=0;
           chi_worst=chi_target;
           for(theta=0.0;theta<2.0*pi;theta+=dtheta){
               for(i=0;i<dim;i++){
                   trial.set(i,true_center.get_data(i));
               }
               trial.set(ix,curvature_center.get_data(0)+radiusOfCurvature*cos(theta));
               trial.set(iy,curvature_center.get_data(1)+radiusOfCurvature*sin(theta));
               mu=chisq(trial);
               
               if(mu<=chi_target){
                   for(sgn=-1.0;sgn<1.1;sgn+=2.0){
                       for(i=0;i<dim;i++){
                           highball.set(i,trial.get_data(i));
                       }
                       mu=-2.0*exception_value;
                       while(mu<chi_target){
                           highball.add_val(ix,sgn*cos(theta));
                           highball.add_val(iy,sgn*sin(theta));
                           mu=chisq(highball);
                       }
                       
                       simple_bisection(trial,highball,chi_target,&chisq,best_pt);
                       mu=chisq(best_pt);
                       if(fabs(mu-chi_target)>fabs(chi_worst-chi_target)){
                           chi_worst=mu;
                       }
                       boundary.set(ct,0,best_pt.get_data(ix));
                       boundary.set(ct,1,best_pt.get_data(iy));
                       ct++;
                   
                       if(best_pt.get_data(ix)>globalXmax){
                           globalXmax=best_pt.get_data(ix);
                       }
                       if(best_pt.get_data(ix)<globalXmin){
                           globalXmin=best_pt.get_data(ix);
                       }
                       
                       if(best_pt.get_data(iy)>globalYmax){
                           globalYmax=best_pt.get_data(iy);
                       }
                       if(best_pt.get_data(iy)<globalYmin){
                           globalYmin=best_pt.get_data(iy);
                       }
                   
                   }
               }

           }
           printf("chi_worst on 0,1 %e\n",chi_worst);
           
       }
       else{
           chi_worst=25.0;
           ct=0;
           if(ix==0){
               iMargin=1;
               minval=globalXmin;
               maxval=globalXmax;
           }
           else if(ix==1){
               iMargin=0;
               minval=globalYmin;
               maxval=globalYmax;
           }
           
           for(i=0;i<dim;i++){
               trial.set(i,true_center.get_data(i));
           }
           

           
           dtheta=0.01;
           for(theta=0.0;theta<2.0*pi;theta+=dtheta){

               for(i=0;i<dim;i++){
                   highball.set(i,true_center.get_data(i));
               }
               
               mu=-2.0*exception_value;
               while(mu<chi_target){
                   highball.add_val(ix,cos(theta));
                   highball.add_val(iy,sin(theta));
                   mu=marginalized_call(highball,iMargin,&chisq);
               }
               
               marginalized_bisection(true_center,highball,chi_target,&chisq,best_pt,iMargin);
               
               mu=marginalized_call(best_pt,iMargin,&chisq);
               if(fabs(mu-chi_target)>fabs(chi_worst-chi_target)){
                   chi_worst=mu;
                   printf("    chi_worst %e\n",mu);
               }
               
               boundary.set(ct,0,best_pt.get_data(ix));
               boundary.set(ct,1,best_pt.get_data(iy));
               
               if((best_pt.get_data(iy)-true_center.get_data(iy))*sin(theta)<0.0){
                   printf("WARNING sin %e but d %e\n",
                   sin(theta),best_pt.get_data(iy)-true_center.get_data(iy));
                   exit(1);
               }
               
               ct++;

               
           }
           
           printf("chi_worst on %d %d %e\n",ix,iy,chi_worst);
           
       
       }
       
       printf("got boundary points %d\n",boundary.get_rows());
       
       min.set(0,2.0*exception_value);
       max.set(0,-2.0*exception_value);
       min.set(1,2.0*exception_value);
       max.set(1,-2.0*exception_value);
       
       for(i=0;i<boundary.get_rows();i++){
           if(boundary.get_data(i,0)<min.get_data(0)){
               min.set(0,boundary.get_data(i,0));
           }
           
           if(boundary.get_data(i,0)>max.get_data(0)){
               max.set(0,boundary.get_data(i,0));
           }
           
           if(boundary.get_data(i,1)<min.get_data(1)){
               min.set(1,boundary.get_data(i,1));
           }
           
           if(boundary.get_data(i,1)>max.get_data(1)){
                max.set(1,boundary.get_data(i,1));
           }
       }
       
       
       boundaryPointsTree=new kd_tree(boundary,min,max);
       origin.set(0,boundaryPointsTree->get_pt(0,0));
       origin.set(1,boundaryPointsTree->get_pt(0,1));
       lastpt.set(0,origin.get_data(0));
       lastpt.set(1,origin.get_data(1));
       
       lastdex=0;
       sprintf(filename,"%s_%d_%d_control.txt",outNameRoot,ix,iy);
       output=fopen(filename,"w");
       while(boundaryPointsTree->get_pts()>1){
           fprintf(output,"%e %e\n",lastpt.get_data(0),lastpt.get_data(1));
           boundaryPointsTree->remove(lastdex);
           
           boundaryPointsTree->nn_srch(lastpt,1,neigh,dist);
           lastdex=neigh.get_data(0);
           lastpt.set(0,boundaryPointsTree->get_pt(lastdex,0));
           lastpt.set(1,boundaryPointsTree->get_pt(lastdex,1));
       }
       fprintf(output,"%e %e\n",lastpt.get_data(0),lastpt.get_data(1));
       fprintf(output,"%e %e\n",origin.get_data(0),origin.get_data(1));
       
       fclose(output);
       printf("done printing\n");
    
}

for(i=0;i<dim;i++){
    printf("center %d %e\n",i,true_center.get_data(i));
}


}
