#include "chisq.h"
#include "kd.h"

int main(){

int dim=5;


double angularWidth,radiusOfCurvature;

angularWidth=1.0;
radiusOfCurvature=20.0;

jellyBean chisq(dim,angularWidth,radiusOfCurvature);

array_1d<double> true_center,widths,curvature_center,radial_direction;
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
array_1d<double> trial;

array_1d<double> chival,chival_sorted;
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

double xval,yval,mu,dx,dy;
double theta,volume,dtheta,dr,rval;
double xmin,xmax,ymin,ymax;
double globalXmin,globalXmax,globalYmin,globalYmax;
double globalDx,globalDy,xMargin;
int rowct;

min.set_dim(2);
max.set_dim(2);
xx.set_cols(2);

kd_tree *allPointsTree,*boundaryPointsTree;

allPointsTree=NULL;
boundaryPointsTree=NULL;

char filename[3*letters];
FILE *output;

for(ix=0;ix<dim;ix++){
    for(iy=ix+1;iy<dim;iy++){
       if(allPointsTree!=NULL){
           delete allPointsTree;
           allPointsTree=NULL;
       } 
       if(boundaryPointsTree!=NULL){
           delete boundaryPointsTree;
           boundaryPointsTree=NULL;
       }
       total=0.0;
       xx.reset_preserving_room();
       xxTree.reset_preserving_room();
       boundary.reset_preserving_room();
       dexes.reset_preserving_room();
       chival.reset_preserving_room();
       chival_sorted.reset_preserving_room();
       
       for(i=0;i<dim;i++){
           trial.set(i,true_center.get_data(i));
       }
       

       dx=4.0*widths.get_data(ix)/double(dd/2);
       dy=4.0*widths.get_data(iy)/double(dd/2);
       
       xmin=true_center.get_data(ix)-(dd/2)*dx;
       xmax=true_center.get_data(ix)+(dd/2)*dx;
       
       ymin=true_center.get_data(iy)-(dd/2)*dy;
       ymax=true_center.get_data(iy)+(dd/2)*dy;

       if(ix==0 && iy==1){
           dx=widths.get_data(1)*4.0/double(dd/2);
           dy=widths.get_data(1)*4.0/double(dd/2);
           
           xmin=2.0*exception_value;
           xmax=-2.0*exception_value;
           ymin=2.0*exception_value;
           ymax=-2.0*exception_value;
           
           for(theta=-4.0*angularWidth;theta<4.0*angularWidth;theta+=0.01){
               for(rval=radiusOfCurvature-(dd/2)*dx;rval<radiusOfCurvature+(dd/2)*dx;rval+=dx){
                   xval=curvature_center.get_data(0)+rval*(
                        radial_direction.get_data(0)*cos(theta)-
                        radial_direction.get_data(1)*sin(theta));
               
                   yval=curvature_center.get_data(1)+rval*(
                        radial_direction.get_data(0)*sin(theta)+
                        radial_direction.get_data(1)*cos(theta));
               
                   if(xval<xmin)xmin=xval;
                   if(xval>xmax)xmax=xval;
                   if(yval<ymin)ymin=yval;
                   if(yval>ymax)ymax=yval;
               }
           }
       
           dy=(ymax-ymin)/double(400);
           dx=(xmax-xmin)/double(400);
       
           globalXmin=xmin;
           globalXmax=xmax;
           globalYmin=ymin;
           globalYmax=ymax;
           
           globalDx=dx;
           globalDy=dy;
       
       }
       else if(ix==0){
           xmin=globalXmin;
           xmax=globalXmax;
           dx=globalDx;
       }
       else if(ix==1){
           xmin=globalYmin;
           xmax=globalYmax;
           dx=globalDy;
       }
       
       rowct=0;
       total=0.0;
       
       
       printf("\nx %e %e %e y %e %e %e\n",xmin,xmax,dx,ymin,ymax,dy);

       for(xval=xmin;xval<xmax;xval+=dx){
           trial.set(ix,xval);
           for(yval=ymin;yval<ymax;yval+=dy){
               trial.set(iy,yval);
               
               
               if(ix>=2 || (ix==0 && iy==1)){
                   mu=chisq(trial);
               }
               else if(ix==0){
                   mu=0.0;
                   for(xMargin=globalYmin;xMargin<globalYmax;xMargin+=globalDy){
                       trial.set(1,xMargin);
                       mu+=exp(-0.5*chisq(trial))*globalDy;
                   }
                   
                   mu=-2.0*log(mu);
                   trial.set(1,true_center.get_data(1));
               }
               else if(ix==1){
                   mu=0.0;
                   for(xMargin=globalXmin;xMargin<globalXmax;xMargin+=globalDx){
                       trial.set(0,xMargin);
                       mu+=exp(-0.5*chisq(trial))*globalDx;
                   }
                   
                   mu=-2.0*log(mu);
                   trial.set(0,true_center.get_data(0));
               }
               
               total+=exp(-0.5*mu);
               chival.add(mu);
               xx.set(rowct,0,xval);
               xx.set(rowct,1,yval);
               dexes.add(rowct);
               rowct++;
           
           }
       }

       printf("rowct %d\n",rowct);
       sort_and_check(chival,chival_sorted,dexes);
       sum=0.0;
       for(i=0;i<xx.get_rows() && sum<0.95*total;i++){
           j=dexes.get_data(i);
           sum+=exp(-0.5*chival.get_data(j));
           xxTree.add_row(xx(j)[0]);
       }
       
       if(i>=xx.get_rows())i=xx.get_rows()-1;
       
       printf("%d %d total %e sum %e -- %e %d %e\n",
       ix,iy,total,sum,0.95*total,xxTree.get_rows(),
       chival_sorted.get_data(i));
      
       min.set(0,0.0);
       min.set(1,0.0);
       max.set(0,dx);
       max.set(1,dy);
       
       allPointsTree=new kd_tree(xxTree,min,max);
       for(i=0;i<xxTree.get_rows();i++){
           allPointsTree->nn_srch(xxTree(i)[0],5,neigh,dist);
           if(dist.get_data(4)>1.001){
               boundary.add_row(xxTree(i)[0]);
           }
       }
       
       printf("got boundary points %d\n",boundary.get_rows());
       
       boundaryPointsTree=new kd_tree(boundary);
       origin.set(0,boundaryPointsTree->get_pt(0,0));
       origin.set(1,boundaryPointsTree->get_pt(0,1));
       lastpt.set(0,origin.get_data(0));
       lastpt.set(1,origin.get_data(1));
       
       lastdex=0;
       sprintf(filename,"controls/jellyBean_%d_%d_control.txt",ix,iy);
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
       
    }
}


}
