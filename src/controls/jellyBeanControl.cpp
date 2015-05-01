#include "chisq.h"
#include "kd.h"

int main(){

int dim=5;


double angularWidth,radiusOfCurvature;

angularWidth=2.0;
radiusOfCurvature=2.0;

jellyBean chisq(dim,angularWidth,radiusOfCurvature;);

array_1d<double> true_center,widths,curvature_center;
true_center.set_name("true_center");
curvature_center.set_name("curvature_center");
widths.set_name("widths");

int ix,iy;
for(ix=0;ix<dim;ix++){
    true_center.set(ix,chisq.get_real_center(0,ix));
    widths.set(ix,chisq.get_width(0,ix));
}
chisq.get_curvature_center(curvature_center);

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
       
       rowct=0;
       total=0.0;
       for(i=-dd/2;i<dd/2;i++){
           xval=true_center.get_data(ix)+i*dx;
           trial.set(ix,xval);
           for(j=-dd/2;j<dd/2;j++){
               yval=true_center.get_data(iy)+j*dy;
               trial.set(iy,yval);
               
               mu=chisq(trial);
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
