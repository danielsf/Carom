#include "chisq.h"

int main(){

int dim=5;

jellyBean chisq(dim,2.0,2.0);

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

}
