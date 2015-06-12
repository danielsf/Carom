#include "jellyBean.h"

jellyBean::~jellyBean(){}

jellyBean::jellyBean(int id, double ww, double rr) : chisquared(id){
    //ww is a width in radians
    //rr is the radius of curvature
    make_bases(22,0);
    
    int ix,iy;
    
    for(ix=0;ix<_dim;ix++){
        for(iy=0;iy<_dim;iy++){
            if(ix==iy){
                _bases.set(ix,iy,1.0);
            }
            else{
                _bases.set(ix,iy,0.0);
            }
        }
    }
    
    _time_spent=0.0;
    _called=0;
    
    //these are in basis coordinates
    _curvature_center.set_name("jellyBean_curvature_center");
    _radial_direction.set_name("jellyBean_radial_direction");
    
    for(ix=0;ix<_dim;ix++){
        _centers.set(0,ix,-10.0+20.0*_dice->doub());
        if(ix>0){
            _widths.set(0,ix,_dice->doub()*0.5+0.1);
        }
    }
    
    _widths.set(0,0,ww*rr);
    

    
    double theta=_dice->doub()*2.0*pi;
    double dx,dy;
    
    _curvature_radius=rr;
    
    dx=_curvature_radius*cos(theta);
    dy=_curvature_radius*sin(theta);
    
    array_1d<double> projected_center;
    projected_center.set_name("jellyBean_constructor_projected_center");
    
    for(ix=0;ix<_dim;ix++){
        projected_center.set(ix,project_to_basis(ix,_centers(0)[0]));
    }
    
    for(ix=0;ix<_dim;ix++){
        _curvature_center.set(ix,projected_center.get_data(ix));
    }
    
    _curvature_center.add_val(0,dx);
    _curvature_center.add_val(1,dy);
    
    _radial_direction.set_dim(_dim);
    _radial_direction.zero();
    _radial_direction.set(0,-1.0*dx);
    _radial_direction.set(1,-1.0*dy);
    
    _radial_direction.normalize();
    
}

void jellyBean::get_curvature_center(array_1d<double> &out){
    int ix;
    for(ix=0;ix<_dim;ix++){
        out.set(ix,project_to_basis(ix,_curvature_center));
    }
}

double jellyBean::operator()(array_1d<double> &pt){
    double before=double(time(NULL));
    array_1d<double> projected_point;
    
    projected_point.set_name("jellyBean_operator_projected_point");
    
    int ix;
    for(ix=0;ix<_dim;ix++){
        projected_point.set(ix, project_to_basis(ix,pt));
    }

    double chisq=0.0;
    for(ix=2;ix<_dim;ix++){
        chisq+=power((_centers.get_data(0,ix)-projected_point.get_data(ix))/_widths.get_data(0,ix),2);
    }
    
    if(isnan(chisq)){
        printf("chisq is nan after simple dims\n");
        exit(1);
    }
    
    double rr;
    array_1d<double> dir;
    dir.set_name("jellyBean_operator_dir");
    dir.set_dim(_dim);
    dir.zero();
    dir.set(0,projected_point.get_data(0)-_curvature_center.get_data(0));
    dir.set(1,projected_point.get_data(1)-_curvature_center.get_data(1));
    rr=dir.normalize();
    
    chisq+=power((rr-_curvature_radius)/_widths.get_data(0,1),2);
    
    //printf("chisq is %e after rr\n",chisq);
    
    if(isnan(chisq)){
        printf("chisq is nan after radial %e\n",rr);
        exit(1);
    }
    
    double dot=0.0;
    for(ix=0;ix<_dim;ix++){
        dot+=dir.get_data(ix)*_radial_direction.get_data(ix);
    }
    
    double theta=acos(dot);
    
    if(isnan(theta) && dot>0.0){
        theta=0.0;
    }
    else if(isnan(theta) && dot<0.0){
        theta=pi;
    }
    
    chisq+=power(theta*_curvature_radius,2)/_widths.get_data(0,0);
    
    //printf("chisq is %e after angle\n",chisq);
    
    _called++;
    _time_spent+=double(time(NULL))-before;
    
    return chisq;

}
