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

////////////////////////////////jellyBeanData

chiSquaredData::chiSquaredData(int dd, int cc, int nData, double sigma) : chisquared(dd, cc){
    printf("_dim %d\n",_dim);

    _x_values.set_name("jellyBeanData_x");
    _y_values.set_name("jellyBeanData_y");
    _sigma.set_name("jellyBeanData_sigma");
    _mean_parameters.set_name("jellyBeanData_mean_parameters");
    _aux_params.set_name("jellyBeanData_aux_parameters");
    _param_buffer.set_name("jellyBeanData_param_buffer");
    
    make_bases(22,0);
    
    int ii,jj,goon,i;
    double nn,normerr,ortherr,theta;
    
    if(_dim>4){
        
        theta=_dice->doub()*2.0*pi;
        for(ii=0;ii<_dim;ii++){
            _bases.set(0,ii,0.0);
            _bases.set(1,ii,0.0);
        }
        
        _bases.set(0,0,cos(theta));
        _bases.set(0,1,sin(theta));
        _bases.set(1,0,-sin(theta));
        _bases.set(1,1,cos(theta));
        
        for(ii=2;ii<_dim;ii++){
            goon=1;
	    while(goon==1){
                goon=0;
	        for(i=0;i<_dim;i++)_bases.set(ii,i,_dice->doub()-0.5);
	        for(jj=0;jj<ii;jj++){
	            nn=0.0;
		    for(i=0;i<_dim;i++)nn+=_bases.get_data(ii,i)*_bases.get_data(jj,i);
		    for(i=0;i<_dim;i++)_bases.subtract_val(ii,i,nn*_bases.get_data(jj,i));
	        }
	    
	        nn=0.0;
	        for(i=0;i<_dim;i++){
		    nn+=power(_bases.get_data(ii,i),2);
	        }
	        if(nn<1.0e-20)goon=1;
	        nn=sqrt(nn);
	        for(i=0;i<_dim;i++){
	            _bases.divide_val(ii,i,nn);
	        }
	    }
        }
   
        for(ii=0;ii<_dim;ii++){
            nn=0.0;
	    for(i=0;i<_dim;i++)nn+=power(_bases.get_data(ii,i),2);
	    nn=fabs(1.0-nn);
	    if(ii==0 || nn>normerr)normerr=nn;
	
	    for(jj=ii+1;jj<_dim;jj++){
	       nn=0.0;
	       for(i=0;i<_dim;i++)nn+=_bases.get_data(ii,i)*_bases.get_data(jj,i);
	       nn=fabs(nn);
	       if((ii==0 && jj==1) || nn>ortherr)ortherr=nn;
	    }
        }
    
        printf("normerr %e ortherr %e\n",normerr,ortherr);
        if(normerr>1.0e-3 || ortherr>1.0e-3){
            death_knell("normerr or ortherr too large");
        }
    }
    
    
    _ndata=nData;
    _sig=sigma;
    
    int ix,ic;
        
    for(ic=0;ic<_ncenters;ic++){ 
        for(ix=0;ix<_dim;ix++){
            _widths.set(ic,ix,fabs(normal_deviate(_dice,0.1,0.5))+0.2);
        }
    }    
        
    for(ix=0;ix<3;ix++){
        _mean_parameters.set(ix,1.0);
    }
    
    int aux_ct;
    double norm;
    aux_ct=0;
    for(;ix<_dim;ix++){
        _mean_parameters.set(ix,_dice->doub()*4.0-2.0);
        
        norm=exp(log(10.0)*(aux_ct/2));
        
        _aux_params.set(aux_ct,_dice->doub()*norm);
        aux_ct++;
        _aux_params.set(aux_ct,_dice->doub()*20.0+26.0);
        aux_ct++;
        
 
    }
    //printf("aux_ct %d %e %e\n",aux_ct,_mean_parameters.get_data(3),_mean_parameters.get_data(4));
    for(ix=0;ix<aux_ct;ix++){
        printf("aux %e\n",_aux_params.get_data(ix));
    }
}


void chiSquaredData::write_data(){
    FILE *output;
    output=fopen("data_scratch.txt", "w");
    int ix;
    for(ix=0;ix<_x_values.get_dim();ix++){
        fprintf(output,"%e %e %e\n",_x_values.get_data(ix),_y_values.get_data(ix),_sigma.get_data(ix));
    }
    fclose(output);
}

void chiSquaredData::print_mins(){
    array_1d<double> trial,params;
    trial.set_name("jellyBeanData_print_mins_trial");
    params.set_name("jellyBeanData_print_mins_params");
    int ic,ix;
    for(ic=0;ic<_ncenters;ic++){
        for(ix=0;ix<_dim;ix++){
            trial.set(ix,_centers.get_data(ic,ix));
        } 
        printf("center %d val %e\n",ic,this[0](trial));
        convert_params(trial,params,ic);
        for(ix=0;ix<_dim;ix++){
            printf("    %e -- width %e -- mean %e \n",
            params.get_data(ix),_widths.get_data(0,ix),_mean_parameters.get_data(ix));
        }
    }
    printf("x_vals %d\n",_x_values.get_dim());
    for(ix=0;ix<_aux_params.get_dim();ix++){
        printf("    aux %d %e\n",ix,_aux_params.get_data(ix));
    }
    
    printf("\nbases\n");
    for(ix=0;ix<_dim;ix++){
        for(ic=0;ic<_dim;ic++){
            printf("%.3e ",_bases.get_data(ic,ix));
        }
        printf("\n");
    }

}

void chiSquaredData::initialize_data(){

    array_1d<double> params;
    params.set_name("jellyBeanData_initiailize_params");
    int ix,ic;
    
    for(ix=0;ix<_dim;ix++)params.set(ix,0.0);
    
    array_1d<double> samples;
    samples.set_name("jellyBeanData_initialize_samples");

    double xx,yy,mean,var;
    for(ix=0, xx=0.0;xx<3.0;xx+=0.03, ix++){
        mean=data_function(params,xx);
        for(ic=0;ic<_ndata;ic++){
            samples.set(ic,normal_deviate(_dice,mean,_sig));
        }
        mean=0.0;
        for(ic=0;ic<samples.get_dim();ic++){
            mean+=samples.get_data(ic);
        }
        mean=mean/double(samples.get_dim());
        _x_values.set(ix,xx);
        _y_values.set(ix,mean);
        
        var=0.0;
        for(ic=0;ic<samples.get_dim();ic++){
            var+=power(samples.get_data(ic)-mean,2);
        }
        var=var/double((samples.get_dim()-1)*samples.get_dim());
        _sigma.set(ix,sqrt(var));
    }
    
    write_data();

}

void chiSquaredData::convert_params(array_1d<double> &pt, array_1d<double> &out, int ic){
    printf("Called void convert_params");
    exit(1);
}

double chiSquaredData::data_function(array_1d<double> &params, double xx){
    
    double ans=0.0;

    ans=(_mean_parameters.get_data(0)+params.get_data(0))*xx*xx;
    ans+=(_mean_parameters.get_data(1)+params.get_data(1))*xx;
    ans+=_mean_parameters.get_data(2)+params.get_data(2);

    double mu,amp;
    int ix=3,i1,i2;
    i1=0;
    i2=1;
    for(;ix<_dim;ix++){
        amp=_mean_parameters.get_data(ix)+params.get_data(ix);
        mu=amp*sin(_aux_params.get_data(i2)*(xx-_aux_params.get_data(i1)));
        ans+=mu;
        i1+=2;
        i2+=2;
    }
    
    return ans;

}

double chiSquaredData::operator()(array_1d<double> &pt){

    double before=double(time(NULL));

    if(_x_values.get_dim()==0){
        initialize_data();
    }

    double chisq,chisq_min;
    double yy;
    int ic,ix;
    for(ic=0;ic<_ncenters;ic++){
        chisq=0.0;
        convert_params(pt,_param_buffer, ic);
        for(ix=0;ix<_x_values.get_dim();ix++){
            yy=data_function(_param_buffer,_x_values.get_data(ix));
            chisq+=power((yy-_y_values.get_data(ix))/_sigma.get_data(ix),2);
        }
        
        if(ic==0 || chisq<chisq_min){
            chisq_min=chisq;
        }
    }
    
    _called++;
    _time_spent+=double(time(NULL))-before;
    
    if(isnan(chisq_min)){
        printf("chisq is nan\n");
        for(ix=0;ix<_param_buffer.get_dim();ix++){
            printf("%e\n",_param_buffer.get_data(ix));
        }
        exit(1);
    }
    
    return chisq_min;
}


jellyBeanData::jellyBeanData(int dd, int cc, int nData, double sigma, double ww, double radial, double rr) : 
chiSquaredData(dd, cc, nData, sigma){

    /*
    dim
    ncenters
    nData points
    spread of data points
    angular width
    radial width (fraction)
    curvature radius
    */

    _curvature_centers.set_name("jellyBeanData_curvature_center");
    _radial_directions.set_name("jellyBeanData_radial_directions");
    _radii.set_name("jellyBeanData_radii");
    _dir.set_name("jellyBeanData_global_dir");
    _planar_dir.set_name("jellyBeanData_planar_dir");

    _curvature_centers.set_dim(_ncenters,_dim);
    _radial_directions.set_dim(_ncenters,_dim);
    _planar_dir.set_dim(_dim);
    
    array_1d<double> dir,trial;
    dir.set_name("jellyBeanData_constructor_dir");
    trial.set_name("jellyBeanData_constructor_trial");
    
    double curvature_radius=rr;
    
    double mu,dot_check;
    int ix,iy,ic;
    for(ic=0;ic<_ncenters;ic++){
        for(ix=0;ix<_dim;ix++){
            mu=-4.0+_dice->doub()*8.0;
            _centers.set(ic,ix,mu);
            _curvature_centers.set(ic,ix,mu);
        }
        
        dir.set_dim(_dim);
        dir.zero();
        
        for(ix=0;ix<2;ix++){
            mu=normal_deviate(_dice,0.0,1.0);
            for(iy=0;iy<_dim;iy++){
                dir.add_val(iy,_bases.get_data(ix,iy));
            }
        }
        
        dir.normalize();
        
        for(ix=0;ix<_dim;ix++){
            _curvature_centers.add_val(ic,ix,curvature_radius*dir.get_data(ix));
        }
        
        for(ix=0;ix<_dim;ix++){
            _radial_directions.set(ic,ix,_centers.get_data(ic,ix)-_curvature_centers.get_data(ic,ix));
        }
        
        _radii.set(ic,_radial_directions(ic)->normalize());
        
        _widths.set(ic,0,ww);
        
        _widths.set(ic,1,radial*curvature_radius);
        
        
        for(ix=2;ix<_dim;ix++){
            dot_check=0.0;
            for(iy=0;iy<_dim;iy++){
                dot_check+=_radial_directions.get_data(ic,iy)*_bases.get_data(ix,iy);
            }
            if(fabs(dot_check)>0.001){
                printf("WARNING dot_check %d %e\n",ix,dot_check);
                exit(1);
            }
        }
        

    }
    
    
}


void jellyBeanData::convert_params(array_1d<double> &pt, array_1d<double> &out, int ic){

    int ix,iy;
    double mu;
    
    for(ix=0;ix<_dim;ix++){
        _dir.set(ix,pt.get_data(ix)-_curvature_centers.get_data(ic,ix));
    }
    
    for(ix=0;ix<_dim;ix++)_planar_dir.set(ix,0.0);
    for(ix=0;ix<2;ix++){
        mu=project_to_basis(ix,_dir);
        for(iy=0;iy<_dim;iy++){
            _planar_dir.add_val(iy,mu*_bases.get_data(ix,iy));
        }
    }
    
    double planar_radius=_planar_dir.normalize();
    
    double cos_theta=0.0,sin_theta;
    for(ix=0;ix<_dim;ix++){
        cos_theta+=_planar_dir.get_data(ix)*_radial_directions.get_data(ic,ix);
    }
    
    if(fabs(cos_theta)>1.0){
        if(fabs(cos_theta)>1.001){
            printf("WARNING somehow cos_theta>1.0 -- %e\n",cos_theta);
            exit(1);
        }
        sin_theta=0.0;
    }
    else{
        sin_theta=sqrt(1.0-cos_theta*cos_theta);
    }
    
    double x_shldbe,y_shldbe,rsq_shldbe,r_shldbe;
    rsq_shldbe=power(_radii.get_data(ic),2)/(cos_theta*cos_theta+sin_theta*sin_theta*16.0);
    r_shldbe=sqrt(rsq_shldbe);
    
    x_shldbe=r_shldbe*cos_theta;
    y_shldbe=r_shldbe*sin_theta;
    
    double x_is,y_is,d_radius;
    x_is=planar_radius*cos_theta;
    y_is=planar_radius*sin_theta;
    d_radius=sqrt(power(x_is-x_shldbe,2)+power(y_is-y_shldbe,2));
    
    
    out.set(0,(1.0-cos_theta)/_widths.get_data(ic,0));
    out.set(1,d_radius/_widths.get_data(ic,1));
    
    for(ix=0;ix<_dim;ix++){
        _dir.set(ix,pt.get_data(ix)-_centers.get_data(ic,ix));
    }
    
    for(ix=2;ix<_dim;ix++){
        out.set(ix,project_to_basis(ix,_dir)/_widths.get_data(ic,ix));
    }

    for(ix=0;ix<_dim;ix++){
        out.multiply_val(ix,0.01);
    }


}

//////////////////////ellipse classes////////////

ellipseData::ellipseData(int dd, int cc, int nData, double sigma) :
chiSquaredData(dd, cc, nData, sigma){

    _dir.set_name("ellipseData_dir");
    _projected.set_name("ellipseData_projected");
    
    int ix,ic;
    for(ic=0;ic<_ncenters;ic++){
        for(ix=0;ix<_dim;ix++){
            _centers.set(ic,ix,-30.0+60.0*_dice->doub());
        }
    }
}

void ellipseData::convert_params(array_1d<double> &pt, array_1d<double> &out, int ic){
    
    int ix;
    for(ix=0;ix<_dim;ix++){
        _dir.set(ix,pt.get_data(ix)-_centers.get_data(ic,ix));
    }
    
    project_to_basis(_dir,_projected);
    
    for(ix=0;ix<_dim;ix++){
        out.set(ix,_projected.get_data(ix)/_widths.get_data(ic,ix));
    }
}

nonGaussianEllipseData::nonGaussianEllipseData(int i1, int i2, int i3, double d1) : 
ellipseData(i1, i2, i3, d1){}

void nonGaussianEllipseData::convert_params(array_1d<double> &pt, array_1d<double> &out, int ic){

    int ix;
    for(ix=0;ix<_dim;ix++){
        _dir.set(ix,pt.get_data(ix)-_centers.get_data(ic,ix));
    }
    
    project_to_basis(_dir,_projected);
    
    
    for(ix=1;ix<_dim;ix++){
        out.set(ix,_projected.get_data(ix)/_widths.get_data(ic,ix));
    }
    
    double mu,xx;
    xx=_projected.get_data(0)/_widths.get_data(ic,0);
    if(fabs(xx)<fabs(xx-0.05)){
        out.set(0,xx);
    }
    else{
        out.set(0,xx-0.05);
    }
    
}
