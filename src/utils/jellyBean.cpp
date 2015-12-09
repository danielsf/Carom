#include "jellyBean.h"

////////////////////////////////jellyBeanData

chiSquaredData::chiSquaredData(int dd, int cc, double ww, int nData, double sigma) : chisquared(dd, cc, ww){
    printf("_dim %d\n",_dim);

    if(_dim<4){
        printf("WARNING must have at least 4 params in chisquaredData\n");
        exit(1);
    }

    if(_dim%4!=0){
        printf("WARNING dimensionality must be multiple of four\n");
        exit(1);
    }

    _x_values.set_name("jellyBeanData_x");
    _y_values.set_name("jellyBeanData_y");
    _sigma.set_name("jellyBeanData_sigma");
    _wave_phase.set_name("jellyBeanData_phase");
    _wave_lambda.set_name("jellyBeanData_lambda");
    _wave_amp.set_name("jellBeanData_amp");
    _env_center.set_name("jellyBeanData_env_center");
    _env_width.set_name("jellyBeanData_env_width");
    _param_buffer.set_name("jellyBeanData_param_buffer");

    Ran local_dice(nData+dd+cc);

    _ndata=nData;
    _sig=sigma;
    _xmax=3.0;
    _dx=0.03;

    int ix;
    double ll;
    double ll_min,ll_max;
    ll_min=log(0.05*_xmax);
    ll_max=log(0.6*_xmax);
    for(ix=0;ix<_dim;ix+=4){
        _wave_phase.add(local_dice.doub()*_xmax);
        ll=local_dice.doub()*(ll_max-ll_min)+ll_min;
        _wave_lambda.add(exp(ll));
        _wave_amp.add(local_dice.doub()*5.0);
        _env_center.add(local_dice.doub()*_xmax);
        _env_width.add((local_dice.doub()*0.5+0.05)*_xmax);
    }

}

void chiSquaredData::set_width(int ic, int id, double dd){
    _widths.set(ic, id, dd);
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
            printf("    %e %e -- width %e\n",
            _centers.get_data(ic,ix),
            params.get_data(ix),_widths.get_data(0,ix));
        }
    }

    printf("\nbases\n");
    for(ix=0;ix<_dim;ix++){
        for(ic=0;ic<_dim;ic++){
            printf("%.3e ",_bases.get_data(ic,ix));
        }
        printf("\n");
    }

    printf("\nparams\n");
    for(ix=0;ix<_dim/4;ix++){
        printf("amp %e\n",_wave_amp.get_data(ix));
        printf("phase %e\n",_wave_phase.get_data(ix));
        printf("lambda %e\n",_wave_lambda.get_data(ix));
        printf("env_x %e\n",_env_center.get_data(ix));
        printf("env_width %e\n",_env_width.get_data(ix));
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
    for(ix=0, xx=0.0;xx<_xmax;xx+=_dx, ix++){
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
    print_mins();

}

void chiSquaredData::convert_params(array_1d<double> &pt, array_1d<double> &out, int ic){
    printf("Called void convert_params");
    exit(1);
}

double chiSquaredData::data_function(array_1d<double> &params, double xx){

    double ans=0.0;

    double envelope;
    double env_x;
    double wave;
    double amp;
    double lambda;
    double phase;
    int ix,i_param;
    for(ix=0,i_param=0;ix<_dim;ix+=4,i_param++){
        env_x=(xx-params.get_data(ix)-_env_center.get_data(i_param))/(params.get_data(ix+1)+_env_width.get_data(i_param));
        envelope=0.1+exp(-0.5*power(env_x,2));
        amp=params.get_data(ix+2)+_wave_amp.get_data(i_param);
        lambda=params.get_data(ix+3)*0.5*_xmax+_wave_lambda.get_data(i_param);
        wave=sin(2.0*pi*(xx-_wave_phase.get_data(i_param))/lambda);
        ans+=envelope*amp*wave;
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

    if(_with_logging==1){
        _log_point(pt, chisq_min);
    }

    return chisq_min;
}


jellyBeanData::jellyBeanData(int dd, int cc, double wc, int nData, double sigma, double ww, double radial, double rr) :
chiSquaredData(dd, cc, wc, nData, sigma){

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
            _curvature_centers.set(ic,ix,_centers.get_data(ic,ix));
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
chiSquaredData(dd, cc, 1.0, nData, sigma){

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
