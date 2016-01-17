#ifndef EXAMPLE_LIKELIHOODS_H
#define EXAMPLE_LIKELIHOODS_H

#include "jellyBean.h"

class integrableJellyBean : public jellyBeanData{

    public:
        integrableJellyBean() :
        jellyBeanData(4,1,1.0,100,0.4){
            _widths.set(0,0,100.0);
            _widths.set(0,1,5.0);
            _widths.set(0,2,5.0);
            _widths.set(0,3,40.0);
            _parabola_curvature=4.0;
        }

    private:

        virtual void convert_params(array_1d<double> &pt_in, array_1d<double> &out, int ic){

            array_1d<double> pt;
            pt.set_name("convert_params_pt");
            project_to_basis(pt_in,pt);

            double x_is,y_is;
            x_is=(pt.get_data(0)-_parabola_centers.get_data(ic,0))*_parabola_x.get_data(ic,0)
                 +(pt.get_data(1)-_parabola_centers.get_data(ic,1))*_parabola_x.get_data(ic,1);

            y_is=(pt.get_data(0)-_parabola_centers.get_data(ic,0))*_parabola_y.get_data(ic,0)
                 +(pt.get_data(1)-_parabola_centers.get_data(ic,1))*_parabola_y.get_data(ic,1);

            double rr=sqrt(x_is*x_is+y_is*y_is);

            double cos_theta,sin_theta;
            cos_theta=x_is/rr;
            sin_theta=y_is/rr;

            if(fabs(cos_theta)>1.0){
                if(fabs(cos_theta)>1.001){
                    printf("WARNING somehow cos_theta>1.0 -- %e\n",cos_theta);
                    exit(1);
                }
            }

            if(fabs(sin_theta)>1.0){
                if(fabs(sin_theta)>1.001){
                    printf("WARNING somehow sin_theta>1.0 -- %e\n",sin_theta);
                    exit(1);
                }
            }

            double r_shldbe;

            if(fabs(cos_theta)<1.0e-5 && sin_theta>0.0){
                r_shldbe=1.0e10/_parabola_curvature;
            }
            else if(fabs(cos_theta)<1.0e-10 && sin_theta<0.0){
                r_shldbe=0.25/_parabola_curvature;
            }
            else{
                r_shldbe=(1.0+sin_theta)/(2.0*_parabola_curvature*cos_theta*cos_theta);
            }

            double d_radius;
            d_radius=fabs(rr-r_shldbe);

            double y_distance;
            y_distance=(y_is+0.25/_parabola_curvature)/_widths.get_data(ic,0);
            out.set(0,y_distance);

            double x_shldbe,dx;
            if(y_distance<0.0){
                out.set(1,d_radius/_widths.get_data(ic,1));
            }
            else{
                x_shldbe=sqrt((y_is+0.25/_parabola_curvature)/_parabola_curvature);
                if(fabs(x_is+x_shldbe)<fabs(x_is-x_shldbe)){
                    dx=fabs(x_is+x_shldbe);
                }
                else{
                    dx=fabs(x_is-x_shldbe);
                }

                out.set(1,dx/_widths.get_data(ic,1));
            }

            if(out.get_data(1)<0.0){
                printf("WARNING out 1 %e\n",out.get_data(1));
                exit(1);
            }

            int ix;
            for(ix=2;ix<_dim;ix++){
                if(ix%4==2){
                    out.set(ix,exp((pt.get_data(ix)-_centers.get_data(ic,ix))/_widths.get_data(ic,ix)));
                }
                else  if(ix%4==3){
                    if(pt.get_data(ix)>_centers.get_data(ic,ix)){
                        out.set(ix,log(1.0+(pt.get_data(ix)-_centers.get_data(ic,ix))/_widths.get_data(ic,ix)));
                    }
                    else{
                        out.set(ix,0.3*(pt.get_data(ix)-_centers.get_data(ic,ix))/_widths.get_data(ic,ix));
                    }
                }
                else{
                    out.set(ix,(pt.get_data(ix)-_centers.get_data(ic,ix))/_widths.get_data(ic,ix));
                }
            }

            for(ix=0;ix<_dim;ix++){
                out.multiply_val(ix,0.01);
            }


        }

};

class gaussianJellyBean4 : public jellyBeanData{

    public:
        gaussianJellyBean4() :
        jellyBeanData(4,1,1.0,100,0.4){
            _widths.set(0,0,250.0);
            _widths.set(0,1,2.0);
            _widths.set(0,2,1.0);
            _widths.set(0,3,3.0);

            _parabola_curvature=8.0;
        }

};

class gaussianJellyBean12 : public jellyBeanData{

    public:
        gaussianJellyBean12() :
        jellyBeanData(12,1,1.0,100,0.4){
            Ran constructor_dice(44);
            int i;

            _widths.set(0,0,250.0);
            _widths.set(0,1,2.0);
            _parabola_curvature=8.0;

            for(i=2;i<12;i++){
                _widths.set(0,i,0.5+3.0*constructor_dice.doub());
            }
        }
};


class gaussianJellyBean24 : public jellyBeanData{

    public:
        gaussianJellyBean24() :
        jellyBeanData(24,1,1.0,100,0.4){
            Ran constructor_dice(99);
            int i;

            _widths.set(0,0,250.0);
            _widths.set(0,1,2.0);
            _parabola_curvature=8.0;

            for(i=2;i<12;i++){
                _widths.set(0,i,0.5+3.0*constructor_dice.doub());
            }
        }

};

#endif
