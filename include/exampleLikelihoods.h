#ifndef EXAMPLE_LIKELIHOODS_H
#define EXAMPLE_LIKELIHOODS_H

#include "jellyBean.h"

class integrableJellyBean : public jellyBeanData{

    public:
        integrableJellyBean() :
        jellyBeanData(4,1,1.0,100,0.4){
            _widths.set(0,0,10.0);
            _widths.set(0,1,2.0);
            _widths.set(0,2,10.0);
            _widths.set(0,3,10.0);

            _parabola_curvature=0.5;

            _bases.set(0,0,0.3);
            _bases.set(0,1,sqrt(1.0-0.3*0.3));
            _bases.set(0,2,0.0);
            _bases.set(0,3,0.0);

            _bases.set(1,0,-1.0*_bases.get_data(0,1));
            _bases.set(1,1,_bases.get_data(0,0));
            _bases.set(1,2,0.0);
            _bases.set(1,3,0.0);

            _bases.set(2,0,0.0);
            _bases.set(2,1,0.0);
            _bases.set(2,2,0.4);
            _bases.set(2,3,sqrt(1.0-0.4*0.4));

            _bases.set(3,0,0.0);
            _bases.set(3,1,0.0);
            _bases.set(3,2,_bases.get_data(2,3));
            _bases.set(3,3,-1.0*_bases.get_data(2,2));
        }
};

class integrableJellyBean12 : public jellyBeanData{

    public:
        integrableJellyBean12() : jellyBeanData(12,1,1.0,100,0.4){
            Ran local_dice(1132);
            int i;
            for(i=0;i<8;i++){
                _d8_center.set(i, -10.0+20.0*local_dice.doub());
                _d8_widths.set(i,local_dice.doub()*2.0);
            }
            _d8_bases.set_cols(8);
            array_1d<double> trial;
            int j;
            double component;
            while(_d8_bases.get_rows()<8){
                for(i=0;i<8;i++){
                    trial.set(i,normal_deviate(&local_dice, 0.0, 1.0));
                }
                for(i=0;i<_d8_bases.get_rows();i++){
                    component=0.0;
                    for(j=0;j<8;j++){
                        component+=trial.get_data(j)*_d8_bases.get_data(i,j);
                    }
                    for(j=0;j<8;j++){
                        trial.subtract_val(j,component*_d8_bases.get_data(i,j));
                    }
                }
                component=trial.normalize();
                if(component>1.0e-20){
                    _d8_bases.add_row(trial);
                }
            }
        }

        virtual double operator()(array_1d<double> &pt){
            array_1d<double> four_d_pt;
            int i;
            for(i=0;i<4;i++){
                four_d_pt.set(i,pt.get_data(i));
            }
            double d4_val = _d4_chisq(four_d_pt);
            double d8_val = 0.0;
            double component;
            int ix;
            for(ix=0;ix<8;ix++){
                component=0.0;
                for(i=0;i<8;i++){
                    component+=(pt.get_data(4+i)-_d8_center.get_data(i))*_d8_bases.get_data(ix,i);
                }
                d8_val += power(component/_d8_widths.get_data(ix),2);
            }

            return d4_val + d8_val;
        }

    private:
        integrableJellyBean _d4_chisq;
        array_1d<double> _d8_center;
        array_2d<double> _d8_bases;
        array_1d<double> _d8_widths;

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

            for(i=2;i<24;i++){
                _widths.set(0,i,0.5+3.0*constructor_dice.doub());
            }
        }

};

#endif
