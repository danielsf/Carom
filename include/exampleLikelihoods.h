#ifndef EXAMPLE_LIKELIHOODS_H
#define EXAMPLE_LIKELIHOODS_H

#include "jellyBean.h"

class integrableJellyBean : public jellyBeanData{

    public:
        integrableJellyBean() :
        jellyBeanData(4,1,1.0,100,0.4){
            _widths.set(0,0,100.0);
            _widths.set(0,1,8.0);
            _widths.set(0,2,8.0);
            _widths.set(0,3,8.0);
            _parabola_curvature=4.0;
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
