#ifndef EXAMPLE_LIKELIHOODS_H
#define EXAMPLE_LIKELIHOODS_H

#include "jellyBean.h"

class integrableJellyBean : public jellyBeanData{

    public:
        integrableJellyBean() :
        jellyBeanData(4,1,1.0,100,0.4){
            _widths.set(0,0,100.0);
            _widths.set(0,1,2.0);
            _widths.set(0,2,10.0);
            _widths.set(0,3,10.0);

            _parabola_curvature=4.0;

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
