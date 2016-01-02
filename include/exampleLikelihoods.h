#ifndef EXAMPLE_LIKELIHOODS_H
#define EXAMPLE_LIKELIHOODS_H

#include "jellyBean.h"

class gaussianJellyBean : public jellyBeanData{

    public:
        gaussianJellyBean() :
        jellyBeanData(4,1,1.0,100,0.4){
            _widths.set(0,0,250.0);
            _widths.set(0,1,1.0);
            _widths.set(0,2,1.0);
            _widths.set(0,3,3.0);

            _parabola_curvature=8.0;
        }

};


#endif
