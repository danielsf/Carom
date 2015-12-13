#ifndef EXAMPLE_LIKELIHOODS_H
#define EXAMPLE_LIKELIHOODS_H

#include "jellyBean.h"

class gaussianJellyBean : public jellyBeanData{

    public:
        gaussianJellyBean() :
        jellyBeanData(4,1,1.0,100,0.4,0.4,0.02,20.0){
            _widths.set(0,2,1.0);
            _widths.set(0,3,5.0);
        }

};


#endif
