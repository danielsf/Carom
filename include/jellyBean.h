#ifndef JELLY_BEAN_H
#define JELLY_BEAN_H

#include "chisq.h"

class jellyBean : public chisquared{

    public:
        ~jellyBean();
        jellyBean(int,double, double);
        virtual double operator()(array_1d<double>&);
        void get_curvature_center(array_1d<double>&);

    private:
        array_1d<double> _curvature_center,_radial_direction;
        double _curvature_radius;

};

#endif
