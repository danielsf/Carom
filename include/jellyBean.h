#ifndef JELLY_BEAN_H
#define JELLY_BEAN_H

#include "chisq.h"

class chiSquaredData : public chisquared{

    public:
        chiSquaredData(int, int, double, int, double);
        ~chiSquaredData(){};

        void set_width(int, int, double);

        virtual double operator()(array_1d<double>&);
        void write_data();
        void print_mins();
        
    protected:
        int _ndata;
        double _sig;
        array_1d<double> _x_values,_y_values,_sigma;
        array_1d<double> _mean_parameters,_aux_params;
        array_1d<double> _param_buffer;

        virtual void convert_params(array_1d<double>&, array_1d<double>&, int);

        double data_function(array_1d<double>&, double);
        
        void initialize_data();

};

class jellyBeanData : public chiSquaredData{

    public:
        ~jellyBeanData(){}
        jellyBeanData(int,int,double,int,double,double,double,double);
        

    protected:
        array_1d<double> _radii,_dir,_planar_dir;
        array_2d<double> _curvature_centers,_radial_directions;

        virtual void convert_params(array_1d<double>&, array_1d<double>&, int);
};

class ellipseData : public chiSquaredData{

    public:
        ~ellipseData(){}
        ellipseData(int, int, int, double);
    
    protected:
        array_1d<double> _dir,_projected;
        
        virtual void convert_params(array_1d<double>&, array_1d<double>&, int);

};

class nonGaussianEllipseData : public ellipseData{
    public:
        ~nonGaussianEllipseData(){}
        nonGaussianEllipseData(int, int, int, double);
    
    protected:
        virtual void convert_params(array_1d<double>&, array_1d<double>&, int);
};

#endif
