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
        double _xmax,_dx;
        array_1d<double> _x_values,_y_values,_sigma;
        array_1d<double> _wave_phase,_wave_lambda,_env_d_center,_env_width;
        array_1d<double> _wave_amp;
        array_1d<double> _param_buffer;

        virtual void convert_params(array_1d<double>&, array_1d<double>&, int);

        double data_function(array_1d<double>&, double);
        
        void initialize_data();

};

class jellyBeanData : public chiSquaredData{

    public:
        ~jellyBeanData(){}
        jellyBeanData(int,int,double,int,double);
        

    protected:
        double _parabola_curvature;
        array_2d<double> _parabola_centers,_parabola_x,_parabola_y;

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
