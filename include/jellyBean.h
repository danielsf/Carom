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


class chiSquaredData : public chisquared{

    public:
        chiSquaredData(int, int, int, double);
        ~chiSquaredData(){};

        virtual double operator()(array_1d<double>&);
        void write_data();
        void print_mins();
        
    protected:
        int _ndata;
        double _sig;
        array_1d<double> _x_values,_y_values,_sigma;
        array_1d<double> _mean_parameters,_aux_params;
        array_1d<double> _param_buffer,_projected_pt,_dir;

        virtual void convert_params(array_1d<double>&, array_1d<double>&, int);

        double data_function(array_1d<double>&, double);
        
        void initialize_data();

};

class jellyBeanData : public chiSquaredData{

    public:
        ~jellyBeanData(){}
        jellyBeanData(int,int,int,double,double,double,double);
        

    protected:
        int _ndata;
        double _sig;
        array_1d<double> _radii;
        array_2d<double> _curvature_centers,_radial_directions;

        void convert_params(array_1d<double>&, array_1d<double>&, int);
};


#endif
