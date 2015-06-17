#ifndef CONTROL_INTEGRATOR_H
#define CONTROL_INTEGRATOR_H

#include "goto_tools.h"
#include "containers.h"
#include "wrappers.h"
#include "kd.h"
#include "simplex.h"

class control_integrator{

    public:
        ~control_integrator(){}
        
        control_integrator(function_wrapper&,array_1d<double>&,array_1d<double>&,array_1d<double>&,char*);
        
        void set_chiLim3d(double);
        
        void run_analysis();
        void convert_to_boundary(array_2d<double>&,double,double,int,array_2d<double>&);
        int get_dex(double, double, double, int);
        void write_output(int, int,
                        array_1d<double>&,
                        array_1d<double>&,
                        array_1d<int>&, array_2d<double>&,
                        double,
                        array_1d<double>&, array_1d<double>&, array_1d<double>&,
                        double,
                        char*);
    
        void find_chi_min();
    
    private:
        array_1d<double> _min,_max,_dx;
        double _chiLim3d,_chi_min;
        function_wrapper *_chisq;
        char _name_root[letters];

};

#endif
