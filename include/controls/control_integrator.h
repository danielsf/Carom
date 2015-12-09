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
        
        void set_chi_lim_freq(double);
        
        void run_analysis(double);
        void run_analysis();

        int get_dex(double, double, double, int);
        void write_output(int, int,
                        array_1d<double>&,
                        array_1d<double>&,
                        array_1d<int>&, array_2d<int>&,
                        double,
                        double,
                        char*);
    
        void find_chi_min(array_1d<double>&, array_1d<double>&);
    
    private:
        array_1d<double> _min,_max,_dx;
        double _chi_lim_freq,_delta_chi_bayes,_chi_min;
        double _confidence_limit;
        function_wrapper *_chisq;
        char _name_root[letters];

};

#endif
