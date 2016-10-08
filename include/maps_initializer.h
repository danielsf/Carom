#ifndef MAPS_INIT_H
#define MAPS_INIT_H

#include "chisq_wrapper.h"
#include "simplex.h"

class maps_initializer{

    public:

        ~maps_initializer(){}

        maps_initializer(){
            _chifn=NULL;
            _transform.set_name("maps_init_transform");
            _bases.set_name("maps_init_bases");
            _particles.set_name("maps_init_particles");
            _radii.set_name("maps_init_radii");
            _center.set_name("maps_init_center");
            _fn.set_name("maps_init_fn");
            _dex.set_name("maps_init_dex");
            _dir.set_name("maps_init_dir");
            _pt_shell.set_name("maps_init_pt_shell");
            _pt_true.set_name("maps_init_pt_true");

        }

        void set_chifn(chisq_wrapper *cc){
            _chifn=cc;
        }

        void safety_check(){
            if(_chifn==NULL){
                printf("WARNING _chifn is null\n");
                exit(1);
            }
        }

        void initialize();
        void search();
        void set_bases();
        void convert_to_truth(const array_1d<double>&,array_1d<double>&);
        void convert_to_shell(const array_1d<double>&,array_1d<double>&);
        void evaluate(array_1d<double>&, double*, int*);
        void sample();

    private:
        chisq_wrapper *_chifn;
        array_1d<double> _transform;
        array_2d<double> _particles;
        array_1d<double> _fn;
        array_1d<int> _dex;
        array_2d<double> _bases;
        array_1d<double> _radii;
        array_1d<double> _center;
        array_1d<double>_pt_true;
        array_1d<double>_pt_shell;
        array_1d<double> _dir;
        double _fn_max;
        int _max_dex;
        int _ct_replace;
        double _vol;
};


#endif
