#ifndef DALEX_H
#define DALEX_H

#include "containers.h"
#include "goto_tools.h"
#include "chisq_wrapper.h"
#include "simplex.h"

class dalex{

    public:
        ~dalex(){};
        dalex(){
            _chifn=NULL;
            _particles.set_name("dalex_particles");
            _origins.set_name("dalex_origins");
            _target_factor=1.0;
            _simplex_mindex=-1;
        };

        void build(chisq_wrapper*);
        void search();
        void simplex_search();

        void calculate_gradient(int, array_1d<double>&);
        void propagate(int);

        int bisection(int, int, double, double);
        int bisection(int, array_1d<double>&, double, double);
        int bisection(array_1d<double>&, array_1d<double>&, double, double);

        void set_target_factor(double tt){
            _target_factor=tt;
        }

    private:

        void safety_check(char *word){
            if(_chifn==NULL){
                printf("ERROR: dalex called %s but _chifn is NULL\n", word);
                exit(1);
            }
        }

        double target(){
            safety_check("target");
            return _target_factor*_chifn->target();
        }

        array_1d<int> _particles,_origins;
        chisq_wrapper *_chifn;
        double _target_factor;
        int _simplex_mindex;

        void _propagate_bisection(int);
        void _propagate_ricochet(int);
        void _propagate_midpt(int);

};

#endif
