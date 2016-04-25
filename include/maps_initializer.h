#ifndef MAPS_INIT_H
#define MAPS_INIT_H

#include "chisq_wrapper.h"
#include "simplex.h"

class maps_initializer{

    public:
    
        ~maps_initializer(){}
    
        maps_initializer(){
            _chifn=NULL;
            _particles.set_name("maps_init_particles");
            _local_min.set_name("maps_init_local_min");
            _abs_min.set_name("maps_init_abs_min");
            _since_min.set_name("maps_init_since_min");
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

        double evaluate(array_1d<double> &pt, int *i_found){
            double mu;
            _chifn->evaluate(pt, &mu, i_found);
            return mu;
        }


        void search();

    private:
        chisq_wrapper *_chifn;
        array_1d<int> _particles,_local_min,_abs_min,_since_min;
};


#endif
