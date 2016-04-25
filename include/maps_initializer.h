#ifndef MAPS_INIT_H
#define MAPS_INIT_H

#include "chisq_wrapper.h"
#include "simplex.h"

class maps_initializer{

    public:
    
        ~maps_initializer(){}
    
        maps_initializer(){
            _chifn=NULL;
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
};


#endif
