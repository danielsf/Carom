#ifndef DALEX_PLANCK_H
#define DALEX_PLANCK_H

#include "chisq.h"
#include "planck_wrapper.h"

class DalexPlanckLikelihood: public chisquared {

    public:

        DalexPlanckLikelihood(){
            _dim = 33;
            _called = 0;
            _time_spent = 0.0;
            _params = new double[_dim];
        }

        ~DalexPlanckLikelihood(){
            delete [] _params;
        }

        virtual double operator(const array_1d<double> &params){
            double before = double(time(NULL));
            int i;
            _called++;
            int is_valid = 1;
            for(i=0;i<_dim && is_valid==1;i++){
                if(params.get_data(i)<_mins.get_data(i)){
                    is_valid = 0;
                }
                if(params.get_data(i)>_maxs.get_data(i)){
                    is_valid = 0;
                }
            }

            if(is_valid == 0){
                _time_spent += double(time(NULL))-before;
                return 2.0*exception_value;
            }

            for(i=0;i<_dim;i++){
                _params[i] = params.get_data(i);
            }
            return _planck_wrapper(_params);
        }

    private:

        PlanckWrapper _planck_wrapper;
        int _dim;
        double *_params;

};

#endif
