#ifndef COST_FN_H
#define COST_FN_H

#include "chisq_wrapper.h"

class cost_fn : public function_wrapper{
    public:
        cost_fn(chisq_wrapper*, array_1d<int>&);
        ~cost_fn(){};
        virtual double operator()(const array_1d<double>&);
        virtual int get_called();
        double nn_distance(const array_1d<double>&);
        void set_envelope(double dd){
            _envelope=dd;
        }

        void use_median(){
            _just_median=1;
        }


        void copy_bases(array_2d<double> &out){
            int i;
            out.reset_preserving_room();
            for(i=0;i<_bases.get_rows();i++){
                out.add_row(_bases(i));
            }
        }

        double get_norm(int ii){
            return _norm.get_data(ii);
        }

    private:
        array_1d<int> _associates;
        array_1d<double> _median_associate;
        chisq_wrapper *_chifn;
        int _called;
        int _just_median;
        double _envelope;

        array_2d<double> _bases;

        void _principal_set_bases();
        void _random_set_bases();

        void _set_bases(){
            if(_associates.get_dim()<_chifn->get_dim()){
                _random_set_bases();
            }
            else{
                _principal_set_bases();
            }
        }


        void _set_norm();
        array_1d<double> _norm;
        array_1d<double> _cardinal_norm;

};

#endif
