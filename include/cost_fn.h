#ifndef COST_FN_H
#define COST_FN_H

#include "chisq_wrapper.h"

class cost_fn : public function_wrapper{
    public:
        cost_fn(chisq_wrapper*, array_1d<int>&);
        cost_fn(chisq_wrapper*, array_1d<int>&, int);
        cost_fn(){_chifn=NULL;};
        ~cost_fn(){};
        void build(chisq_wrapper*, array_1d<int>&, int);
        virtual double operator()(const array_1d<double>&);
        virtual int get_called();
        void multiply_norm(double dd){
            _scalar_norm*=dd;
        }
        double nn_distance(const array_1d<double>&);
        void set_envelope(double dd){
            _envelope=dd;
        }

        int get_cached_values(const int dex, double *fn){

            int i;
            for(i=0;i<_pt_cache.get_dim();i++){
                if(_pt_cache.get_data(i)==dex){
                    fn[0]=_fn_cache.get_data(i);
                    return 1;
                }
            }
            return 0;

        }

        void set_bases(const array_2d<double> &b_in){
            int i,j;
            _bases.reset();
            _bases.set_dim(b_in.get_rows(),b_in.get_cols());
            for(i=0;i<b_in.get_rows();i++){
                for(j=0;j<b_in.get_cols();j++){
                    _bases.set(i,j,b_in.get_data(i,j));
                }
            }
            _set_d_params();
        }

        double get_norm(int i){
            return _relative_norm.get_data(i);
        }

        double basis(int i, int j){
            return _bases.get_data(i,j);
        }

    private:
        array_1d<int> _associates;
        double _scalar_norm;
        chisq_wrapper *_chifn;
        int _called;
        double _envelope;
        int _min_or_med;

        array_1d<int> _pt_cache;
        array_1d<double> _fn_cache;
        array_1d<double> _relative_norm;

        array_2d<double> _bases;
        array_2d<double> _projected_associates;
        array_1d<double> _pt_projected;

        void _project_associates();
        void _set_d_params();
        void _set_scalar_norm();

};

#endif
