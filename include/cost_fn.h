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

        void set_mask(array_1d<int> &mm){
            int i;
            if(mm.get_dim()!=_associates.get_dim()){
                printf("WARNING trying to set %d masks for %d associates\n",
                mm.get_dim(),_associates.get_dim());
                printf("in cost_fn\n");
                exit(1);
            }
            for(i=0;i<mm.get_dim();i++){
                _mask.set(i,mm.get_data(i));
            }

            array_1d<double> max,min;
            int j;
            for(i=0;i<_associates.get_dim();i++){
                if(_mask.get_data(i)==1){
                    for(j=0;j<_chifn->get_dim();j++){
                        if(j>=max.get_dim() || _chifn->get_pt(_associates.get_data(i),j)>max.get_data(j)){
                            max.set(j,_chifn->get_pt(_associates.get_data(i),j));
                        }
                        if(j>=min.get_dim() || _chifn->get_pt(_associates.get_data(i),j)<min.get_data(j)){
                            min.set(j,_chifn->get_pt(_associates.get_data(i),j));
                        }
                    }
                }
            }

            for(i=0;i<_chifn->get_dim();i++){
                if(max.get_data(i)-min.get_data(i)>1.0e-20 && max.get_data(i)-min.get_data(i)<_scalar_norm){
                    _scalar_norm=max.get_data(i)-min.get_data(i);
                }
            }

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

        void set_bases(){
            if(_associates.get_dim()<_chifn->get_dim()){
                _random_set_bases();
            }
            else{
                _principal_set_bases();
            }
        }


        double get_hyper_norm(int ii){
            return _hyper_norm.get_data(ii);
        }

    private:
        array_1d<int> _associates;
        array_1d<int> _mask;
        array_1d<double> _median_associate;
        double _scalar_norm;
        chisq_wrapper *_chifn;
        int _called;
        int _just_median;
        double _envelope;

        array_2d<double> _bases;

        void _principal_set_bases();
        void _random_set_bases();

        void _set_hyper_ellipse();
        double _hyper_ellipse_distance(array_1d<double>&);
        array_1d<double> _hyper_center;
        array_1d<double> _hyper_norm;

};

#endif
