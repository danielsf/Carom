#ifndef DCHI_SIMPLEX_H
#define DCHI_SIMPLEX_H

#include "chisq_wrapper.h"

class dchi_simplex_base : public function_wrapper{

    public:
        dchi_simplex_base(chisq_wrapper*, array_1d<int>&);
        ~dchi_simplex_base(){};

        virtual int get_called();
        double associate_distance(array_1d<double>&);
        virtual double operator()(array_1d<double>&);

    protected:
        array_1d<int> _associates;
        array_1d<double> _norm;
        chisq_wrapper *_chifn;
        int _called;
        double _min_0;
};

class dchi_boundary_simplex : public dchi_simplex_base{

    public:
        dchi_boundary_simplex(chisq_wrapper*, array_1d<int>&);
        ~dchi_boundary_simplex(){};

        virtual double operator()(array_1d<double>&);
};


class dchi_multimodal_simplex : public dchi_simplex_base{

    public:
        dchi_multimodal_simplex(chisq_wrapper*, array_1d<int>&);
        ~dchi_multimodal_simplex(){};

        void set_norm(array_1d<double>&);
        virtual double operator()(array_1d<double>&);
};

class dchi_interior_simplex : public function_wrapper{
    public:
        dchi_interior_simplex(chisq_wrapper*, array_1d<int>&);
        ~dchi_interior_simplex(){};
        virtual double operator()(array_1d<double>&);
        virtual int get_called();
        double nn_distance(array_1d<double>&);
        void set_envelope(double dd){
            _envelope=dd;
        }

        void set_mask(array_1d<int> &mm){
            int i;
            if(mm.get_dim()!=_associates.get_dim()){
                printf("WARNING trying to set %d masks for %d associates\n",
                mm.get_dim(),_associates.get_dim());
                printf("in dchi_interior_simplex\n");
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
                if(max.get_data(i)-min.get_data(i)>1.0e-20 && max.get_data(i)-min.get_data(i)<_norm){
                    _norm=max.get_data(i)-min.get_data(i);
                }
            }

        }

    private:
        array_1d<int> _associates;
        array_1d<int> _mask;
        double _norm;
        chisq_wrapper *_chifn;
        int _called;
        double _envelope;
};

#endif
