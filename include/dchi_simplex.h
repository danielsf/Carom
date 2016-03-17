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

    private:
        array_1d<int> _associates;
        double _norm;
        chisq_wrapper *_chifn;
        int _called;
        double _envelope;
};

#endif
