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
        chisq_wrapper *_chisq;
        int _called;
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

#endif
