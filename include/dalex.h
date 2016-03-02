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
            _particle_log.set_name("dalex_particle_log");
            _target_factor=1.0;
            _simplex_mindex=-1;

            _basis_associates.set_name("dalex_basis_associates");
            _basis_mm.set_name("dalex_basis_mm");
            _basis_vv.set_name("dalex_basis_vv");
            _basis_lengths.set_name("dalex_basis_lengths");
            _basis_vectors.set_name("dalex_basis_vectors");
            _basis_ddsq.set_name("dalex_basis_ddsq");

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

        double get_norm(int);

    private:

        double chimin(){
            safety_check("chimin");
            return _chifn->chimin();
        }

        int mindex(){
            safety_check("mindex");
            return _chifn->mindex();
        }

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

        array_1d<int> _particles,_origins,_particle_log;
        chisq_wrapper *_chifn;
        double _target_factor;
        int _simplex_mindex;

        void _propagate_bisection(int);
        void _propagate_ricochet(int);
        void _propagate_midpt(int);

        ////////code related to finding basis vectors
        array_1d<int> _basis_associates;
        array_1d<double> _basis_mm,_basis_bb,_basis_model,_basis_vv;
        array_1d<double> _basis_lengths;
        array_2d<double> _basis_vectors,_basis_ddsq;

        double basis_error(array_2d<double>&, array_1d<double>&);

};

#endif
