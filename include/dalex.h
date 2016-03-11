#ifndef DALEX_H
#define DALEX_H

#include "containers.h"
#include "goto_tools.h"
#include "chisq_wrapper.h"
#include "simplex.h"
#include "eigen_wrapper.h"
#include "search_types.h"
#include "gp_lin.h"
#include "dchi_simplex_gp.h"

class dalex{

    public:
        ~dalex(){};
        dalex(){
            _chifn=NULL;
            _good_points.set_name("dalex_good_points");
            _explorers.set_name("dalex_explorers");
            _explorer_temp=1.0;
            _explorer_step=1.0;
            _target_factor=1.0;
            _simplex_mindex=-1;
            _last_checked_good=0;
            _log=NULL;

            _charges.set_name("dalex_charges");

            _basis_chimin=2.0*exception_value;
            _basis_associates.set_name("dalex_basis_associates");
            _basis_mm.set_name("dalex_basis_mm");
            _basis_vv.set_name("dalex_basis_vv");
            _basis_lengths.set_name("dalex_basis_lengths");
            _basis_vectors.set_name("dalex_basis_vectors");
            _basis_ddsq.set_name("dalex_basis_ddsq");

        };

        void build(chisq_wrapper*);

        void set_log(asymm_array_2d<int> *_ll){
            _log=_ll;
        }

        void search();
        void simplex_search();
        void simplex_search(int);
        void simplex_boundary_search();
        void simplex_boundary_search(int);
        void explore();

        int bisection(int, int, double, double);
        int bisection(int, array_1d<double>&, double, double);
        int bisection(array_1d<double>&, array_1d<double>&, double, double);

        double get_basis(int i, int j){
            return _basis_vectors.get_data(i,j);
        }

        void set_target_factor(double tt){
            _target_factor=tt;
        }

        double get_norm(int);

    private:

        void evaluate(array_1d<double> &pt, double *mu_out, int *i_out){
            _chifn->evaluate(pt,mu_out,i_out);
            if(mu_out[0]<target() && _good_points.contains(i_out[0])==0){
                _good_points.add(i_out[0]);
            }
        }

        void _add_good_points(){
            _add_good_points(_last_checked_good);
        }

        void _add_good_points(int i_start){
            safety_check("_add_good_points");
            int i;
            for(i=i_start;i<_chifn->get_pts();i++){
                if(_chifn->get_fn(i)<target() && _good_points.contains(i)==0){
                    _good_points.add(i);
                }
            }
            _last_checked_good=_chifn->get_pts();
        }

        void add_external_good_points(){
            safety_check("add_good_points");
            int i,j;
            array_1d<double> trial;
            trial.set_name("dalex_add_good_trial");
            double mu;
            int i_found;
            for(i=_last_checked_good;i<_chifn->get_pts();i++){
                if(_chifn->get_fn(i)<target() && _good_points.contains(i)==0){
                    for(j=0;j<_chifn->get_dim();j++){
                        trial.set(j,0.5*(_chifn->get_pt(_chifn->mindex(),j)+
                                          _chifn->get_pt(i,j)));
                    }

                    evaluate(trial,&mu,&i_found);
                    if(mu<target() && _good_points.contains(i)==0){
                        _good_points.add(i);
                    }
                }
            }

            _last_checked_good=_chifn->get_pts();
        }

        void assess_good_points(){
            int i;
            for(i=0;i<_good_points.get_dim();i++){
                if(_chifn->get_fn(_good_points.get_data(i))>target()){
                    _good_points.remove(i);
                    i--;
                }
            }
        }

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

        asymm_array_2d<int> *_log;
        chisq_wrapper *_chifn;
        double _target_factor;
        int _simplex_mindex;
        array_1d<int> _good_points;

        ////////code related to finding basis vectors
        array_1d<int> _basis_associates;
        array_1d<double> _basis_mm,_basis_bb,_basis_model,_basis_vv;
        array_1d<double> _basis_lengths;
        array_2d<double> _basis_vectors,_basis_ddsq;
        double _basis_chimin;

        double basis_error(array_2d<double>&, array_1d<double>&);
        void find_trial_bases(int, array_1d<double>&, array_2d<double> &out_bases);
        void validate_bases(array_2d<double>&, char*);
        void guess_bases(array_2d<double>&);
        void find_covariance_matrix(int, array_2d<double>&);
        void find_bases();


        ///////code related to explorers
        array_1d<int> _explorers;
        double _explorer_temp,_explorer_step;
        int _last_checked_good;


        //////code related to tendrils
        void get_gradient(int,array_1d<double>&,array_1d<double>&);
        void tendril_search();
        array_1d<int> _charges;

        void assess_charges(){
            int i;
            for(i=0;i<_charges.get_dim();i++){
                if(_chifn->get_fn(_charges.get_data(i))>target()){
                    _charges.remove(i);
                    i--;
                }
            }
        }

        void add_charge(int ii){
            if(_charges.contains(ii)==0){
                _charges.add(ii);
            }
        }

};

#endif
