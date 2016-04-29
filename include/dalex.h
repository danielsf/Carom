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
            _end_points.set_name("dalex_end_points");
            _explorers.set_name("dalex_explorers");
            _tendril_walkers.set_name("dalex_tendril_walkers");
            _explorer_temp=1.0;
            _explorer_step=1.0;
            _target_factor=1.0;
            _simplex_mindex=-1;
            _tendril_origin=-1;
            _last_checked_good=0;
            _log=NULL;

            _charges.set_name("dalex_charges");

            _basis_chimin=2.0*exception_value;
            _basis_associates.set_name("dalex_basis_associates");
            _basis_mm.set_name("dalex_basis_mm");
            _basis_vv.set_name("dalex_basis_vv");
            _basis_norm.set_name("dalex_basis_norm");
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
        void simplex_search(array_1d<int>&);
        void simplex_gp_search();
        void simplex_boundary_search();
        void simplex_boundary_search(int,array_1d<double>&);
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

       double distance(int i1, int i2){
           int i,j;
           double component,dd;
           dd=0.0;
           for(i=0;i<_chifn->get_dim();i++){
               component=0.0;
               for(j=0;j<_chifn->get_dim();j++){
                   component+=(_chifn->get_pt(i1,j)-_chifn->get_pt(i2,j))*_basis_vectors.get_data(i,j);
               }
               dd+=power(component/_basis_norm.get_data(i),2);
           }
           return sqrt(dd);
       }

       void add_good_point(int ii, int aa){
           if(_chifn->get_fn(ii)<target() && _good_points.contains(ii)==0){
               _good_points.add(ii);
               _good_point_origins.add(aa);
           }
       }

        void evaluate(array_1d<double> &pt, double *mu_out, int *i_out){
            evaluate(pt, mu_out, i_out, -1);
        }

        void evaluate(array_1d<double> &pt, double *mu_out, int *i_out, int i_origin){
            _chifn->evaluate(pt,mu_out,i_out);
            if(mu_out[0]<target() && _good_points.contains(i_out[0])==0){
                add_good_point(i_out[0], i_origin);
            }
        }

        void _update_good_points(){
            _update_good_points(_last_checked_good, -1, -1);
        }

        void _update_good_points(int ii){
            _update_good_points(ii, -1, -1);
        }

        void _update_good_points(int i_start, int io1, int io2){
            safety_check("_update_good_points");

            int i;
            double dd1,dd2;
            int origin_dex;
            for(i=i_start;i<_chifn->get_pts();i++){
                if(_chifn->get_fn(i)<target() && _good_points.contains(i)==0){
                    if(io1<0 && io2>=0){
                        origin_dex=io2;
                    }
                    else if(io2<0 && io1>=0){
                         origin_dex=io1;
                    }
                    else if(io1>=0 && io2>=0){
                        dd1=distance(i,io1);
                        dd2=distance(i,io2);
                        if(dd1<dd2){
                            origin_dex=io1;
                        }
                        else{
                            origin_dex=io2;
                        }
                    }
                    else{
                        origin_dex=-1;
                    }
                    add_good_point(i, origin_dex);
                }
            }
            _last_checked_good=_chifn->get_pts();
        }

        void update_good_points_external(){
            safety_check("update_good_points_external");
            int i,j;
            array_1d<double> trial;
            trial.set_name("dalex_update_good_external_trial");
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
                        add_good_point(i, -1);
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
                    _good_point_origins.remove(i);
                    i--;
                }
            }
            assess_good_point_origins();
        }

        void assess_good_point_origins();

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

            if(_good_points.get_dim()!=_good_point_origins.get_dim()){
                printf("ERROR: good points %d good point origins %d\n",
                _good_points.get_dim(),_good_point_origins.get_dim());
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
        array_1d<int> _good_point_origins;

        ////////code related to finding basis vectors
        array_1d<int> _basis_associates;
        array_1d<double> _basis_mm,_basis_bb,_basis_model,_basis_vv;
        array_1d<double> _basis_norm;
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
        void tendril_seed(function_wrapper*, int, array_2d<double>&);
        array_1d<int> _charges;
        array_1d<int> _end_points;
        array_2d<double> _tendril_walkers;
        int _tendril_origin;

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
            if(_charges.contains(ii)==0 && _chifn->get_fn(ii)<target()){
                _charges.add(ii);
            }
        }


        int check_association(int i1, int i2){
            array_1d<double> trial;
            trial.set_name("dalex_check_association_trial");
            double wgt,mu;
            int i,i_found,i_origin;
            for(wgt=0.25;wgt<0.77;wgt+=0.25){

                if(wgt<0.5){
                    i_origin=i1;
                }
                else{
                    i_origin=i2;
                }

                for(i=0;i<_chifn->get_dim();i++){
                    trial.set(i,wgt*_chifn->get_pt(i1,i)+(1.0-wgt)*_chifn->get_pt(i2,i));
                }

                evaluate(trial,&mu,&i_found,i_origin);

                if(mu>target()){
                    return 0;
                }
            }
            return 1;
        }


        void create_mask(int i_pt, array_1d<int> &mask){
            mask.reset_preserving_room();

            array_1d<int> considered,is_associate;
            considered.set_name("dalex_create_mask_considered");
            is_associate.set_name("dalex_create_mask_is_associate");
            int i,j,i_origin,n_0;
            n_0=_good_points.get_dim();
            for(i=0;i<n_0;i++){
                i_origin=_good_point_origins.get_data(i);
                if(i_origin<0){
                    mask.set(i,1);
                }
                else if(considered.contains(i_origin)==1){
                    for(j=0;considered.get_data(j)!=i_origin;j++);
                    if(considered.get_data(j)!=i_origin){
                        printf("WARNING you did the indexing wrong in create_mask\n");
                        exit(1);
                    }
                    mask.set(i,is_associate.get_data(j));
                }
                else{
                    j=check_association(i_pt, i_origin);
                    considered.add(i_origin);
                    is_associate.add(j);
                    mask.set(i,j);
                }
            }
        }

};

#endif
