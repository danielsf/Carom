#ifndef MCMC_H
#define MCMC_H

#include "containers.h"
#include "goto_tools.h"
#include "eigen_wrapper.h"
#include "wrappers.h"
#include "time.h"

class mcmc_sampler{

public:

    mcmc_sampler(){
        _chisq_fn=NULL;
        _dice=NULL;
        _n_chains=0;
        _write_every=1000;
        _do_learning=0;
        _dim=0;
        _has_written=0;
        _seed=-1;

        _points.set_name("mcmc_points");
        _point_dexes.set_name("mcmc_point_dexes");
        _bases.set_name("mcmc_bases");
        _radii.set_name("mcmc_radii");
        _chisq_arr.set_name("mcmc_chisq_arr");

        sprintf(_name_root,"mcmc_chain");
    }

    ~mcmc_sampler(){
        if(_dice!=NULL){
            delete _dice;
        }
    }

    void set_seed(int ii){
        if(_dice!=NULL){
            printf("WARNING setting seed, but _dice not NULL\n");
            exit(1);
        }
        if(ii>0){
            _seed=ii;
        }
        else{
            _seed=int(time(NULL));
        }
        _dice=new Ran(_seed);
    }

    void set_write_every(int ii){
        _write_every=ii;
    }

    void set_name_root(char *word){
        sprintf(_name_root,"%s",word);
    }

    void do_learning(){
        _do_learning=1;
        printf("WARNING learning of bases not implemented yet\n");
        exit(1);
    }
    void do_not_do_learning(){
        _do_learning=0;
    }

    void set_dim(int ii){
        _dim=ii;
    }

    void set_n_chains(int ii){
        _n_chains=ii;
    }

    void set_bases(const array_2d<double> &b_in, const array_1d<double> &r_in){
        if(_dim>0){
            if(_dim!=b_in.get_rows() || _dim!=b_in.get_cols()){
                printf("WARNING setting bases with wrong dim\n");
                exit(1);
            }
            if(_dim!=r_in.get_dim()){
                printf("WARNING setting radii with wrong dim\n");
                exit(1);
            }
        }
        _dim = b_in.get_cols();
        int i;
        for(i=0;i<_dim;i++){
            _bases.add_row(b_in(i));
            _radii.add(r_in.get_data(i));
        }
    }

    void set_chisq_fn(function_wrapper *fn){
        if(_chisq_fn!=NULL){
            printf("WARNING; overwriting chisq_fn\n");
            exit(1);
        }
        _chisq_fn=fn;
    }

    void set_start_pts(const array_2d<double> &pts_in){
        if(_n_chains>0){
            if(_n_chains!=pts_in.get_rows()){
                printf("WARNING setting start points with wrong number of points\n");
                exit(1);
            }
        }
        _n_chains=pts_in.get_rows();
        if(_dim>0){
            if(_dim!=pts_in.get_cols()){
                printf("WARNING setting start points with wrong dimensionality\n");
                exit(1);
            }
        }
        _dim=pts_in.get_cols();
        _points.reset_preserving_room();
        _point_dexes.reset_preserving_room();
        _chisq_arr.reset_preserving_room();
        int i;
        double mu;
        for(i=0;i<_n_chains;i++){
            _points.add_row(pts_in(i));
            _chisq_arr.add(_chisq_fn[0](pts_in(i)));
            _point_dexes.set(i,0,_points.get_rows()-1);
        }
        if(_points.get_rows()!=_chisq_arr.get_dim()){
            printf("WARNING after setting start points, %d points but %d chisq\n",
                   _points.get_rows(),_chisq_arr.get_dim());
            exit(1);
        }
        _has_written=0;
    }

    void sample(int);
    void write_chains();

private:

    function_wrapper *_chisq_fn;
    Ran *_dice;
    int _n_chains;
    int _write_every;
    int _do_learning;
    int _has_written;
    int _dim;
    int _seed;
    char _name_root[letters];
    array_2d<double> _points;
    array_1d<double> _chisq_arr;
    asymm_array_2d<int> _point_dexes;
    array_2d<double> _bases;
    array_1d<double> _radii;

};

#endif
