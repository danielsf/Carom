#ifndef MAPS_H
#define MAPS_H

#include <stdio.h>
#include <time.h>
#include "eigen_wrapper.h"
#include "chisq_wrapper.h"
#include "simplex.h"
#include "gp_lin.h"
#include "dchi_simplex.h"
#include "dchi_simplex_gp.h"
#include "search_types.h"
#include "dalex.h"
#include "maps_initializer.h"

class maps{

public:

    maps();
    ~maps();

    void initialize(int);

    void set_seed(int);
    void set_min(array_1d<double>&);
    void set_max(array_1d<double>&);
    void set_characteristic_length(int,double);
    void set_deltachi(double);
    void set_target(double);
    void set_write_every(int);
    void set_outname(char*);
    void set_timingname(char*);

    void set_chisquared(chisquared*);

    void search(int);

    void mcmc_init();

    void explore();
    void write_pts();
    void write_log();

    int get_called();
    double get_chimin();

    void set_confidence_limit(double);
    void set_dof(int);

    double evaluate(array_1d<double>&, int*);
    void assess_good_points();
    void assess_good_points(int);
    void assess_good_points(int,int);

private:

    chisq_wrapper _chifn;
    gp_lin _interpolator;
    int _write_every,_last_written;
    int _ct_simplex_boundary,_ct_simplex_min,_calls_to_simplex_boundary;
    int _ct_dalex;
    int _last_wrote_log;
    int _simplex_mindex;

    dalex _cloud;
    double _mcmc_basis_min;

    asymm_array_2d<int> _log;

    double _time_started;

    char _outname[letters],_timingname[letters];

    array_1d<int> _good_points,_duds,_duds_for_min;
    array_1d<int> _failed_mins;
    array_1d<int> _explorers;

    int _last_did_min;
    double _explorer_temp;

};

#endif
