#ifndef NODE_H
#define NODE_H

#include "chisq_wrapper.h"
#include "eigen_wrapper.h"
#include "simplex.h"
#include "dchi_simplex.h"

class node{

public:

    node();
    node(const node&);
    node& operator=(const node&);
    ~node();
    void copy(const node&);
    void merge(const node&);
    void deactivate();

    void initialize_ricochet();

    void set_chisquared(chisq_wrapper*);
    void set_center(int);
    void set_basis(int,int,double);
    void evaluate(array_1d<double>&, double*, int*);
    int get_center();
    void find_bases();
    void compass_search();
    void search();
    void ricochet();

    void swarm_shoot(int);
    void swarm_evaluate(array_1d<double>&,double*);
    void swarm_search();

    void _shift_ricochet(int);
    int _ricochet(int, array_1d<double>&);
    int _adaptive_ricochet(int, array_1d<double>&);
    void simplex_search();
    double volume();
    double get_norm(int);
    double get_projected_norm(int);
    double get_transform(int);
    double projected_volume();
    int get_activity();
    int get_ct_ricochet();
    int get_n_particles();
    int get_n_candidates();

    int get_n_boundary();
    int get_boundary(int);
    int get_n_associates();
    int get_associate(int);

    void project_to_bases(array_1d<double>&,array_1d<double>&);
    void recalibrate_max_min();

    int is_this_an_associate(int);
    int is_this_an_associate_gross(int);
    int find_local_center();

    int get_highball_calls();
    int get_total_bisections();
    int get_bisection_calls();
    int get_total_ricochets();
    int get_ricochet_calls();
    int get_ricochet_bisection_calls();
    int get_ricochet_bisections();
    int get_gradient_calls();

    int get_convergence_ct();
    int get_total_trimmed();
    int get_total_kicks();
    int get_failed_kicks();
    int get_successful_kicks();
    int get_good_shots();
    int get_bad_shots();
    int get_strikeouts();
    int get_successful_ricochets();
    int get_proper_ricochets();
    void set_id_dex(int);
    void write_node_log(char*);
    int get_swarm_expand();
    int get_swarm_outside();
    void deactivate_simplex();

    double get_ricochet_growth();
    double get_swarm_growth();
    double get_mcmc_growth();

private:
    int _do_simplex;
    int _id_dex,_last_wrote_log;
    int _first_centerdex;
    int _centerdex,_geo_centerdex,_centerdex_basis,_active,_found_bases,_ellipse_center;
    int _min_changed,_allowed_ricochet_strikes,_failed_simplexes;
    int _ct_ricochet,_ct_simplex;
    double _chimin,_chimin_bases,_chimin_ricochet;
    double _volume_of_last_geom;

    int _compass_calls;
    int _total_bisections,_bisection_calls,_total_ricochets,_ricochet_calls;
    int _highball_calls;
    int _ricochet_bisection_calls,_gradient_calls,_ricochet_bisections;
    int _failed_kicks,_successful_kicks,_total_kicks,_total_trimmed;
    int _proper_ricochets;

    double _min_basis_error;
    double _node_dd_tol;
    int _since_expansion,_min_basis_error_changed,_convergence_ct;
    int _successful_ricochets,_good_shots,_bad_shots;

    array_1d<int> _compass_points,_basis_associates,_off_center_compass_points;
    asymm_array_2d<int> _transform_associates;
    array_1d<int> _firework_centers;
    array_1d<int> _ricochet_candidates;
    array_1d<int> _avg_pts;
    array_1d<int> _off_center_origins,_off_center_candidates;
    array_1d<int> _associates;

    array_1d<double> _basis_mm,_basis_bb,_basis_model,_basis_vv;
    array_1d<double> _basis_lengths;
    array_1d<double> _max_found,_min_found;
    array_1d<double> _projected_min,_projected_max;
    array_2d<double> _basis_vectors,_basis_ddsq;

    array_2d<double> _swarm;
    array_1d<double> _swarm_center;
    array_1d<double> _swarm_norm,_swarm_max,_swarm_min;
    array_1d<int> _swarm_associates;
    int _swarm_acceptances,_swarm_rejections;
    double _swarm_step,_f_swarm_center;
    int _swarm_outsiders,_swarm_expanders;

    double _mcmc_step;
    int _mcmc_acceptances,_mcmc_rejections;

    array_1d<double> _true_min,_true_max,_transform;

    array_2d<double> _ricochet_velocities,_ricochet_candidate_velocities;
    array_1d<int> _ricochet_particles;
    array_1d<int> _ricochet_origins;
    array_1d<int> _boundary_points;
    array_1d<int> _ricochet_strikes;
    double _v0;
    int _since_culled;

    array_1d<int> _ricochet_log;
    array_1d<int> _ricochet_candidate_log;
    array_1d<int> _compass_log;
    array_1d<int> _wander_log;

    double _ricochet_growth,_mcmc_growth,_swarm_growth;

    chisq_wrapper *_chisquared;

    int node_bisection(array_1d<double>&,double,array_1d<double>&,double,int);
    int node_bisection(array_1d<double>&,double,array_1d<double>&,double,int,double,double);
    int node_bisection(int,int,double,double);
    int node_bisection_origin_dir(int,array_1d<double>&);
    int node_bisection_origin_dir(int,array_1d<double>&,double,double);

    double evaluate_dir(int,array_1d<double>&);

    void perturb_bases(int,array_1d<double>&,array_2d<double>&);
    double basis_error(array_2d<double>&,array_1d<double>&);
    int findAcceptableCenter();
    void findCovarianceMatrix(int,array_2d<double>&);
    void guess_bases(array_2d<double>&);
    void validate_bases(array_2d<double>&, char*);
    void initialize();
    void add_to_boundary(int);
    void add_to_compass(int);

    void is_it_safe(char*);

    void populate_basis_associates();

    void compass_search(int);
    void compass_diagonal(int);
    void compass_search_geometric_center();
    void off_center_compass(int);
    void trim_ricochet(int);
    void cull_ricochet();
    void remove_particle(int);
    void set_particle(int,int,array_1d<double>&);
    double ricochet_model(array_1d<double>&);
    double ricochet_model(array_1d<double>&,array_1d<int>&);
    double ricochet_model(array_1d<double>&,double*);
    double ricochet_model(array_1d<double>&, kd_tree&);
    double ricochet_model(array_1d<double>&, kd_tree&, array_1d<int>&);
    double ricochet_model(array_1d<double>&, kd_tree&, double*);
    double _ricochet_model(array_1d<double>&, kd_tree&, double*, int, array_1d<int>&);
    double apply_quadratic_model(array_1d<double>&);
    double apply_quadratic_model(int);

    void mcmc_walk(int, int*, array_1d<double>&, int, array_1d<int>&);

    int mcmc_kick(int, int*, array_1d<double>&, int);
    int originate_particle_compass(array_1d<double>&);
    int originate_particle_shooting(array_1d<double>&);
    void originate_particle_simplex();
    void _originate_particle_paperwork(int, int);
    void _filter_candidates();
    double _nearest_other_particle(int,int);

    double node_distance(array_1d<double>&, array_1d<double>&);
    double node_distance(int, int);
    double node_distance(int, array_1d<double>&);
    double normalized_node_distance(int,int);
    void node_gradient(int, array_1d<double>&);
    void _node_2sided_gradient(int,array_1d<double>&);
    void _node_1sided_gradient(int,array_1d<double>&);
    int _node_1sided_projected_gradient(int,array_1d<double>&);

    double node_second_derivative(int,int,int);
    double node_second_derivative_different(int,int,int);
    double node_second_derivative_same(int,int);

    int _are_connected(int, int);

    int is_it_a_strike(int, kd_tree&);
    double get_pt(int,int);
    void transform_pt_to_node(array_1d<double>&,array_1d<double>&);
    void transform_pt_to_truth(array_1d<double>&,array_1d<double>&);
    void get_true_pt(int,array_1d<double>&);
    void set_transform();

    int choose_off_center_point();
};

class arrayOfNodes{

public:
    arrayOfNodes();
    ~arrayOfNodes();

    void add(int,chisq_wrapper*);
    void add(chisq_wrapper*,int);
    int get_dim();
    void cull();
    void remove(int);
    node* operator()(int);

private:
    node *_data;
    int _ct,_room;
};


#endif
