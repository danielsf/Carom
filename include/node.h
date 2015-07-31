#ifndef NODE_H
#define NODE_H

#include "chisq_wrapper.h"
#include "eigen_wrapper.h"
#include "simplex.h"

class node{

public:
    
    node();
    node(const node&);
    node& operator=(const node&);
    ~node();
    void copy(const node&);
    
    void set_chisquared(chisq_wrapper*);
    void set_center(int);
    void set_basis(int,int,double);
    void evaluate(array_1d<double>&, double*, int*);
    int get_center();
    void find_bases();
    void search();
    void ricochet();
    void simplex_search();
    double volume();
    double projected_volume();
    double distance_traveled(int);
    int get_activity();
    int get_ct_ricochet();
    int get_n_particles();
    int get_n_candidates();
    
    void project_to_bases(array_1d<double>&,array_1d<double>&);
    void recalibrate_projected_max_min();
    void print_ricochet_discoveries(char*);
    
    int is_this_an_associate(int);
    
    int get_failed_kicks();
    int get_successful_kicks();
    int get_strikeouts();
    int get_successful_ricochets();
    
private:
    int _centerdex,_geo_centerdex,_centerdex_basis,_active,_found_bases,_ellipse_center;
    int _min_changed,_allowed_ricochet_strikes,_failed_simplexes;
    int _ct_ricochet,_ct_simplex,_calls_to_ricochet;
    double _chimin,_chimin_bases;
    double _volume_of_last_geom;
    
    int _failed_kicks,_successful_kicks;
    
    double _min_basis_error;
    int _since_expansion,_min_basis_error_changed;
    int _strikeouts,_successful_ricochets;
    
    array_1d<int> _compass_points,_basis_associates,_off_center_compass_points;
    array_1d<int> _ricochet_candidates,_off_center_origins;
    array_1d<int> _associates;
    array_1d<double> _basis_mm,_basis_bb,_basis_model,_basis_vv;
    array_1d<double> _basis_lengths;
    array_1d<double> _max_found,_min_found,_distance_traveled;
    array_1d<double> _projected_min,_projected_max;
    array_2d<double> _basis_vectors,_basis_ddsq;
    
    array_2d<double> _ricochet_velocities,_ricochet_candidate_velocities;
    array_1d<int> _ricochet_particles;
    array_1d<int> _ricochet_origins;
    array_1d<int> _ricochet_strikes;
    array_1d<int> _boundary_points;
    
    asymm_array_2d<int> _ricochet_discoveries,_ricochet_discovery_time;
    array_1d<int> _ricochet_discovery_dexes;
    asymm_array_2d<double> _ricochet_distances,_ricochet_grad_norm;
    asymm_array_2d<double> _ricochet_dir_norm,_ricochet_mu;
    asymm_array_2d<int> _ricochet_strike_log;
    
    chisq_wrapper *_chisquared;
    
    int node_bisection(array_1d<double>&,double,array_1d<double>&,double,int);
    int node_bisection(array_1d<double>&,double,array_1d<double>&,double,int,double,double);
    
    void perturb_bases(int,array_1d<double>&,array_2d<double>&);
    double basis_error(array_2d<double>&,array_1d<double>&);
    int findAcceptableCenter();
    void findCovarianceMatrix(int,array_2d<double>&);
    void guess_bases(array_2d<double>&);
    void validate_bases(array_2d<double>&, char*);
    void initialize();
    void add_to_boundary(int);
    
    void is_it_safe(char*);
    
    void compass_search();
    void compass_search(int);
    void compass_diagonal(int);
    void compass_search_geometric_center();
    void off_center_compass(int);
    void initialize_ricochet();
    double ricochet_distance(int,int);
    double ricochet_model(array_1d<double>&, kd_tree&);
    double apply_quadratic_model(array_1d<double>&);
    
    int kick_particle(int, array_1d<double>&);
    int step_kick(int, double, array_1d<double>&);
    int t_kick(int,array_1d<double>&);
    int smart_step_kick(int, double, array_1d<double>&);
    void originate_particle_compass(int, array_1d<double>&);
    void originate_particle_shooting(int, array_1d<double>&);
    void _originate_particle_paperwork(int, int);
    void _filter_candidates();
    
    double node_distance(array_1d<double>&, array_1d<double>&);
    double node_distance(int, int);
    double node_distance(int, array_1d<double>&);
    void node_gradient(int, array_1d<double>&);
    double node_second_derivative(int,int,int);
    double node_second_derivative_different(int,int,int);
    double node_second_derivative_same(int,int);

    int _are_connected(int, int);
    
    int is_it_a_strike(int, kd_tree&);
};

class arrayOfNodes{

public:
    arrayOfNodes();
    ~arrayOfNodes();
        
    void add(int,chisq_wrapper*);
    void add(chisq_wrapper*,int);
    int get_dim();
    void remove(int);
    node* operator()(int);
        
private:
    node *_data;
    int _ct,_room;
};


#endif
