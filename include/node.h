#ifndef NODE_H
#define NODE_H

#include "wrappers.h"
#include "eigen_wrapper.h"

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
    void ricochet();
    double volume();
    int get_activity();
    int get_ct_ricochet();
    
private:
    int _centerdex,_active,_found_bases,_ellipse_center;
    int _ricochet_since_expansion,_min_changed;
    int _allowed_ricochet_strikes;
    int _ct_ricochet,_calls_to_ricochet;
    double _chimin,_bisection_tolerance;
    
    array_1d<int> _compass_points,_basis_associates,_off_center_compass_points;
    array_1d<double> _basis_mm,_basis_bb,_basis_model,_basis_vv;
    array_1d<double> _basis_lengths;
    array_1d<double> _max_found,_min_found;
    array_2d<double> _basis_vectors,_basis_ddsq;
    
    array_2d<double> _ricochet_particles,_ricochet_velocities;
    array_1d<int> _ricochet_strikes;
    
    chisq_wrapper *_chisquared;
    
    int bisection(array_1d<double>&,double,array_1d<double>&,double,int);
    
    void perturb_bases(int,array_1d<double>&,array_2d<double>&);
    double basis_error(array_2d<double>&,array_1d<double>&);
    void initialize();
    
    void is_it_safe(char*);
    
    void compass_search();
    void compass_off_diagonal();
    void off_center_compass(int);
    void initialize_ricochet();
    double ricochet_distance(int,int);
    double ricochet_model(array_1d<double>&, kd_tree&);
    double apply_quadratic_model(array_1d<double>&);
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
