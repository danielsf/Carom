#ifndef NODE_H
#define NODE_H

#include "wrappers.h"

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
    
private:
    int _centerdex;
    int _ricochet_since_expansion,_ricochet_since_perturbation,_min_changed;
    double _chimin,_bisection_tolerance;
    
    array_1d<int> _compass_points,_basis_associates;
    array_1d<double> _basis_mm,_basis_bb,_basis_model,_basis_vv;
    array_1d<double> _max_found,_min_found;
    array_2d<double> _basis_vectors,_basis_ddsq;
    
    array_2d<double> _ricochet_particles,_ricochet_velocities;
    
    chisq_wrapper *_chisquared;
    
    int bisection(array_1d<double>&,double,array_1d<double>&,double);
    
    void perturb_bases(int,array_1d<double>&,array_2d<double>&);
    double basis_error(array_2d<double>&,array_1d<double>&);
    void initialize();
    
    void is_it_safe(char*);
    
    void compass_search();
    void compass_off_diagonal();
    void initialize_ricochet();
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
