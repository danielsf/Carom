#ifndef NODE_H
#define NODE_H

#include "wrappers.h"

class node{

public:
    
    node();
    node(const node&);
    node& operator=(const node&);
    ~node();
    
    void set_chisquared(chisquared*);
    void set_center(int);
    void evaluate(array_1d<double>&, double*, int*);
    int get_center();
    
private:
    int _centerdex;
    double _chimin;
    
    array_1d<int> _compass_points,_basis_associates;
    array_1d<double> _basis_mm,_basis_bb,_basis_model;
    array_2d<double> _basis_vectors;
    
    chisq_wrapper *_chisquared;
    
    void bisection(array_1d<double>&,double,array_1d<double>&,double);
    
    void perturb_bases(array_2d<double>&,int,array_1d<double>&);
    double basis_error(array_2d<double>&,array_1d<double>&);
    void initialize();
    void copy(const node&);
    
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
    node *data;
    int _ct,_room;
};


#endif
