#ifndef CHAIN_H
#define CHAIN_H

#include "containers.h"
#include "goto_tools.h"

class chain{

public:

    ~chain();

    chain();
    chain(int);
    chain(int, char*);
    
    void read_chain(char*);
    
    void set_dim(int);
    void set_dice(Ran*);
    void set_output_name(char*);
    
    int is_current_point_valid();
    double get_current_point(int);
    double get_current_chisquared();
    int get_current_degeneracy();
    
    int get_points();
    int get_rows();
    int get_dim();
    double get_point(int,int);
    double get_chisquared(int);
    int get_degeneracy(int);
    
    void add_point(array_1d<double>&, double);
    
    //write routine to get thinned dexes
    
    void get_thinned_indices(int,int,array_1d<int>&,int);
    int get_thinby(double, int, int, int);
    
    void get_thinned_indices(int,int,array_1d<int>&);
    int get_thinby(double, int, int);
    
    void write_chain();
    void write_burnin();
    void copy(const chain&);

private:

    int _dim,_n_written,_n_written_burnin;
    int _current_degeneracy;
    
    double _current_chi;
    
    array_2d<double> _points;
    array_1d<int> _degeneracy;
    array_1d<double> _chisquared,_current_point;

    Ran *_dice;
    
    char _output_name[letters];

    void initialize();
    void write(char*,int);
    void verify_dim(int,char*);
    void verify_points(int,char*);
    void is_dice_safe(char*);

};

class arrayOfChains{

public:

    ~arrayOfChains();
    arrayOfChains();
    arrayOfChains(int,int,Ran*);
    arrayOfChains(array_2d<double>&, array_1d<double>&,Ran*);

    int get_n_chains();
    void get_covariance_matrix(double,int,array_2d<double>&);

    void initialize(int,int,Ran*);

    void add(array_1d<double>&, double);
    void add(array_2d<double>&, array_1d<double>&);
    void remove(int);

    chain* operator()(int);
    
    void calculate_R(array_1d<double>&, array_1d<double>&, array_1d<double>&);

private:
    chain *_data;
    Ran *_dice;
    int _n_chains;
    int _dim;
    
    asymm_array_2d<int> _independent_samples;

    void verify_chains(int,char*);

};

#endif
