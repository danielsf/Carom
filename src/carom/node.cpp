#include "node.h"

node::node(){
    initialize();
}

node::~node(){}

node::node(const node &in){
    initialize();
    copy(in);
}

node& node::operator=(const node&in){
    if(this==&in) return *this;
    initialize();
    copy(in);
    return *this;
}

void node::initialize(){
    _chisquared=NULL;
    _chimin=2.0*exception_value;
    _centerdex=-1;
    _compass_points.set_name("node_compass_points");
    _basis_associates.set_name("node_basis_associates");
    _basis_mm.set_name("node_basis_mm");
    _basis_bb.set_name("node_basis_bb");
    _basis_model.set_name("node_basis_model");
    _basis_vectors.set_name("node_basis_vectors");
}

void node::copy(const node &in){
    _centerdex=in._centerdex;
    _chimin=in._chimin;
    
    int i,j;
    
    _chisquared=in._chisquared;
    
    _compass_points.reset();
    for(i=0;i<in._compass_points.get_dim();i++){
        _compass_points.set(i,in._compass_points.get_data(i));
    }
    
    _basis_associates.reset();
    for(i=0;i<in._basis_associates.get_dim();i++){
        _basis_associates.set(i,in._basis_associates.get_data(i));
    }
    
    _basis_mm.reset();
    for(i=0;i<in._basis_mm.get_dim();i++){
        _basis_mm.set(i,in._basis_mm.get_data(i));
    }
    
    _basis_bb.reset();
    for(i=0;i<in._basis_bb.get_dim();i++){
        _basis_bb.set(i,in._basis_bb.get_data(i));
    }
    
    _basis_model.reset();
    for(i=0;i<in._basis_model.get_dim();i++){
        _basis_model.set(i,in._basis_model.get_data(i));
    }
    
    _basis_vectors.reset();
    _basis_vectors.set_cols(in._basis_vectors.get_cols());
    for(i=0;i<in._basis_vectors.get_rows();i++){
        for(j=0;j<in._basis_vectors.get_cols();j++){
            _basis_vectors.get_data(i,j,in._basis_vectors.get_data(i,j));
        }
    }
}

///////////////arrayOfNodes code below//////////

arrayOfNodes::arrayOfNodes(){
    data=NULL;
    ct=0;
    room=0;
}

arrayOfNodes::~arrayOfNodes(){
    if(data!=NULL){
        delete [] data;
    }
}

void arrayOfNodes::add(int cc, chisq_wrapper *gg){

    node *buffer;
    int i,j;
    
    if(_ct==_room){
        if(_ct>0){
            buffer=new node[ct];
            for(i=0;i<ct;i++){
                buffer[i].copy(data[i]);
            }
            
            delete [] data;
        }
        
        _room+=5;
        data=new node[_room];
        
        if(_ct>0){
            for(i=0;i<_ct;i++){
                data[i].copy(buffer[i]);
            }
            delete [] buffer;
        }
    }
    
    data[ct].set_chisquared(gg);
    data[ct].set_center(cc);

    _ct++;

}

void arrayOfNodes::add(gpWrapper *g, int i){
    add(i,g);
}

int arrayOfNodes::get_dim(){
    return _ct;
}

void arrayOfNodes::remove(int ii){
    int i;
    if(ii>=_ct) return;
    
    for(i=ii+1;i<_ct;i++){
        data[i-1].copy(data[i]);
    }
    _ct--;
}

node* arrayOfNodes::operator()(int dex){
    if(dex<0 || dex>=_ct){
        printf("WARNING asked arrayOfNodes for dex %d but only have %d (operator)\n",
        dex,_ct);
        
        exit(1);
    }
    return &data[dex];
}
