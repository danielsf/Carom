#include "wrappers.h"

function_wrapper::function_wrapper(){}

function_wrapper::~function_wrapper(){}

double function_wrapper::operator()(array_1d<double> &vv){
    printf("WARNING calling un-implemented function_wrapper operator\n");
    exit(1);
}

int function_wrapper::get_called(){
    printf("WARNING calling un-implemented function_wrapper get_called\n");
    exit(1);
}

////////////////////////////

chisq_wrapper::chisq_wrapper(){
    _kptr=NULL;
    _dice=NULL;
    _chifn=NULL;
    _fn.set_name("chisq_wrapper_fn");
    _characteristic_length.set_name("chisq_wrapper_characteristic_length");
    _range_min.set_name("chisq_wrapper_range_min");
    _range_max.set_name("chisq_wrapper_range_max");
    _valid_neigh.set_name("chisq_wrapper_valid_neigh");
    _valid_dd.set_name("chisq_wrapper_valid_dd");
    _ddmin=1.0e-8;
    _adaptive_target=1;
    _delta_chi=-1.0;
    _chimin=2.0*exception_value;
    _target=2.0*exception_value;
    _seed=-1;
    _called=0;
}

chisq_wrapper::~chisq_wrapper(){
    if(kptr!=NULL){
        delete kptr;
    }
    
    if(dice!=NULL){
        delete dice;
    }
}

void chisq_wrapper::set_chisquared(chisquared *xx){
    chifn=xx;
}

void chisq_wrapper::set_seed(int i){
    _seed=i;
}

void chisq_wrapper::set_target(double xx){
    _target=xx;
    _adaptive_target=0;
}

void chisq_wrapper::set_ddmin(double dd){
    _ddmin=dd;
}

void chisq_wrapper::set_min(array_1d<double> &vv){
    if(kptr!=NULL){
        printf("WARNING in chisq_wrapper setting min but kptr not null\n");
        exit(1);
    }
    
    int i;
    _range_min.reset();
    for(i=0;i<vv.get_dim();i++){
        _range_min.set(i,vv.get_data(i));
    }
}

void chisq_wrapper::set_max(array_1d<double> &vv){
    if(kptr!=NULL){
        printf("WARNINGin chisq_wrapper setting max but kptr not null\n");
        exit(1);
    }

    int i;
    _range_max.reset();
    for(i=0;i<vv.get_dim();i++){
        _range_max.set(i,vv.get_data(i));
    }
}

void chisq_wrapper::set_characteristic_length(int dex, double xx){
    if(kptr!=NULL){
        printf("WARNING in chisq_wrapper setting characteristic_length but kptr not null\n");
        exit(1);
    }

    int i;
    if(dex>_characteristic_length.get_dim()){
        for(i=_characteristic_length.get_dim();i<dex;i++){
            _characteristic_length.set(i,-1.0);
        }
    }
    
    _characteristic_length.set(dex,xx);
}

void chisq_wrapper::set_deltachi(double xx){
    if(_adaptive_target!=1){
        printf("WARNING chisq_wrapper trying to set detlachi, but not an adaptive target\n");
        exit(1);
    }
    _deltachi=xx;
}

void chisq_wrapper::initialize(int npts, int dim){
    if(chifn==NULL){
        printf("WARNING calling chisq_wrapper_initialize with null chifn\n");
        exit(1);
    }
    
    if(kptr!=NULL){
        printf("WARNING calling chisq_wrapper_initialize even though kptr not null\n");
        exit(1);
    }
    
    if(_range_min.get_dim()!=dim){
        printf("WARNING chisq_wrapper_initialize dim %d range_min %d\n",
        dim,_range_min.get_dim());
        
        exit(1);
    }
    
    if(_range_max.get_dim()!=dim){
        printf("WARNING chisq_wrapper_initialize dim %d range_max %d\n",
        dim,_range_max.get_dim());
        
        exit(1);
    }
    
    if(dice==NULL){
        if(_seed<0)_seed=int(time(NULL));
        dice=new Ran(_seed);
    }
    
    array_2d<double> data;
    array_1d<double> vv;
    vv.set_name("chisq_wrapper_initialize_vv");
    data.set_name("chisq_wrapper_data");
    int i,j;
    double mu;
    
    _fn.reset();
    data.set_cols(dim);
    for(i=0;i<npts;i++){
        for(j=0;j<dim;j++){
            vv.set(j,_range_min.get_data(j)+dice->doub()*(_range_max.get_data(j)-_range_min.get_data(j)));
        }
        mu=chifn[0](vv);
        _called++;
        if(mu<_chimin){
            _chimin=mu;
            _mindex=i;
        }
        _fn.add(mu);
        data.add_row(vv);
    }
    
    array_1d<double> temp_max,temp_min;
    temp_max.set_name("chisq_wrapper_initialize_temp_max");
    temp_min.set_name("chisq_wrapper_initialize_temp_min");
    for(i=0;i<dim;i++){
        if(_characteristic_length.get_dim()>i && _characteristic_length.get_data(i)>0.0){
            temp_min.set(i,0.0);
            temp_max.set(i,_characteristic_length.get_data(i));
        }
        else{
            temp_min.set(i,_range_min.get_data(i));
            temp_max.set(i,_range_max.get_data(i));
        }
    }

    kptr=new kd_tree(data,temp_min,temp_max);
    
    if(_adaptive_target==1){
        if(_deltachi<0.0){
            printf("WARNING when initializing chisq_wrapper deltachi %e\n",_deltachi);
            exit(1);
        }
        _target=_chimin+_deltachi;
    }
}

void chisq_wrapper::is_it_safe(char *word){
    if(kptr==NULL){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("kptr is null\n");
        exit(1);
    }
    
    if(dice==NULL){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("dice is null\n");
        exit(1);
    }
    
    if(chifn==NULL){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("chifn is null\n");
        exit(1);
    }
    
    if(_adaptive_target==1 && _delta_chi<0){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("adaptive target but delta_chi %e\n",_deltachi);
        exit(1);
    }
    
    if(_fn.get_dim()!=kptr->get_pts()){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("fn dim %d kptr pts %d\n",_fn.get_dim(),kptr->get_pts());
        exit(1);
    }
}

double chisq_wrapper::target(){
    return _target;
}

int chisq_wrapper::is_valid(array_1d<double> &pt, int *neighdex){
    is_it_safe("is_valid");
    
    int i;
    neighdex[0]=-1;
    for(i=0;i<pt.get_dim();i++){
        if(pt.get_data(i)>chifn->get_max(i)){
            return 0;
        }
        if(pt.get_data(i)<chifn->get_min(i)){
            return 0;
        }
    }

    kptr->nn_srch(pt,1,_valid_neigh,_valid_dd);
    if(_valid_dd.get_data(0)<_ddmin){
        neighdex[0]=_valid_neigh.get_data(0);
        return 0;
    }
    
    return 1;
    
}

void chisq_wrapper::evaluate(array_1d<double> &pt, double *value, int *dex){
    is_it_safe("evaluate");
    
    int validity,neighdex;
    validity=is_valid(pt,&neighdex);
    
    if(validity!=1){
        if(neighdex>=0){
            value[0]=_fn.get_data(neighdex);
            dex[0]=neighdex;
        }
        else{
            value[0]=2.0*exception_value;
            dex[0]=-1;
        }
        return;
    }
    
    double mu;
    mu=chifn[0](pt);
    _called++;
    
    if(mu<exception_value){
        kptr->add(pt);
        _fn.add(mu);
    
        if(mu<_chimin){
            _chimin=mu;
            _mindex=kptr->get_pts();
            if(_adaptive_target==1){
                _target=_chimin+_deltachi;
            }
        }
    }
}

int chisq_wrapper::get_pts(){
    is_it_safe("get_pts");
    return kptr->get_pts();
}

double chisq_wrapper::random_double(){
    if(dice==NULL){
        printf("WARNING chisq_wrapper random_double dice is null\n");
        exit(1);
    }
    
    return dice->doub();
}

int chisq_wrapper::random_int(){
    if(dice==NULL){
        printf("WARNING chisq_wrapper random_int dice is null\n");
        exit(1);
    }
    
    return dice->int32();
}

double chisq_wrapper::get_fn(int dex){
    if(dex<0 || dex>_fn.get_dim()){
        printf("WARNING asking for fn %d but %d\n",dex,_fn.get_dim());
    }
    
    return _fn.get_data(dex);
}

void chisq_wrapper::nn_srch(array_1d<double> &vv, int kk, array_1d<int> &neigh, array_1d<double> &dd){
    is_it_safe("nn_srch");
    kptr->nn_srch(vv,kk,neigh,dd);
}
