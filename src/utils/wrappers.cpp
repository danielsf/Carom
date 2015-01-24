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
    _ct_where.set_name("chisq_wrapper_ct_where");
    _where_log.set_name("chisq_wrapper_where_log");
    _ddmin=1.0e-8;
    _adaptive_target=1;
    _deltachi=-1.0;
    _chimin=2.0*exception_value;
    _target=2.0*exception_value;
    _seed=-1;
    _called=0;
    _iWhere=-1;
    
    _ct_where.set(iSimplex,0);
    _ct_where.set(iRicochet,0);
    _ct_where.set(iCompass,0);
}

chisq_wrapper::~chisq_wrapper(){
    if(_kptr!=NULL){
        delete _kptr;
    }
    
    if(_dice!=NULL){
        delete _dice;
    }
}

void chisq_wrapper::set_chisquared(chisquared *xx){
    _chifn=xx;
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
    if(_kptr!=NULL){
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
    if(_kptr!=NULL){
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
    if(_kptr!=NULL){
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

void chisq_wrapper::initialize(int npts){
    if(_chifn==NULL){
        printf("WARNING calling chisq_wrapper_initialize with null chifn\n");
        exit(1);
    }
    
    if(_kptr!=NULL){
        printf("WARNING calling chisq_wrapper_initialize even though kptr not null\n");
        exit(1);
    }
    
    if(_range_min.get_dim()!=_chifn->get_dim()){
        printf("WARNING chisq_wrapper_initialize dim %d range_min %d\n",
        _chifn->get_dim(),_range_min.get_dim());
        
        exit(1);
    }
    
    if(_range_max.get_dim()!=_chifn->get_dim()){
        printf("WARNING chisq_wrapper_initialize dim %d range_max %d\n",
        _chifn->get_dim(),_range_max.get_dim());
        
        exit(1);
    }
    
    if(_dice==NULL){
        if(_seed<0)_seed=int(time(NULL));
        _dice=new Ran(_seed);
    }
    
    array_2d<double> data;
    array_1d<double> vv;
    vv.set_name("chisq_wrapper_initialize_vv");
    data.set_name("chisq_wrapper_data");
    int i,j;
    double mu;
    
    _fn.reset();
    data.set_cols(_chifn->get_dim());
    for(i=0;i<npts;i++){
        for(j=0;j<_chifn->get_dim();j++){
            vv.set(j,_range_min.get_data(j)+_dice->doub()*(_range_max.get_data(j)-_range_min.get_data(j)));
        }
        mu=_chifn[0](vv);
        _called++;
        if(mu<_chimin){
            _chimin=mu;
            _mindex=i;
        }
        _fn.add(mu);
        data.add_row(vv);
        _where_log.set(i,iSimplex);
    }
    
    array_1d<double> temp_max,temp_min;
    temp_max.set_name("chisq_wrapper_initialize_temp_max");
    temp_min.set_name("chisq_wrapper_initialize_temp_min");
    for(i=0;i<_chifn->get_dim();i++){
        if(_characteristic_length.get_dim()>i && _characteristic_length.get_data(i)>0.0){
            temp_min.set(i,0.0);
            temp_max.set(i,_characteristic_length.get_data(i));
        }
        else{
            temp_min.set(i,_range_min.get_data(i));
            temp_max.set(i,_range_max.get_data(i));
        }
    }

    _kptr=new kd_tree(data,temp_min,temp_max);
    
    if(_adaptive_target==1){
        if(_deltachi<0.0){
            printf("WARNING when initializing chisq_wrapper deltachi %e\n",_deltachi);
            exit(1);
        }
        _target=_chimin+_deltachi;
    }
}

void chisq_wrapper::is_it_safe(char *word){
    if(_kptr==NULL){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("kptr is null\n");
        exit(1);
    }
    
    if(_dice==NULL){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("dice is null\n");
        exit(1);
    }
    
    if(_chifn==NULL){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("chifn is null\n");
        exit(1);
    }
    
    if(_adaptive_target==1 && _deltachi<0){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("adaptive target but deltachi %e\n",_deltachi);
        exit(1);
    }
    
    if(_fn.get_dim()!=_kptr->get_pts()){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("fn dim %d kptr pts %d\n",_fn.get_dim(),_kptr->get_pts());
        exit(1);
    }
    
    if(_chifn->get_dim()!=_kptr->get_dim()){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("chifn dim %d kptr %d\n",_chifn->get_dim(),_kptr->get_dim());
        exit(1);
    }
}

double chisq_wrapper::target(){
    return _target;
}

double chisq_wrapper::get_deltachi(){
    return _deltachi;
}

int chisq_wrapper::is_valid(array_1d<double> &pt, int *neighdex){
    is_it_safe("is_valid");
    
    int i;
    neighdex[0]=-1;
    for(i=0;i<pt.get_dim();i++){
        if(_chifn->get_max(i)>-1.0*exception_value && pt.get_data(i)>_chifn->get_max(i)){
            return 0;
        }
        if(_chifn->get_min(i)<exception_value && pt.get_data(i)<_chifn->get_min(i)){
            return 0;
        }
    }

    _kptr->nn_srch(pt,1,_valid_neigh,_valid_dd);
    if(_valid_dd.get_data(0)<_ddmin){
        neighdex[0]=_valid_neigh.get_data(0);
        return 0;
    }
    
    return 1;
    
}

double chisq_wrapper::operator()(array_1d<double> &pt){
    double mu;
    int dex;
    int i;
    evaluate(pt,&mu,&dex);
    return mu;
}

void chisq_wrapper::set_iWhere(int ii){
    _iWhere=ii;
}

int chisq_wrapper::get_ct_where(int ii){
    return _ct_where.get_data(ii);
}

int chisq_wrapper::get_where_log(int ii){
    return _where_log.get_data(ii);
}

void chisq_wrapper::evaluate(array_1d<double> &pt, double *value, int *dex){
    is_it_safe("evaluate");
    
    int validity,neighdex;
    validity=is_valid(pt,&neighdex);
    
    dex[0]=-1;
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
    mu=_chifn[0](pt);
    value[0]=mu;
    _called++;
    
    if(_iWhere>=0){
        _ct_where.add_val(_iWhere,1);
    }
    
    if(mu<exception_value){
        _kptr->add(pt);
        _fn.add(mu);
        dex[0]=_kptr->get_pts()-1;
        _where_log.set(dex[0],_iWhere);
    
        if(mu<_chimin){
            _chimin=mu;
            _mindex=_kptr->get_pts();
            if(_adaptive_target==1){
                _target=_chimin+_deltachi;
            }
        }
    }
}

int chisq_wrapper::get_called(){
    is_it_safe("get_called");
    return _chifn->get_called();
}

int chisq_wrapper::get_pts(){
    is_it_safe("get_pts");
    return _kptr->get_pts();
}

int chisq_wrapper::get_dim(){
    is_it_safe("get_dim");
    return _kptr->get_dim();
}

double chisq_wrapper::random_double(){
    if(_dice==NULL){
        printf("WARNING chisq_wrapper random_double dice is null\n");
        exit(1);
    }
    
    return _dice->doub();
}

int chisq_wrapper::random_int(){
    if(_dice==NULL){
        printf("WARNING chisq_wrapper random_int dice is null\n");
        exit(1);
    }
    
    return _dice->int32();
}

double chisq_wrapper::get_fn(int dex){
    if(dex<0 || dex>=_fn.get_dim()){
        printf("WARNING asking for fn %d but %d\n",dex,_fn.get_dim());
    }
    
    return _fn.get_data(dex);
}


double chisq_wrapper::get_pt(int dex, int idim){
    is_it_safe("get_pt");

    if(dex<0 || dex>=_fn.get_dim()){
        printf("WARNING asking for pt %d but only have %d \n",
        dex,_fn.get_dim());
        
        exit(1);
    }
    
    if(idim<0 || idim>=_kptr->get_dim()){
        printf("WARNING asking for pt dim %d but only have %d\n",
        idim,_kptr->get_dim());
        
        exit(1);
    }
    
    return _kptr->get_pt(dex,idim);
}

void chisq_wrapper::nn_srch(array_1d<double> &vv, int kk, array_1d<int> &neigh, array_1d<double> &dd){
    is_it_safe("nn_srch");
    _kptr->nn_srch(vv,kk,neigh,dd);
}

Ran* chisq_wrapper::get_dice(){
    is_it_safe("get_dice");
    return _dice;
}

void chisq_wrapper::get_min(array_1d<double> &vv){
    int i;
    vv.set_dim(_range_min.get_dim());
    for(i=0;i<_range_min.get_dim();i++){
        vv.set(i,_range_min.get_data(i));
    }
}

void chisq_wrapper::get_max(array_1d<double> &vv){
    int i;
    vv.set_dim(_range_max.get_dim());
    for(i=0;i<_range_max.get_dim();i++){
        vv.set(i,_range_max.get_data(i));
    }
}

double chisq_wrapper::distance(array_1d<double> &p, int dex){
    is_it_safe("distance(arr,i)");
    return _kptr->distance(p,dex);
}

double chisq_wrapper::distance(int i1, int i2){
    is_it_safe("distance(i,i)");
    return _kptr->distance(i1,i2);
}

void chisq_wrapper::find_gradient(array_1d<double> &pt, array_1d<double> &grad){
    is_it_safe("find_gradient");

    double fCenter;
    int iCenter;
    evaluate(pt,&fCenter,&iCenter);
    
    int ix,i;
    double x2,dx;
    array_1d<double> trial;
    trial.set_name("chisq_wrapper_gradient_trial");
    
    for(i=0;i<_kptr->get_dim();i++){
        trial.set(i,pt.get_data(i));
    }
    
    for(ix=0;ix<_kptr->get_dim();ix++){
        if(_range_max.get_dim()==0 || _range_max.get_data(ix)-_range_min.get_data(ix)>1.0){
            dx=0.01;
        }
        else{
            dx=0.01*(_range_max.get_data(ix)-_range_min.get_data(ix));
        }
        
        trial.add_val(ix,dx);
        evaluate(trial,&x2,&i);
        grad.set(ix,(fCenter-x2)/(pt.get_data(ix)-trial.get_data(ix)));
        trial.set(ix,pt.get_data(ix));
    }
    
}

void chisq_wrapper::copy(chisq_wrapper &in){
    if(in._chifn==NULL){
        printf("WARNING cannot copy chisq_wrapper chifn is null\n");
        exit(1);
    }
    
    if(in._kptr==NULL){
        printf("WARNING cannot copy chisq_wrapper kptr is null\n");
        exit(1);
    }
    
    printf("copying chisq\n");
    
    if(_kptr!=NULL){
        delete _kptr;
        _kptr=NULL;
    }
    _chifn=in._chifn;
    _chimin=in._chimin;
    _target=in._target;
    _deltachi=in._deltachi;
    _ddmin=in._ddmin;
    _adaptive_target=in._adaptive_target;
    _seed=in._seed;
    _mindex=in._mindex;
    
    
    int i,j;
    
    _range_max.reset();
    _range_min.reset();
    _characteristic_length.reset();
    for(i=0;i<in._characteristic_length.get_dim();i++){
        _characteristic_length.set(i,in._characteristic_length.get_data(i));
    }
    for(i=0;i<in._range_min.get_dim();i++){
        _range_min.set(i,in._range_min.get_data(i));
    }
    for(i=0;i<in._range_max.get_dim();i++){
        _range_max.set(i,in._range_max.get_data(i));
    }
    
    _fn.reset();
    
    for(i=0;i<in.get_pts();i++){
        _fn.set(i,in.get_fn(i));
    }
    
    printf("making kptr\n");
    _kptr=new kd_tree(in._kptr[0]);
    if(_seed<0)_seed=int(time(NULL));
    _dice=new Ran(_seed);
    printf("done copying\n");
}
