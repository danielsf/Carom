#include "mcmc/chain.h"

chain::~chain(){}

void chain::is_dice_safe(char *routine){
    if(_dice==NULL){
        printf("WARNING in chain::%s\n",routine);
        printf("_dice is null\n");
        exit(1);
    }
}

void chain::initialize(){
    _dice=NULL;
    _dim=0;
    _n_written=0;
    _current_chi=2.0*exception_value;
    _output_name[0]=0;
    _points.set_name("chain_points");
    _degeneracy.set_name("chain_degeneracy");
    _chisquared.set_name("chain_chisquared");
}

void chain::set_dim(int ii){
    _dim=ii;
    _points.reset();
    _degeneracy.reset();
    _chisquared.reset();
    
    _points.set_cols(_dim);
}

void chain::set_dice(Ran *dd){
   _dice=dd;
}

void chain::set_output_name(char *nn){
    int i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _output_name[i]=nn[i];
    }
    _output_name[i]=0;
}

chain::chain(){
    initialize();    
}

chain::chain(int ii){
    initialize();
    set_dim(ii);
}

chain::chain(int ii, char *input_name){
    initialize();
    set_dim(ii);
    
    int i,ct;
    double chisq,mu;
    array_1d<double> vv;
    vv.set_name("chain_constructor_vv");
    
    FILE *input;
    input=fopen(input_name,"r");
    while(fscanf(input,"%d",&ct)){
        _degeneracy.add(ct);
        fscanf(input,"%le",&chisq);
        _chisquared.add(chisq);
        
        for(i=0;i<_dim;i++){
            fscanf(input,"%le",&mu);
            vv.set(i,mu);
        }
        _points.add_row(vv);
    }
    
    fclose(input);
}

void chain::verify_dim(int dex, char *routine){
    if(dex<0 || dex>=_dim){
        printf("WARNING chain::%s\n",routine);
        printf("%d but dim %d\n",dex,_dim);
        exit(1);
    }
}

void chain::verify_points(int dex, char *routine){
    if(dex<0 || dex>=_points.get_rows()){
        printf("WARNING chain::%s\n",routine);
        printf("%d but points %d\n",dex,_points.get_rows());
        exit(1);
    }
}

double chain::get_current_point(int dex){
    verify_dim(dex,"get_current_point");
    
    return _points.get_data(_points.get_rows()-1,dex);
}

double chain::get_current_chisquared(){
    return _chisquared.get_data(_chisquared.get_dim()-1);
}

int chain::get_current_degeneracy(){
    return _degeneracy.get_data(_degeneracy.get_dim()-1);
}

int chain::get_points(){
    int i,ans;
    ans=0;
    for(i=0;i<_degeneracy.get_dim();i++){
        ans+=_degeneracy.get_data(i);
    }
    return ans;
}

int chain::get_rows(){
    return _points.get_rows();
}

int chain::get_dim(){
    return _dim;
}

double chain::get_point(int dex, int dim){
    verify_dim(dim,"get_point");
    verify_points(dex,"get_point");
    return _points.get_data(dex,dim);
}

double chain::get_chisquared(int dex){
    verify_points(dex,"get_chisquared");
    return _chisquared.get_data(dex);
}

int chain::get_degeneracy(int dex){
    verify_points(dex,"get_degeneracy");
    return _degeneracy.get_data(dex);
}

void chain::add_point(array_1d<double> &pt, double mu){
    is_dice_safe("add_point");
    
    int add_it;
    double roll,ratio;
    
    add_it=0;

    if(mu<_current_chi){
        add_it=1;
    }
    else{
        roll=_dice->doub();
        ratio=exp(-0.5*(mu-_current_chi));
        if(ratio>roll){
            add_it=1;
        }
    }
    
    if(add_it){
        _points.add_row(pt);
        _chisquared.add(mu);
        _degeneracy.add(1);
        _current_chi=mu;
    }
    else{
        _degeneracy.add_val(_degeneracy.get_dim()-1,1);
    }
}

void chain::write_chain(){
    if(_output_name[0]==0){
        printf("WARNING asked to write chain but have no name\n");
        exit(1);
    }
    
    FILE *output;
    if(_n_written==0){
        output=fopen(_output_name,"w");
    }
    else{
        output=fopen(_output_name,"a");
    }
    
    int i,j;
    for(i=0;i<_points.get_rows();i++){
        fprintf(output,"%d %e ",_degeneracy.get_data(i),_chisquared.get_data(i));
        for(j=0;j<_dim;j++){
            fprintf(output,"%e ",_points.get_data(i,j));
        }
        fprintf(output,"\n");
    }
    fclose(output);
    
    _n_written+=_points.get_rows();
    
    _points.reset_preserving_room();
    _degeneracy.reset_preserving_room();
    _chisquared.reset_preserving_room();
}

void chain::copy(const chain &in){

    if(this==&in){
        return;
    }

    _points.reset();
    _degeneracy.reset();
    _chisquared.reset();
    _dice=in._dice;
    _current_chi=in._current_chi;
    
    _dim=in._dim;
    _n_written=in._n_written;
    int i,j;
    for(i=0;i<letters;i++){
        _output_name[i]=in._output_name[i];
    }
    
    _points.set_cols(_dim);
    for(i=0;i<in._points.get_rows();i++){
        _degeneracy.set(i,in._degeneracy.get_data(i));
        _chisquared.set(i,in._chisquared.get_data(i));
        for(j=0;j<_dim;j++){
            _points.set(i,j,in._points.get_data(i,j));
        }
    }
}

void chain::get_thinned_indices(int thinby, int burnin, array_1d<int> &output){
    output.reset_preserving_room();
    int ct,i,currentDegen,ctStart;
    
    ct=0;
    for(i=0;i<_degeneracy.get_dim() && ct+_degeneracy.get_data(i)<burnin;i++){
        ct+=_degeneracy.get_data(i);
    }
       
    ctStart=ct;
    ct=thinby-(burnin-ctStart);
    ctStart=ct;
       
    for(;i<_points.get_rows();i++){
        if(ct+_degeneracy.get_data(i)>=thinby){
            currentDegen=_degeneracy.get_data(i);
            while(ct+currentDegen>=thinby){
                output.add(i);
                currentDegen-=(thinby-ct);
                ct=0;
            }
            ct=currentDegen;
        }
        else{
            ct+=_degeneracy.get_data(i);
        }
    }
}

int chain::get_thinby(double threshold, int burnin, int step){
    ///find the amount to thinby to achieve the specified threshold
    ///burnin is the number of points to discard
    
    array_1d<double> means,covars;
    int i,j,ix,total;
    means.set_name("chain_get_thinby_means");
    covars.set_name("chain_get_thinby_covars");
    
    array_1d<int> dexes;
    dexes.set_name("chain_get_thinby_dexes");
    
    total=0;
    for(i=0;i<_degeneracy.get_dim();i++){
        total+=_degeneracy.get_data(i);
    }
    
    int thinby,thinbyBest,i1,i2;
    double covarMax,covarMaxBest,lastVal;
    
    thinbyBest=-1;
    covarMaxBest=2.0*exception_value;
    
    means.set_dim(_dim);
    
    for(thinby=step;covarMaxBest>threshold && thinby>(total-burnin)/10; thinby+=step){
       get_thinned_indices(thinby,burnin,dexes);
       means.zero();
       for(i=0;i<dexes.get_dim();i++){
           for(j=0;j<_dim;j++){
               means.add_val(j,_points.get_data(dexes.get_data(i),j));
           }
       }
       
       for(i=0;i<_dim;i++){
           means.divide_val(i,double(dexes.get_dim()));
       }
       
       covars.zero();
       for(i=0;i<dexes.get_dim()-1;i++){
           i1=dexes.get_data(i);
           i2=dexes.get_data(i+1);
           for(j=0;j<_dim;j++){
               covars.add_val(j,(means.get_data(j)-_points.get_data(i1,j))*(means.get_data(j)-_points.get_data(i2,j)));
           }
       }
       
       covarMax=-1.0;
       for(i=0;i<_dim;i++){
           if(dexes.get_dim()>2){
               covars.divide_val(i,double(dexes.get_dim()-2));
           }
           
           if(fabs(covars.get_data(i))>covarMax){
               covarMax=fabs(covars.get_data(i));
           }
       }
       
       if(thinbyBest<0 || covarMax<covarMaxBest){
           thinbyBest=thinby;
           covarMaxBest=covarMax;
       }
       
    }
    printf("thinby %d %e\n",thinby,covarMaxBest);
    return thinbyBest;
    
}

void chain::get_thinned_samples(int thinby, int burnin, array_2d<double> &samples){
    samples.reset();
    samples.set_cols(_dim);
    array_1d<int> dexes;
    dexes.set_name("chain_get_thinned_samples_dexes");
    get_thinned_indices(thinby,burnin,dexes);
    int i;
    for(i=0;i<dexes.get_dim();i++){
        samples.add_row(_points(dexes.get_data(i))[0]);
    }
}


////////////////////array of chains

arrayOfChains::~arrayOfChains(){
    if(_data!=NULL){
        delete [] _data;
    }
}

arrayOfChains::arrayOfChains(){
    _data=NULL;
    _n_chains=0;
    _dim=0;
}

void arrayOfChains::initialize(int nChains, int dim, Ran *dice){
    _dim=dim;
    _n_chains=nChains;
    _dice=dice;
    
    _data = new chain[_n_chains];
    
    int i;
    for(i=0;i<_n_chains;i++){
        _data[i].set_dim(_dim);
        _data[i].set_dice(_dice);
    }
    
}

arrayOfChains::arrayOfChains(int nChains, int dim, Ran *dice){
    initialize(nChains,dim,dice);
}

int arrayOfChains::get_n_chains(){
    return _n_chains;
}

arrayOfChains::arrayOfChains(array_2d<double> &pts, array_1d<double> &chisq, Ran *dice){
    initialize(pts.get_rows(),pts.get_cols(),dice);
    
    int i;
    for(i=0;i<_n_chains;i++){
        _data[i].add_point(pts(i)[0],chisq.get_data(i));
    }
}

void arrayOfChains::verify_chains(int dex, char *routine){
    if(dex<0 || dex>=_n_chains){
        printf("WARNING arrayOfChains::%s",routine);
        printf("asked for %d but _n_chains %d\n",dex,_n_chains);
        exit(1);
    }
}

chain* arrayOfChains::operator()(int dex){
    verify_chains(dex,"operator");
    return &_data[dex];
}

void arrayOfChains::add(array_1d<double> &pt, double mu){
    chain *buffer;
    int i;
    
    buffer=NULL;
    
    if(_n_chains>0){
        buffer=new chain[_n_chains];
        for(i=0;i<_n_chains;i++){
            buffer[i].copy(_data[i]);
        }
        delete [] _data;
    }
    else{
        if(_data!=NULL){
            printf("WARNING _n_chains is zero but data is not null\n");
            exit(1);
        }
        
        _dim=pt.get_dim();
    }
    
    _data=new chain[_n_chains+1];
    
    if(buffer!=NULL){
        for(i=0;i<_n_chains;i++){
            _data[i].copy(buffer[i]);
        }
        delete [] buffer;
    }
    
    _data[_n_chains].set_dim(_dim);
    _data[_n_chains].set_dice(_dice);
    _data[_n_chains].add_point(pt,mu);
    _n_chains++;
}

void arrayOfChains::add(array_2d<double> &pts, array_1d<double> &mu){
    chain *buffer;
    int i;
    
    buffer=NULL;
    
    if(_n_chains>0){
        buffer=new chain[_n_chains];
        for(i=0;i<_n_chains;i++){
            buffer[i].copy(_data[i]);
        }
        delete [] _data;
    }
    else{
        if(_data!=NULL){
            printf("WARNING _n_chains is zero but data is not null\n");
            exit(1);
        }
        
        _dim=pts.get_cols();
    }
    
    _data=new chain[_n_chains+pts.get_rows()];
    
    if(buffer!=NULL){
        for(i=0;i<_n_chains;i++){
            _data[i].copy(buffer[i]);
        }
        delete [] buffer;
    }
    
    int j;
    for(i=0;i<pts.get_rows();i++){
        j=_n_chains+i;
        _data[j].set_dim(_dim);
        _data[j].set_dice(_dice);
        _data[j].add_point(pts(i)[0],mu.get_data(i));
    }
    
    _n_chains+=pts.get_rows();
}

void arrayOfChains::remove(int dex){
    verify_chains(dex,"remove");
    chain *buffer;
    buffer=NULL;
    
    int i,j;
    
    if(_n_chains>1){
        buffer=new chain[_n_chains-1];
        for(i=0,j=0;i<_n_chains;i++){
            if(i!=dex){
                buffer[j].copy(_data[i]);
                j++;
            }
        }
    }
    
    delete [] _data;
    _data=NULL;
    
    _n_chains--;
    
    if(buffer!=NULL){
        _data=new chain[_n_chains];
        for(i=0;i<_n_chains;i++){
            _data[i].copy(buffer[i]);
        }
        delete [] buffer;
    }
    
}

void arrayOfChains::get_covariance_matrix(double threshold, int burninDenom, array_2d<double> &covar){

    int i,j,thinby,thinbyMax,burnin;
    
    thinbyMax=-1;
    for(i=0;i<_n_chains;i++){
        if(burninDenom>1){
            burnin=_data[i].get_points()/burninDenom;
        }
        else{
            burnin=0;
        }
        thinby=_data[i].get_thinby(threshold,burnin,10);
        
        if(thinby>thinbyMax){
            thinbyMax=thinby;
        }
    }

    array_1d<int> dexes,tags,temp_dexes;
    dexes.set_name("arrayOfChains_get_covar_dexes");
    tags.set_name("arrayOfChains_get_covar_tags");
    temp_dexes.set_name("arrayOfChains_get_covar_temp_dexes");
    
    for(i=0;i<_n_chains;i++){
        if(burninDenom>1){
            burnin=_data[i].get_points()/burninDenom;
        }
        else{
            burnin=0;
        }
        _data[i].get_thinned_indices(thinbyMax,burnin,temp_dexes);
        for(j=0;j<temp_dexes.get_dim();j++){
            dexes.add(temp_dexes.get_data(j));
            tags.add(i);
        }
    }
    
    array_1d<double> means;
    means.set_name("arrayOfchains_get_covar_means");
    means.set_dim(_dim);
    means.zero();
    for(i=0;i<dexes.get_dim();i++){
        for(j=0;j<_dim;j++){
            means.add_val(j,_data[tags.get_data(i)].get_point(dexes.get_data(i),j));
        }
    }
    
    for(i=0;i<_dim;i++){
        means.divide_val(i,double(dexes.get_dim()));
    }
    
    covar.reset();
    covar.set_dim(_dim,_dim);
    covar.zero();
    
    int iChain,iPt,ix,iy;
    for(i=0;i<dexes.get_dim();i++){
        iChain=tags.get_data(i);
        iPt=dexes.get_data(i);
        for(ix=0;ix<_dim;ix++){
            for(iy=ix;iy<_dim;iy++){
                covar.add_val(ix,iy,(means.get_data(ix)-_data[iChain].get_point(iPt,ix))*
                                    (means.get_data(iy)-_data[iChain].get_point(iPt,iy)));
            
            }
        }
    }
    
    for(ix=0;ix<_dim;ix++){
        for(iy=ix;iy<_dim;iy++){
            covar.divide_val(ix,iy,double(dexes.get_dim()-1));
            if(ix!=iy){
                covar.set(iy,ix,covar.get_data(ix,iy));
            }
        }
    }

}
