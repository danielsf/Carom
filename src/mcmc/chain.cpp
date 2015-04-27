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
    _iteration=0;
    _chain_label=0;
    _n_written=0;
    _current_chi=2.0*exception_value;
    _current_degeneracy=0;
    _output_name_root[0]=0;
    _points.set_name("chain_points");
    _degeneracy.set_name("chain_degeneracy");
    _chisquared.set_name("chain_chisquared");
    _current_point.set_name("chain_current_point");
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

void chain::set_output_name_root(char *nn){
    int i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _output_name_root[i]=nn[i];
    }
    _output_name_root[i]=0;
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
}

void chain::read_chain(char *input_name){
    
    printf("reading %s\n",input_name);
    
    int i,ct;
    double chisq,mu;
    array_1d<double> vv;
    vv.set_name("chain_constructor_vv");
    
    FILE *input;
    input=fopen(input_name,"r");
    while(fscanf(input,"%d",&ct)>0){
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

int chain::is_current_point_valid(){
    if(_current_point.get_dim()==_dim){
        return 1;
    }
    
    return 0;
}

double chain::get_current_point(int dex){
    verify_dim(dex,"get_current_point");
    
    return _current_point.get_data(dex);
}

double chain::get_current_chisquared(){
    return _current_chi;
}

int chain::get_current_degeneracy(){
    return _current_degeneracy;
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
    
    int i,add_it;
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
        _current_degeneracy=1;
        for(i=0;i<_dim;i++){
            _current_point.set(i,pt.get_data(i));
        }
    }
    else{
        if(_degeneracy.get_dim()>0){
            _degeneracy.add_val(_degeneracy.get_dim()-1,1);
        }
        else{
            //This means that we recently wrote out the chain and all of the
            //storage arrays are empty

            _points.add_row(pt);
            _chisquared.add(mu);
            _degeneracy.add(1);
        }
        
        _current_degeneracy++;
    }
}

void chain::increment_iteration(){
    _iteration++;
    _n_written=0;
}

void chain::set_chain_label(int ii){
    _chain_label=ii;
}

void chain::write_chain(){
    if(_output_name_root[0]==0){
        printf("WARNING asked to write chain but have no name\n");
        exit(1);
    }

    char output_name[2*letters];
    sprintf(output_name,"%s_%d_%d.txt",_output_name_root,_iteration,_chain_label);

    if(_n_written==0){
        write(output_name,0);
    }
    else{
        write(output_name,1);
    }
    
    _n_written+=_points.get_rows();

    _points.reset_preserving_room();
    _degeneracy.reset_preserving_room();
    _chisquared.reset_preserving_room();

}    
    

void chain::write(char *name, int append){  
    
    FILE *output;
    if(append==0){
        output=fopen(name,"w");
    }
    else{
        output=fopen(name,"a");
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
    _current_degeneracy=in._current_degeneracy;
    _n_written=in._n_written;
    _iteration=in._iteration;
    _chain_label=in._chain_label;
    
    _dim=in._dim;
    _n_written=in._n_written;
    int i,j;
    for(i=0;i<letters;i++){
        _output_name_root[i]=in._output_name_root[i];
    }
    
    _current_point.reset();
    for(i=0;i<in._current_point.get_dim();i++){
        _current_point.set(i,in._current_point.get_data(i));
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
    get_thinned_indices(thinby, burnin, output, -1);
}

void chain::get_thinned_indices(int thinby, int burnin, array_1d<int> &output, int limit){
    output.reset_preserving_room();
    int ct,i,currentDegen,ctStart,total;
    
    if(output.get_dim()!=0){
        printf("WARNING failed to reset output\n");
        exit(1);
    }
    
    ct=0;
    for(i=0;i<_degeneracy.get_dim() && ct+_degeneracy.get_data(i)<burnin;i++){
        ct+=_degeneracy.get_data(i);
    }
    
    total=ct;
       
    ctStart=ct;
    ct=thinby-(burnin-ctStart);
    ctStart=ct;
     
    for(;i<_points.get_rows() && (limit<=0 || total<limit);i++){
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
        total+=_degeneracy.get_data(i);
    }
}

int chain::get_thinby(double threshold, int burnin, int step){
    return get_thinby(threshold, burnin, step, -1);
}

int chain::get_thinby(double threshold, int burnin, int step, int limit){
    ///find the amount to thinby to achieve the specified threshold
    ///burnin is the number of points to discard
    
    array_1d<double> means,covars,vars;
    int i,j,ix,total;
    means.set_name("chain_get_thinby_means");
    covars.set_name("chain_get_thinby_covars");
    vars.set_name("chain_get_thinby_vars");
    
    array_1d<int> dexes;
    dexes.set_name("chain_get_thinby_dexes");

    total=0;
    for(i=0;i<_degeneracy.get_dim();i++){
        total+=_degeneracy.get_data(i);
    }
   
    if(limit>0 && limit<total){
        total=limit;
    } 
    
    int thinby,thinbyBest,i1,i2,bestPts,maxDex,bestDex,repeats;
    double covarMax,covarMaxBest,mu,varMax;
    
    thinbyBest=-1;
    covarMaxBest=2.0*exception_value;
    bestPts=0;
    
    means.set_dim(_dim);
    vars.set_dim(_dim);
    covars.set_dim(_dim);
    
    printf("in get thinby total %d\n",total);

    for(thinby=step;(fabs(threshold-covarMaxBest)>0.1*threshold && covarMaxBest>threshold)
                   && thinby<(total-burnin)/10; thinby+=step){

       get_thinned_indices(thinby,burnin,dexes,limit);
       means.zero();
       for(i=0;i<dexes.get_dim();i++){
           for(j=0;j<_dim;j++){
               means.add_val(j,_points.get_data(dexes.get_data(i),j));
           }
       }
       
       for(i=0;i<_dim;i++){
           means.divide_val(i,double(dexes.get_dim()));
       }
       
       vars.zero();
       for(i=0;i<dexes.get_dim();i++){
           for(j=0;j<_dim;j++){
               vars.add_val(j,power(means.get_data(j)-_points.get_data(dexes.get_data(i),j),2));
           }
       }
       
       for(i=0;i<_dim;i++){
           vars.divide_val(i,double(dexes.get_dim()-1));
       }
       
       repeats=0;
       covars.zero();
       for(i=0;i<dexes.get_dim()-1;i++){
           i1=dexes.get_data(i);
           i2=dexes.get_data(i+1);
           if(i1==i2)repeats++;
           for(j=0;j<_dim;j++){
               covars.add_val(j,(means.get_data(j)-_points.get_data(i1,j))*(means.get_data(j)-_points.get_data(i2,j)));
           }
       }
       
       covarMax=-1.0;
       for(i=0;i<_dim;i++){
           if(dexes.get_dim()>2){
               covars.divide_val(i,double(dexes.get_dim()-2));
           }
           
           covars.divide_val(i,vars.get_data(i));
           if(!isnan(covars.get_data(i))){
               mu=fabs(covars.get_data(i));
               if(mu>covarMax){
                   varMax=sqrt(fabs(vars.get_data(i)))/means.get_data(i);
                   covarMax=mu;
                   maxDex=i;
               }
           }
       }
       
       //printf("    thinby %d covar %e -- %d %d %e %d\n",thinby,covarMax,dexes.get_dim(),repeats,varMax,maxDex);
       
       if(thinbyBest<0 || covarMax<covarMaxBest){
           thinbyBest=thinby;
           covarMaxBest=covarMax;
           bestPts=dexes.get_dim();
           bestDex=maxDex;
       }
       
    }
    
    printf("thinby %d best %e pts %d dex %d\n",thinbyBest,covarMaxBest,bestPts,bestDex);
    
    if(thinbyBest<0){
        thinbyBest = (total-burnin)/10;
    }
    
    return thinbyBest;
    
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
    
    _independent_sample_dexes.set_name("arrayOfChains_independent_sample_dexes");
    _independent_samples.set_name("arrayOfChains_independent_samples");
    
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

int arrayOfChains::get_points(){
    int i,ans;
    ans=0;
    for(i=0;i<_n_chains;i++){
        ans+=_data[i].get_points();
    }
    return ans;
}

int arrayOfChains::get_thinby(double threshold, double burninDenom){

    int i,thinby,thinbyMax,burnin;
    int step,total;
    
    thinbyMax=-1;
    for(i=0;i<_n_chains;i++){
        if(burninDenom>1){
            burnin=_data[i].get_points()/burninDenom;
        }
        else{
            burnin=0;
        }
        
        total=_data[i].get_points();
        if(total-burnin<=100){
            step=(total-burnin)/20;
        }
        else{
            step=10;
        }
        
        if(step==0)step=10;
        
        thinby=_data[i].get_thinby(threshold,burnin,step);
        
        if(thinby>thinbyMax){
            thinbyMax=thinby;
        }
    }
    
    return thinbyMax;
}

void arrayOfChains::get_covariance_matrix(double threshold, int burninDenom, array_2d<double> &covar){

    int i,j,thinbyMax,burnin;
    int step;
    
    thinbyMax=get_thinby(threshold,burninDenom);
 
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
    means.set_name("arrayOfChains_get_covar_means");
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

void arrayOfChains::get_independent_samples(double threshold, int limit){
    
    _independent_sample_dexes.reset();
    _independent_samples.reset();
    
    printf("getting independent samples\n");
    
    int ic;
    int thinby,thinbyMax;
    int total, step;
    
    thinbyMax=-1;
    for(ic=0;ic<_n_chains;ic++){
        total=_data[ic].get_points();
        if(total<=100){
           step=total/20;
        }
        else{
            step=10;
        }
        
        if(step==0){
            step=10;
        }
    
        thinby=_data[ic].get_thinby(threshold,0,step,limit);
        if(thinby>thinbyMax){
            thinbyMax=thinby;
        }
    }
    
    printf("thinning by %d\n",thinbyMax);

    total=0;    
    array_1d<int> temp_dexes;
    temp_dexes.set_name("arrayOfChains_temp_dexes");
    for(ic=0;ic<_n_chains;ic++){
        _data[ic].get_thinned_indices(thinbyMax, 0, temp_dexes, limit);
        _independent_sample_dexes.add_row(temp_dexes);
        total+=temp_dexes.get_dim();
    }
    
    printf("%d independent samples\n",total);
}

void arrayOfChains::_get_full_independent_samples(){
    if(_independent_sample_dexes.get_rows()==0){
        printf("WARNING cannot get full independent samples; there are no dexes\n");
        exit(1);
    }
    
    _independent_samples.set_cols(_dim);
    
    int ic,ip,ix,row,dex;
    row=0;
    for(ic=0;ic<_n_chains;ic++){
        for(ip=0;ip<_independent_sample_dexes.get_cols(ic);ip++){
            dex=_independent_sample_dexes.get_data(ic,ip);
            for(ix=0;ix<_dim;ix++){
                _independent_samples.set(row,ix,_data[ic].get_point(dex,ix));
            }
            row++;
        }
    }
    
    printf("set row %d\n",row);
    _density.set_data(&_independent_samples);
}

void arrayOfChains::calculate_R(array_1d<double> &R, array_1d<double> &V, array_1d<double> &W){
    if(_independent_sample_dexes.get_rows()==0){
        printf("WARNING cannot calculate R; you have no independent samples\n");
        exit(1);
    }

    R.reset();
    V.reset();
    W.reset();
    
    R.set_dim(_dim);
    V.set_dim(_dim);
    W.set_dim(_dim);

    array_1d<double> BoverN,totalMean;
    BoverN.set_name("arrayOfChains_calculate_R_BoverN");
    totalMean.set_name("arrayOfChains_calculate_R_totalMean");
    
    array_2d<double> chainMean;
    chainMean.set_name("arrayOfChains_calculate_R_chainMean");
    
    int ix,ic,ip;
    for(ix=0;ix<_dim;ix++){
        BoverN.set(ix,0.0);
        totalMean.set(ix,0.0);
        W.set(ix,0.0);
    }
    
    chainMean.set_dim(_n_chains,_dim);
    for(ic=0;ic<_n_chains;ic++){
        for(ix=0;ix<_dim;ix++){
            chainMean.set(ic,ix,0.0);
        }
    }
    
    int totalPts,dex;
    
    
    totalPts=0;
    
    for(ic=0;ic<_n_chains;ic++){
        for(ip=0;ip<_independent_sample_dexes.get_cols(ic);ip++){
            totalPts++;
            dex=_independent_sample_dexes.get_data(ic,ip);
            for(ix=0;ix<_dim;ix++){
                chainMean.add_val(ic,ix,_data[ic].get_point(dex,ix));
                totalMean.add_val(ix,_data[ic].get_point(dex,ix));
            }
        }
        for(ix=0;ix<_dim;ix++){
            chainMean.divide_val(ic,ix,double(_independent_sample_dexes.get_cols(ic)));
        }
    }
    
    for(ix=0;ix<_dim;ix++){
        totalMean.divide_val(ix,double(totalPts));
    }
    
    for(ix=0;ix<_dim;ix++){
        for(ic=0;ic<_n_chains;ic++){
            BoverN.add_val(ix,power(chainMean.get_data(ic,ix)-totalMean.get_data(ix),2));
        }
        BoverN.divide_val(ix,double(_n_chains-1));
    }
    
    double mu;
    for(ix=0;ix<_dim;ix++){
        for(ic=0;ic<_n_chains;ic++){
            mu=0.0;
            for(ip=0;ip<_independent_sample_dexes.get_cols(ic);ip++){
                dex=_independent_sample_dexes.get_data(ic,ip);
                mu+=power(_data[ic].get_point(dex,ix)-chainMean.get_data(ic,ix),2);
            }
            mu=mu/double(_independent_sample_dexes.get_cols(ic)-1);
            W.add_val(ix,mu);
        }
        W.divide_val(ix,double(_n_chains));
    }
    
    double nn;
    nn=0.0;
    for(ic=0;ic<_n_chains;ic++){
        nn+=double(_independent_sample_dexes.get_cols(ic));
    }
    nn=nn/double(_n_chains);
    
    double sigmaPlus;
    for(ix=0;ix<_dim;ix++){
        sigmaPlus=(nn-1.0)*W.get_data(ix)/nn + BoverN.get_data(ix);
        V.set(ix,sigmaPlus+BoverN.get_data(ix)/double(_n_chains));
        R.set(ix,V.get_data(ix)/W.get_data(ix));
    }

}

void arrayOfChains::plot_contours(int ix, int iy, double fraction, char *nameRoot){
    if(_independent_sample_dexes.get_rows()==0){
        printf("WARNING cannot plot contours; no independent samples\n");
        exit(1);
    }
    
    if(_independent_samples.get_rows()==0){
        _get_full_independent_samples();
    }
    
    double dx,dy;
    double xmax,xmin,ymax,ymin;
    int i;
    
    for(i=0;i<_independent_samples.get_rows();i++){
        if(i==0 || _independent_samples.get_data(i,ix)<xmin){
            xmin=_independent_samples.get_data(i,ix);
        }
        
        if(i==0 || _independent_samples.get_data(i,ix)>xmax){
            xmax=_independent_samples.get_data(i,ix);
        }
        
        if(i==0 || _independent_samples.get_data(i,iy)<ymin){
            ymin=_independent_samples.get_data(i,iy);
        }
        
        if(i==0 || _independent_samples.get_data(i,iy)>ymax){
            ymax=_independent_samples.get_data(i,iy);
        }
    }
    
    dx=(xmax-xmin)/100.0;
    dy=(ymax-ymin)/100.0;
    
    char boundaryName[2*letters],scatterName[2*letters];
    sprintf(boundaryName,"%s_%d_%d_contour.txt",nameRoot,ix,iy);
    sprintf(scatterName,"%s_%d_%d_scatter.txt",nameRoot,ix,iy);
    
    
    
    _density.plot_density(ix,dx,iy,dy,fraction,scatterName,3);
    _density.plot_boundary(ix,dx,iy,dy,fraction,boundaryName,3);
    
}

int arrayOfChains::get_n_samples(){
    return _independent_samples.get_rows();
}

double arrayOfChains::get_sample(int dex, int ix){
    return _independent_samples.get_data(dex,ix);
}

double arrayOfChains::acceptance_rate(){
    int ct,rows;
    int i,j;
    
    ct=0;
    rows=0;
    for(i=0;i<_n_chains;i++){
        for(j=0;j<_data[i].get_rows();j++){
            rows++;
            ct+=_data[i].get_degeneracy(j);
        }
    }
    
    return double(rows)/double(ct);
}
