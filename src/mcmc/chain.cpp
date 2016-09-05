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
    _total_written=0;
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
    _total_written=0;
}

void chain::set_chain_label(int ii){
    _chain_label=ii;
}

void chain::write_chain(int dump){
    if(_output_name_root[0]==0){
        printf("WARNING asked to write chain but have no name\n");
        exit(1);
    }

    char output_name[2*letters];
    sprintf(output_name,"%s_%d_%d.txt",_output_name_root,_iteration,_chain_label);

    if(dump==0){
        if(_degeneracy.get_data(_degeneracy.get_dim()-1)>1){
            _degeneracy.subtract_val(_degeneracy.get_dim()-1,1);
        }
        else{
            _degeneracy.remove(_degeneracy.get_dim()-1);
            _chisquared.remove(_chisquared.get_dim()-1);
            _points.remove_row(_points.get_rows()-1);
        }
    }

    if(_total_written==0){
        write(output_name,0);
    }
    else{
        write(output_name,1);
    }

    _total_written+=_points.get_rows();

    int zeropt;

    if(dump==1){
        _n_written=0;
        _points.reset_preserving_room();
        _degeneracy.reset_preserving_room();
        _chisquared.reset_preserving_room();
    }
    else{
        zeropt=_n_written;
        _n_written+=_points.get_rows()-zeropt;
        _points.add_row(_current_point);
        _chisquared.add(_current_chi);
        _degeneracy.add(1);
    }
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
    for(i=_n_written;i<_points.get_rows();i++){
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
    _total_written=in._total_written;
    _n_written=in._n_written;
    _iteration=in._iteration;
    _chain_label=in._chain_label;
    _dim=in._dim;

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

    total=0;
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

    if(get_points()==0){
        return -1;
    }

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

    if(limit>0 && limit+burnin<total){
        total=limit+burnin;
    }

    int thinby,thinbyBest,i1,i2,bestPts,maxDex,bestDex,repeats;
    double covarMax,covarMaxBest,mu,varMax,covRat,meanRat,meanVal;

    thinbyBest=-1;
    covarMaxBest=2.0*exception_value;
    bestPts=0;

    means.set_dim(_dim);
    vars.set_dim(_dim);
    covars.set_dim(_dim);

    //printf("in get thinby total %d\n",total);

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

           mu=sqrt(fabs(covars.get_data(i)))/fabs(means.get_data(i));
           //mu=fabs(covars.get_data(i))/fabs(vars.get_data(i));
           if(mu>covarMax && !isnan(mu)){
               varMax=sqrt(fabs(vars.get_data(i)))/means.get_data(i);
               covRat=covars.get_data(i)/vars.get_data(i);
               meanRat=sqrt(fabs(covars.get_data(i)))/fabs(means.get_data(i));
               meanVal=means.get_data(i);
               covarMax=mu;
               maxDex=i;
           }

       }

       //printf("    thinby %d covar %.3e %.3e -- %d %d %.3e %.3e %.3e\n",thinby,covarMax,meanVal,
       //dexes.get_dim(),repeats,varMax,covRat,meanRat);

       if(thinbyBest<0 || covarMax<covarMaxBest){
           thinbyBest=thinby;
           covarMaxBest=covarMax;
           bestPts=dexes.get_dim();
           bestDex=maxDex;
       }

    }

    //printf("thinby %d best %e pts %d dex %d\n",thinbyBest,covarMaxBest,bestPts,bestDex);

    if(thinbyBest<=0){
        thinbyBest = (total-burnin)/10;
    }
    //printf("    thinbyBest %d covarMaxBest %e\n",thinbyBest,covarMaxBest);
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
    _chisq_guess=-1.0;

}

void arrayOfChains::initialize(int nChains, int dim, Ran *dice){
    _dim=dim;
    _n_chains=nChains;
    _dice=dice;

    _data = new chain[_n_chains];

    _independent_sample_dexes.set_name("arrayOfChains_independent_sample_dexes");
    _independent_samples.set_name("arrayOfChains_independent_samples");
    _independent_sample_weights.set_name("arrayOfChains_independent_sample_weights");

    _contour_maxes.set_name("arrayOfChains_contour_maxes");
    _contour_mins.set_name("arrayOfChains_contour_mins");

    int i;
    for(i=0;i<_n_chains;i++){
        _data[i].set_dim(_dim);
        _data[i].set_dice(_dice);
    }

    _chisq_guess=-1.0;

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

void arrayOfChains::get_covariance_matrix(array_2d<double> &covar){

    int iChain,iRow,total_points,i;

    array_1d<double> means;
    means.set_name("arrayOfChains_get_covar_means");
    means.set_dim(_dim);
    means.zero();

    total_points=0;
    for(iChain=0;iChain<_n_chains;iChain++){
        for(iRow=0;iRow<_data[iChain].get_rows();iRow++){
            for(i=0;i<_dim;i++){
                means.add_val(i,_data[iChain].get_degeneracy(iRow)*_data[iChain].get_point(iRow,i));
            }
            total_points+=_data[iChain].get_degeneracy(iRow);
        }
    }

    for(i=0;i<_dim;i++){
        means.divide_val(i,double(total_points));
    }


    covar.reset();
    covar.set_dim(_dim,_dim);
    covar.zero();

    int ix,iy;
    for(iChain=0;iChain<_n_chains;iChain++){
        for(iRow=0;iRow<_data[iChain].get_rows();iRow++){
            for(ix=0;ix<_dim;ix++){
                for(iy=ix;iy<_dim;iy++){
                    covar.add_val(ix,iy,(means.get_data(ix)-_data[iChain].get_point(iRow,ix))*
                                        (means.get_data(iy)-_data[iChain].get_point(iRow,iy))*
                                        _data[iChain].get_degeneracy(iRow));
                }
            }
        }
    }


    for(ix=0;ix<_dim;ix++){
        for(iy=ix;iy<_dim;iy++){
            covar.divide_val(ix,iy,double(total_points-1));
            if(ix!=iy){
                covar.set(iy,ix,covar.get_data(ix,iy));
            }
        }
    }

}

void arrayOfChains::use_all(int burnin, int limit){
    _independent_sample_dexes.reset();
    _independent_samples.reset();
    _independent_sample_weights.reset();

    int i,total,burned,j;
    int iChain,toburn,touse;

    array_1d<int> temp_dexes;
    temp_dexes.set_name("arrayOfChains_use_all_temp_dexes");

    printf("burn %d lim %d\n",burnin,limit);

    for(iChain=0;iChain<_n_chains;iChain++){
        temp_dexes.reset_preserving_room();
        total=0;
        burned=0;
        for(i=0;i<_data[iChain].get_rows() && burned<burnin;i++){
           // printf("burning %d burned %d %d %d\n",_data[iChain].get_degeneracy(i),burned,i,_data[iChain].get_rows());
            if(burned+_data[iChain].get_degeneracy(i)<burnin){
                burned+=_data[iChain].get_degeneracy(i);
            }
            else{
                toburn=burnin-burned;
                touse=_data[iChain].get_degeneracy(i)-toburn;
                temp_dexes.add(i);
                _independent_sample_weights.add(double(touse));
                burned=burnin;
                total=touse;
            }

            //printf("burned %d\n",burned);
        }

        //printf("total %d %d %d\n",total,i,_data[iChain].get_rows());

        for(;i<_data[iChain].get_rows() && (limit<=0 || total<limit);i++){

            total+=_data[iChain].get_degeneracy(i);
            temp_dexes.add(i);
            _independent_sample_weights.add(double(_data[iChain].get_degeneracy(i)));
        }

        _independent_sample_dexes.add_row(temp_dexes);

    }

    printf("ind rows %d\n",_independent_sample_dexes.get_rows());
}

void arrayOfChains::get_independent_samples(double threshold, int burnin, int limit){

    _independent_sample_dexes.reset();
    _independent_samples.reset();
    _independent_sample_weights.reset();

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

        thinby=_data[ic].get_thinby(threshold,burnin,step,limit);
        if(thinby>thinbyMax){
            thinbyMax=thinby;
        }
    }

    printf("thinning by %d\n",thinbyMax);

    int iw;

    total=0;
    array_1d<int> temp_dexes;
    temp_dexes.set_name("arrayOfChains_temp_dexes");
    for(ic=0;ic<_n_chains;ic++){
        _data[ic].get_thinned_indices(thinbyMax, burnin, temp_dexes, limit);
        _independent_sample_dexes.add_row(temp_dexes);
        for(iw=0;iw<temp_dexes.get_dim();iw++){
            _independent_sample_weights.add(1.0);
        }
        total+=temp_dexes.get_dim();
    }

    array_1d<double> means,vars,covars;
    means.set_name("arrayOfChains_get_independent_samples_means");
    vars.set_name("arrayOfChains_get_independent_samples_vars");
    covars.set_name("arrayOfChains_get_independent_samples_covars");

    means.set_dim(_dim);
    vars.set_dim(_dim);
    covars.set_dim(_dim);

    int i,j,dex;
    for(ic=0;ic<_independent_sample_dexes.get_rows();ic++){
        for(i=0;i<_independent_sample_dexes.get_cols(ic);i++){
            dex=_independent_sample_dexes.get_data(ic,i);
            for(j=0;j<_dim;j++){
                means.add_val(j,_data[ic].get_point(dex,j));
            }
        }
    }

    for(i=0;i<_dim;i++){
        means.divide_val(i,double(total));
    }

    int dex2,covarTotal;

    covarTotal=0;
    for(ic=0;ic<_independent_sample_dexes.get_rows();ic++){
        for(i=0;i<_independent_sample_dexes.get_cols(ic);i++){
            dex=_independent_sample_dexes.get_data(ic,i);
            if(i<_independent_sample_dexes.get_cols(ic)-1){
                dex2=_independent_sample_dexes.get_data(ic,i+1);
                covarTotal++;
            }
            for(j=0;j<_dim;j++){
                vars.add_val(j,power(_data[ic].get_point(dex,j)-means.get_data(j),2));
                if(i<_independent_sample_dexes.get_cols(ic)-1){
                    covars.add_val(j,(_data[ic].get_point(dex,j)-means.get_data(j))*
                                     (_data[ic].get_point(dex2,j)-means.get_data(j)));
                }
            }
        }
    }

    double covarMax=-1.0;
    for(i=0;i<_dim;i++){
        vars.divide_val(i,double(total-1));
        covars.divide_val(i,double(covarTotal-1));
        covars.divide_val(i,vars.get_data(i));
        if(fabs(covars.get_data(i))>covarMax){
            covarMax=fabs(covars.get_data(i));
        }
    }


    printf("%d independent samples -- covar in time %e\n",total,covarMax);
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
    _density.set_data(&_independent_samples, _independent_sample_weights);
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

void arrayOfChains::plot_chisquared_histogram(int limit, double min, double max, double dx, char *nameRoot){
    int total,globalTotal;
    int iChain,i;
    double cc,mu;
    int dex;

    array_1d<int> counts;

    for(cc=min;cc<=max;cc+=dx){
        counts.add(0);
    }

    printf("in histogram limit %d\n",limit);

    globalTotal=0;
    for(iChain=0;iChain<_n_chains;iChain++){
        total=0;
        for(i=0;i<_data[iChain].get_rows() && (limit<=0 || total<limit);i++){
            total+=_data[iChain].get_degeneracy(i);
            globalTotal+=_data[iChain].get_degeneracy(i);
            cc=_data[iChain].get_chisquared(i);

            mu=(cc-min)/dx;
            dex=int(mu);

            if(mu-1.0*dex>0.5){
                dex++;
            }

            if(dex>counts.get_dim()-1){
                dex=counts.get_dim()-1;
            }

            counts.add_val(dex,_data[iChain].get_degeneracy(i));

        }
    }

    char name[2*letters];
    sprintf(name,"%s_histogram.sav",nameRoot);
    FILE *output;
    output=fopen(name,"w");
    for(i=0;i<counts.get_dim();i++){
        fprintf(output,"%e %d\n",min+dx*i,counts.get_data(i));
    }
    fclose(output);

    printf("plotted histogram of %d points\n",globalTotal);


}

double arrayOfChains::get_chimin(){
    int i,j;
    double chimin=2.0*exception_value;
    for(i=0;i<_n_chains;i++){
        for(j=0;j<_data[i].get_rows();j++){
            if(_data[i].get_chisquared(j)<chimin){
                chimin=_data[i].get_chisquared(j);
            }
        }
    }

    return chimin;
}

void arrayOfChains::_get_contour_bounds(){
    if(_independent_samples.get_rows()==0){
        printf("WARNING cannot get contour bounds; need full independent samples\n");
        exit(1);
    }
    array_1d<double> chisq_list,chisq_list_sorted;
    array_1d<int> dexes;

    _contour_mins.reset();
    _contour_maxes.reset();

    chisq_list.set_name("arrayOfChains_get_chisq_guess_chisq_list");
    chisq_list_sorted.set_name("arrayOfChains_get_chisq_guess_sorted");
    dexes.set_name("arrayOfChains_get_chisq_guess_dexes");

    int i,j,dex;

    for(i=0;i<_dim;i++){
        _contour_maxes.set(i,-2.0*exception_value);
        _contour_mins.set(i,2.0*exception_value);
    }

    double total,mu;

    total=0.0;
    dex=0;
    for(i=0;i<_independent_sample_dexes.get_rows();i++){
        for(j=0;j<_independent_sample_dexes.get_cols(i);j++){
            dexes.add(dex);
            mu=_data[i].get_chisquared(_independent_sample_dexes.get_data(i,j));
            total+=exp(-0.5*mu)*_independent_sample_weights.get_data(dex);
            chisq_list.add(mu);
            dex++;
        }
    }

    sort(chisq_list,chisq_list_sorted,dexes);
    double sum;
    int ix;
    sum=0.0;
    for(i=0;i<chisq_list_sorted.get_dim() && sum<0.68*total;i++){
        sum+=exp(-0.5*chisq_list_sorted.get_data(i))*_independent_sample_weights.get_data(dexes.get_data(i));

        dex=dexes.get_data(i);
        for(ix=0;ix<_dim;ix++){
            if(_independent_samples.get_data(dex,ix)<_contour_mins.get_data(ix)){
                _contour_mins.set(ix,_independent_samples.get_data(dex,ix));
            }

            if(_independent_samples.get_data(dex,ix)>_contour_maxes.get_data(ix)){
                _contour_maxes.set(ix,_independent_samples.get_data(dex,ix));
            }
        }
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

    if(_contour_maxes.get_dim()==0){
        _get_contour_bounds();
    }

    double dx,dy;
    double xmax,xmin,ymax,ymin;
    int i,j,dex;



    dx=(_contour_maxes.get_data(ix)-_contour_mins.get_data(ix))/50.0;
    dy=(_contour_maxes.get_data(iy)-_contour_mins.get_data(iy))/50.0;

    printf("printing %d %e %d %e\n",ix,dx,iy,dy);

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

void arrayOfChains::acceptance_statistics(int burnin, int limit){
    int iChain,iRow,degen;
    double mean,var;
    int max,nOver5,nOver10;
    int total,burned,startRow,rowCt;

    for(iChain=0;iChain<_n_chains;iChain++){
        total=0;
        burned=0;
        mean=0.0;
        var=0.0;
        max=0;
        nOver5=0;
        nOver10=0;
        rowCt=0;

        for(iRow=0;iRow<_data[iChain].get_rows() && burned<burnin; iRow++){
            burned+=_data[iChain].get_degeneracy(iRow);
        }

        startRow=iRow;

        for(;iRow<_data[iChain].get_rows() && (limit<=0 || total<limit); iRow++){
            total+=_data[iChain].get_degeneracy(iRow);
            degen=_data[iChain].get_degeneracy(iRow);
            mean+=double(degen);
            if(degen>max){
                max=degen;
            }

            if(degen>5){
                nOver5++;
            }

            if(degen>10){
                nOver10++;
            }

            rowCt++;

        }

        mean=mean/double(rowCt);

        total=0;
        rowCt=0;
        for(iRow=startRow;iRow<_data[iChain].get_rows() && (limit<=0 || total<limit); iRow++){
            degen=_data[iChain].get_degeneracy(iRow);
            var+=power(double(degen)-mean,2);
            rowCt++;
        }
        var=var/double(rowCt-1);

        printf("chain %d\n",iChain);
        printf("mean %e var %e sqrt(var) %e\n",mean,var,sqrt(var));
        printf("max %d nOver5 %d nOver10 %d\n\n",max,nOver5,nOver10);

    }

}
