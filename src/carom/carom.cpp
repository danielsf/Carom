#include "carom.h"

carom::carom(){
    _write_every=3000;
    _last_written=0;
    sprintf(_outname,"output/carom_output.sav");
    sprintf(_timingname,"output/carom_timing.sav");
    _time_started=double(time(NULL));
}

carom::~carom(){}

void carom::set_seed(int ii){
    _chifn.set_seed(ii);
}

void carom::set_min(array_1d<double> &vv){
    _chifn.set_min(vv);
}

void carom::set_max(array_1d<double> &vv){
    _chifn.set_max(vv);
}

void carom::set_characteristic_length(int dex, double vv){
    _chifn.set_characteristic_length(dex,vv);
}

void carom::set_chisquared(chisquared *xx){
    _chifn.set_chisquared(xx);
}

void carom::set_deltachi(double xx){
    _chifn.set_deltachi(xx);
}

void carom::set_target(double tt){
    _chifn.set_target(tt);
}

void carom::set_write_every(int ww){
    _write_every=ww;
}

void carom::set_outname(char *nn){
    int i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _outname[i]=nn[i];
    }
    _outname[i]=0;
}

void carom::set_timingname(char *nn){
    int i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _timingname[i]=nn[i];
    }
    _timingname[i]=0;
}

void carom::initialize(int npts){
    _chifn.initialize(npts);
    write_pts();
}

void carom::write_pts(){
    FILE *output;
    
    int i,j;
    output=fopen(_outname,"w");
    for(i=0;i<_chifn.get_pts();i++){
        for(j=0;j<_chifn.get_dim();j++){
            fprintf(output,"%.18e ",_chifn.get_pt(i,j));
        }
        fprintf(output,"%.18e 0 0 1\n",_chifn.get_fn(i));
    }
    fclose(output);
    
    output=fopen(_timingname,"a");
    fprintf(output,"%d %e %e\n",
        _chifn.get_called(),
        double(time(NULL))-_time_started,
        (double(time(NULL))-_time_started)/double(_chifn.get_called()));
    fclose(output);
    

    _last_written=_chifn.get_called();
}

void carom::simplex_search(){

    simplex_minimizer ffmin;
    ffmin.set_chisquared(&_chifn);
    ffmin.set_dice(_chifn.get_dice());
    array_1d<double> min,max;
    min.set_name("carom_simplex_search_min");
    max.set_name("carom_simplex_search_min");
    _chifn.get_min(min);
    _chifn.get_max(max);
    ffmin.set_minmax(min,max);
    ffmin.use_gradient();
    
    array_1d<double> minpt;
    minpt.set_name("carom_simplex_search_minpt");
    
    array_2d<double> seed;
    seed.set_name("carom_simplex_search_seed");
    
    seed.set_cols(_chifn.get_dim());
    int i,j;
    for(i=0;i<_chifn.get_dim()+1;i++){
        for(j=0;j<_chifn.get_dim();j++){
            seed.set(i,j,_chifn.get_pt(i,j));
        }
    }
    
    ffmin.find_minimum(seed,minpt);
    write_pts();
}

void carom::search(){
    simplex_search();
}

void carom::assess_node(int dex){
    if(dex<0 || dex>_chifn.get_pts()){
        printf("WARNING asking to assess node %d but only have %d pts\n",
        dex,_chifn.get_pts());
        
        exit(1);
    }
    
    if(_nodes.get_dim()==0 && _chifn.get_fn(dex)<_chifn.target()){
        _nodes.add(dex,&_chifn);
        return;
    }
    
    int keep_it,i,ix,so_far_so_good,iFound;
    double ftrial,dx;
    array_1d<double> trial;
    trial.set_name("carom_assess_node_trial");
    
    keep_it=1;
    for(ix=0;ix<_nodes.get_dim() && keep_it==1;ix++){
        so_far_so_good=0;
        for(dx=0.25;dx<0.8 && so_far_so_good==0;dx+=0.25){
            for(i=0;i<_chifn.get_dim();i++){
                trial.set(i,dx*_chifn.get_pt(dex,i)+(1.0-dx)*_chifn.get_pt(_nodes(i)->get_center(),i));
            }
            _chifn.evaluate(trial,&ftrial,&iFound);
            
            if(ftrial>_chifn.target()){
                so_far_so_good=1;
            }
        }

        if(so_far_so_good==0){
            keep_it=0;
        } 
    }
    
    if(keep_it==1){
        _nodes.add(dex,&_chifn);
    }
    
}


///////////////////////gp_cost

gp_cost::gp_cost(){
    _called=0;
    _ell=-1.0;
    _covarin.set_name("gp_cost_covarin");
    _chifn=NULL;
}

gp_cost::~gp_cost(){}

void gp_cost::set_chifn(chisq_wrapper *cc){
    _chifn=cc;
}

void gp_cost::is_it_safe(char *word){
    if(_chifn==NULL){
        printf("WARNING in gp_cost::%s\n",word);
        printf("_chifn is null\n");
        exit(1);
    }
}

int gp_cost::get_called(){
    return _called;
}

double gp_cost::operator()(array_1d<double> &pt){
    _called++;
    return 10.0;
}
