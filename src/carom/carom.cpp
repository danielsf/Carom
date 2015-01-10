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
    chifn.set_seed(ii);
}

void carom::set_min(array_1d<double> &vv){
    chifn.set_min(vv);
}

void carom::set_max(array_1d<double> &vv){
    chifn.set_max(vv);
}

void carom::set_characteristic_length(int dex, double vv){
    chifn.set_characteristic_length(dex,vv);
}

void carom::set_chisquared(chisquared *xx){
    chifn.set_chisquared(xx);
}

void carom::set_deltachi(double xx){
    chifn.set_deltachi(xx);
}

void carom::set_target(double tt){
    chifn.set_target(tt);
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
    chifn.initialize(npts);
    write_pts();
}

void carom::write_pts(){
    FILE *output;
    
    int i,j;
    output=fopen(_outname,"w");
    for(i=0;i<chifn.get_pts();i++){
        for(j=0;j<chifn.get_dim();j++){
            fprintf(output,"%.18e ",chifn.get_pt(i,j));
        }
        fprintf(output,"%.18e 0 0 1\n",chifn.get_fn(i));
    }
    fclose(output);
    
    output=fopen(_timingname,"a");
    fprintf(output,"%d %e %e\n",
        chifn.get_called(),
        double(time(NULL))-_time_started,
        (double(time(NULL))-_time_started)/double(chifn.get_called()));
    fclose(output);
    

    _last_written=chifn.get_called();
}

void carom::simplex_search(){

    simplex_minimizer ffmin;
    ffmin.set_chisquared(&chifn);
    ffmin.set_dice(chifn.get_dice());
    array_1d<double> min,max;
    min.set_name("carom_simplex_search_min");
    max.set_name("carom_simplex_search_min");
    chifn.get_min(min);
    chifn.get_max(max);
    ffmin.set_minmax(min,max);
    ffmin.use_gradient();
    
    array_1d<double> minpt;
    minpt.set_name("carom_simplex_search_minpt");
    
    array_2d<double> seed;
    seed.set_name("carom_simplex_search_seed");
    
    seed.set_cols(chifn.get_dim());
    int i,j;
    for(i=0;i<chifn.get_dim()+1;i++){
        for(j=0;j<chifn.get_dim();j++){
            seed.set(i,j,chifn.get_pt(i,j));
        }
    }
    
    ffmin.find_minimum(seed,minpt);
    write_pts();
}

void carom::search(){
    simplex_search();
}
