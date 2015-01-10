#include "carom.h"

carom::carom(){
    _write_every=3000;
    _last_written=0;
    sprintf(_outname,"output/carom_output.sav");
    sprintf(_timingname,"output/carom_timing.sav");
    _time_started=double(time(NULL));
}

carom::~carom(){}

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
        fprintf(output,"%.18e\n",chifn.get_fn(i));
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
