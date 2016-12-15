#include "maps.h"

maps::maps(){
    _ct_dalex=0;
    _last_did_min=0;
    sprintf(_outname,"output/carom_output.sav");
    sprintf(_timingname,"output/carom_timing.sav");
    _good_points.set_name("maps_good_points");
    _good_points.set_room(100000);
}

maps::~maps(){
    _chifn.write_pts();
}

double maps::evaluate(array_1d<double> &pt, int *dex){
    double mu;
    _chifn.evaluate(pt,&mu,dex);

    if(mu<_chifn.target() && dex[0]>=0){
        if(_good_points.contains(dex[0])==0){
            _good_points.add(dex[0]);
        }
    }

    return mu;
}

void maps::assess_good_points(){
    assess_good_points(0,_chifn.get_pts());
}

void maps::assess_good_points(int i_min){
    assess_good_points(i_min,_chifn.get_pts());
}

void maps::assess_good_points(int i_min, int i_max){
    int i;

    for(i=0;i<_good_points.get_dim();i++){
        if(_chifn.get_fn(_good_points.get_data(i))>_chifn.target()){
            _good_points.remove(i);
            i--;
        }
    }

    for(i=i_min;i<i_max+1 && i<_chifn.get_pts();i++){
        if(_chifn.get_fn(i)<_chifn.target()){
            if(_good_points.contains(i)==0){
                _good_points.add(i);
            }
        }
    }
}

void maps::set_dof(int dd){
    _chifn.set_dof(dd);
}

void maps::set_confidence_limit(double cc){
    _chifn.set_confidence_limit(cc);
}

void maps::set_seed(int ii){
    _chifn.set_seed(ii);
}

void maps::set_min(array_1d<double> &vv){
    _chifn.set_min(vv);
}

void maps::set_max(array_1d<double> &vv){
    _chifn.set_max(vv);
}

void maps::set_characteristic_length(int dex, double vv){
    _chifn.set_characteristic_length(dex,vv);
}

void maps::set_chisquared(chisquared *xx){
    _chifn.set_chisquared(xx);
}

void maps::set_deltachi(double xx){
    _chifn.set_deltachi(xx);
}

void maps::set_target(double tt){
    _chifn.set_target(tt);
}

void maps::set_write_every(int ww){
    _chifn.set_write_every(ww);
}

int maps::get_dim(){
    return _chifn.get_dim();
}

int maps::get_called(){
    return _chifn.get_called();
}

double maps::get_chimin(){
    return _chifn.chimin();
}

void maps::set_outname(char *nn){
    int i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _outname[i]=nn[i];
    }
    _outname[i]=0;
}

void maps::set_timingname(char *nn){
    int i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _timingname[i]=nn[i];
    }
    _timingname[i]=0;
}

void maps::initialize(int npts){
    _chifn.initialize(npts);
    _cloud.build(&_chifn);
    assess_good_points(0);
}

void maps::mcmc_init(){
    maps_initializer initializer;
    initializer.set_chifn(&_chifn);
    initializer.search();
}

void maps::search(int limit){
    int pt_start;

    double min0=_chifn.chimin();
    _chifn.set_outname(_outname);
    _chifn.set_timingname(_timingname);
    printf("before init min %e\n",_chifn.chimin());
    mcmc_init();
    printf("min now %e -> %e\n",min0,_chifn.chimin());
    printf("called %d\n",_chifn.get_pts());

    _cloud.set_limit(limit);

    while(_chifn.get_pts()<limit){

        pt_start=_chifn.get_pts();
        _cloud.search();
        _ct_dalex+=_chifn.get_pts()-pt_start;
    }

    int i;
    printf("minpt -- %e\n",_chifn.get_fn(_chifn.mindex()));
    for(i=0;i<_chifn.get_dim();i++){
        printf("    %.3e\n",_chifn.get_pt(_chifn.mindex(),i));
    }
    printf("\n\n");
    _chifn.write_pts();
}
