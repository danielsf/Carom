#include "dalex_driver.h"

dalex_driver::dalex_driver(){
    _log_file_name[0]=0;
    _ct_dalex=0;
    _last_did_min=0;
    sprintf(_outname,"output/carom_output.sav");
    sprintf(_timingname,"output/carom_timing.sav");
    _good_points.set_name("dale_driver_good_points");
    _good_points.set_room(100000);
}

dalex_driver::~dalex_driver(){
    _chifn.write_pts();
}

void dalex_driver::set_log_file_name(char *in){
    int i;
    for(i=0;i<letters-1 && in[i]!=0;i++){
        _log_file_name[i]=in[i];
    }
    _log_file_name[i]=0;
}

double dalex_driver::evaluate(array_1d<double> &pt, int *dex){
    double mu;
    _chifn.evaluate(pt,&mu,dex);

    if(mu<_chifn.target() && dex[0]>=0){
        if(_good_points.contains(dex[0])==0){
            _good_points.add(dex[0]);
        }
    }

    return mu;
}

void dalex_driver::assess_good_points(){
    assess_good_points(0,_chifn.get_pts());
}

void dalex_driver::assess_good_points(int i_min){
    assess_good_points(i_min,_chifn.get_pts());
}

void dalex_driver::assess_good_points(int i_min, int i_max){
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

void dalex_driver::set_dof(int dd){
    _chifn.set_dof(dd);
}

void dalex_driver::set_confidence_limit(double cc){
    _chifn.set_confidence_limit(cc);
}

void dalex_driver::set_seed(int ii){
    FILE *log_file;
    if(_log_file_name[0]!=0){
        log_file = fopen(_log_file_name, "a");
        fprintf(log_file,"setting seed to %d\n",ii);
        fclose(log_file);
    }
    _chifn.set_seed(ii);
}

void dalex_driver::set_min(array_1d<double> &vv){
    _chifn.set_min(vv);
}

void dalex_driver::set_max(array_1d<double> &vv){
    _chifn.set_max(vv);
}

void dalex_driver::set_characteristic_length(int dex, double vv){
    _chifn.set_characteristic_length(dex,vv);
}

void dalex_driver::set_chisquared(chisquared *xx){
    _chifn.set_chisquared(xx);
}

void dalex_driver::set_deltachi(double xx){
    _chifn.set_deltachi(xx);
}

void dalex_driver::set_target(double tt){
    _chifn.set_target(tt);
}

void dalex_driver::set_write_every(int ww){
    _chifn.set_write_every(ww);
}

int dalex_driver::get_dim(){
    return _chifn.get_dim();
}

int dalex_driver::get_called(){
    return _chifn.get_called();
}

double dalex_driver::get_chimin(){
    return _chifn.chimin();
}

void dalex_driver::set_outname(char *nn){
    int i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _outname[i]=nn[i];
    }
    _outname[i]=0;
}

void dalex_driver::set_timingname(char *nn){
    int i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _timingname[i]=nn[i];
    }
    _timingname[i]=0;
}

void dalex_driver::initialize(int npts){
    _chifn.initialize(npts);
    _cloud.build(&_chifn);
    assess_good_points(0);
}

void dalex_driver::mcmc_init(){
    dalex_initializer initializer;
    initializer.set_chifn(&_chifn);
    initializer.search();
}

void dalex_driver::search(int limit){
    double min0=_chifn.chimin();
    _chifn.set_outname(_outname);
    _chifn.set_timingname(_timingname);
    printf("before init min %e\n",_chifn.chimin());
    _chifn.set_search_type(_type_init);
    mcmc_init();
    printf("min now %e -> %e\n",min0,_chifn.chimin());
    printf("called %d\n",_chifn.get_pts());
    _search(limit);
}

void dalex_driver::warm_start(char *warm_name, int limit){
    _chifn.set_outname(_outname);
    _chifn.set_timingname(_timingname);

    double mu;
    int i,j;
    char word[letters];
    array_1d<double>pt;
    pt.set_name("dalex_driver_warm_pt");
    FILE *in_file;
    in_file=fopen(warm_name, "r");
    for(i=0;i<get_dim()+3;i++){
        fscanf(in_file,"%s",word);
    }
    while(fscanf(in_file,"%le",&mu)>0){
        pt.set(0,mu);
        for(i=1;i<get_dim();i++){
            fscanf(in_file,"%le",&mu);
            pt.set(i,mu);
        }
        fscanf(in_file,"%le",&mu);
        fscanf(in_file,"%d",&j);
        _chifn.set_search_type(j);
        _chifn.add_pt(pt,mu);
    }
    fclose(in_file);
    printf("read in %d pts\n",_chifn.get_pts());
    _chifn.rebalance();
    _search(limit);
}


void dalex_driver::_search(int limit){
    int pt_start;
    _cloud.set_limit(limit);
    _cloud.set_log_file_name(_log_file_name);

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
