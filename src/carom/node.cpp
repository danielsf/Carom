#include "node.h"

node::node(){
    initialize();
}

node::~node(){
}

node::node(const node &in){
    initialize();
    copy(in);
}

node& node::operator=(const node &in){
    if(this==&in) return *this;
    initialize();
    copy(in);
    return *this;
}

void node::deactivate(){
    _active=0;
}

void node::deactivate_simplex(){
    _do_simplex=0;
}

void node::initialize(){
    _log=NULL;
    _chisquared=NULL;
    _ricochet_growth=1.0;
    _baseline_volume=2.0*exception_value;
    _mcmc_growth=1.0;
    _swarm_growth=1.0;
    _simplex_growth=1.0;
    _compass_growth=1.0;
    _chimin=2.0*exception_value;
    _chimin_ricochet=2.0*exception_value;
    _do_simplex=1;
    _centerdex=-1;
    _geo_centerdex=-1;
    _ellipse_center=-1;
    _centerdex_basis=-1;
    _first_centerdex=-1;
    _min_changed=0;
    _active=1;
    _found_bases=0;
    _ct_ricochet=0;
    _ct_simplex=0;
    _allowed_ricochet_strikes=3;
    _v0=-1.0;
    _since_expansion=0;
    _min_basis_error=exception_value;
    _min_basis_error_changed=0;
    _failed_simplexes=0;
    _failed_kicks=0;
    _successful_kicks=0;
    _volume_of_last_geom=0.0;
    _last_geo_center=-1;
    _successful_ricochets=0;
    _id_dex=0;
    _good_shots=0;
    _bad_shots=0;
    _node_dd_tol=1.0e-2;
    _total_kicks=0;
    _total_trimmed=0;
    _convergence_ct=0;
    _compass_calls=0;
    _since_culled=0;

    _total_bisections=0;
    _bisection_calls=0;
    _total_ricochets=0;
    _ricochet_calls=0;
    _ricochet_bisection_calls=0;
    _ricochet_bisections=0;
    _gradient_calls=0;
    _highball_calls=0;
    _proper_ricochets=0;

    _swarm_acceptances=0;
    _swarm_rejections=0;
    _swarm_step=1.0;
    _swarm_outsiders=0;
    _swarm_expanders=0;

    _mcmc_step=0.1;
    _mcmc_acceptances=0;
    _mcmc_rejections=0;

    _compass_points.set_name("node_compass_points");
    _off_center_compass_points.set_name("node_off_center_compass_points");
    _off_center_origins.set_name("node_off_center_origins");
    _off_center_candidates.set_name("node_off_center_candidates");
    _basis_associates.set_name("node_basis_associates");
    _basis_mm.set_name("node_basis_mm");
    _basis_bb.set_name("node_basis_bb");
    _basis_model.set_name("node_basis_model");
    _basis_vectors.set_name("node_basis_vectors");
    _basis_ddsq.set_name("node_basis_ddsq");
    _basis_vv.set_name("node_basis_vv");
    _basis_lengths.set_name("node_basis_lengths");
    _max_found.set_name("node_max_found");
    _min_found.set_name("node_min_found");
    _projected_min.set_name("node_projected_min");
    _projected_max.set_name("node_projected_max");
    _ricochet_particles.set_name("node_ricochet_particles");
    _ricochet_origins.set_name("node_ricochet_origins");
    _associates.set_name("node_associates");
    _boundary_points.set_name("node_boundary_points");
    _true_min.set_name("node_true_min");
    _true_max.set_name("node_true_max");
    _transform.set_name("node_transform");
    _transform_associates.set_name("node_transform_associates");
    _avg_pts.set_name("node_avg_pts");
    _swarm.set_name("node_swarm");
    _swarm_center.set_name("node_swarm_center");
    _swarm_norm.set_name("node_swarm_norm");
    _swarm_associates.set_name("node_swarm_associates");
    _ricochet_strikes.set_name("node_ricochet_strikes");
}

void node::copy(const node &in){
    if(this==&in){
        return;
    }
    _convergence_ct=in._convergence_ct;
    _centerdex=in._centerdex;
    _first_centerdex=in._first_centerdex;
    _geo_centerdex=in._geo_centerdex;
    _centerdex_basis=in._centerdex_basis;
    _chimin=in._chimin;
    _chimin_ricochet=in._chimin_ricochet;
    _chimin_bases=in._chimin_bases;
    _min_changed=in._min_changed;
    _active=in._active;
    _failed_simplexes=in._failed_simplexes;
    _found_bases=in._found_bases;
    _ct_ricochet=in._ct_ricochet;
    _ct_simplex=in._ct_simplex;
    _allowed_ricochet_strikes=in._allowed_ricochet_strikes;
    _ellipse_center=in._ellipse_center;
    _since_expansion=in._since_expansion;
    _min_basis_error=in._min_basis_error;
    _min_basis_error_changed=in._min_basis_error_changed;
    _volume_of_last_geom=in._volume_of_last_geom;
    _last_geo_center=in._last_geo_center;
    _successful_ricochets=in._successful_ricochets;
    _id_dex=in._id_dex;
    _node_dd_tol=in._node_dd_tol;
    _good_shots=in._good_shots;
    _bad_shots=in._bad_shots;
    _total_kicks=in._total_kicks;
    _total_trimmed=in._total_trimmed;
    _highball_calls=in._highball_calls;
    _compass_calls=in._compass_calls;
    _do_simplex=in._do_simplex;
    _since_culled=in._since_culled;

    _total_bisections=in._total_bisections;
    _total_ricochets=in._total_ricochets;
    _bisection_calls=in._bisection_calls;
    _ricochet_calls=in._ricochet_calls;
    _ricochet_bisection_calls=in._ricochet_bisection_calls;
    _ricochet_bisections=in._ricochet_bisections;
    _gradient_calls=in._gradient_calls;
    _proper_ricochets=in._proper_ricochets;

    _ricochet_growth=in._ricochet_growth;
    _mcmc_growth=in._mcmc_growth;
    _swarm_growth=in._swarm_growth;
    _simplex_growth=in._simplex_growth;
    _compass_growth=in._compass_growth;
    _baseline_volume=in._baseline_volume;

    _swarm_acceptances=in._swarm_acceptances;
    _swarm_rejections=in._swarm_rejections;
    _swarm_step=in._swarm_step;
    _swarm_outsiders=in._swarm_outsiders;
    _swarm_expanders=in._swarm_expanders;

    _mcmc_step=in._mcmc_step;
    _mcmc_acceptances=in._mcmc_acceptances;
    _mcmc_rejections=in._mcmc_rejections;

    _v0=in._v0;

    int i,j;

    _chisquared=in._chisquared;
    _log=in._log;

    _associates.reset();
    for(i=0;i<in._associates.get_dim();i++){
        _associates.add(in._associates.get_data(i));
    }

    _boundary_points.reset();
    for(i=0;i<in._boundary_points.get_dim();i++){
        _boundary_points.add(in._boundary_points.get_data(i));
    }

    _compass_points.reset();
    for(i=0;i<in._compass_points.get_dim();i++){
        _compass_points.set(i,in._compass_points.get_data(i));
    }

    _swarm.reset();
    _swarm.set_cols(in._swarm.get_cols());
    for(i=0;i<in._swarm.get_rows();i++){
        for(j=0;j<in._swarm.get_cols();j++){
            _swarm.set(i,j,in._swarm.get_data(i,j));
        }
    }

    _ricochet_strikes.reset();
    for(i=0;i<in._ricochet_strikes.get_dim();i++){
        _ricochet_strikes.set(i,in._ricochet_strikes.get_data(i));
    }

    _avg_pts.reset();
    for(i=0;i<in._avg_pts.get_dim();i++){
        _avg_pts.set(i,in._avg_pts.get_data(i));
    }

    _off_center_compass_points.reset();
    for(i=0;i<in._off_center_compass_points.get_dim();i++){
        _off_center_compass_points.set(i,in._off_center_compass_points.get_data(i));
    }

    _off_center_origins.reset();
    for(i=0;i<in._off_center_origins.get_dim();i++){
        _off_center_origins.set(i,in._off_center_origins.get_data(i));
    }

    _off_center_candidates.reset();
    for(i=0;i<in._off_center_candidates.get_dim();i++){
        _off_center_candidates.set(i,in._off_center_candidates.get_data(i));
    }

    _basis_associates.reset();
    for(i=0;i<in._basis_associates.get_dim();i++){
        _basis_associates.set(i,in._basis_associates.get_data(i));
    }

    _basis_mm.reset();
    for(i=0;i<in._basis_mm.get_dim();i++){
        _basis_mm.set(i,in._basis_mm.get_data(i));
    }

    _basis_bb.reset();
    for(i=0;i<in._basis_bb.get_dim();i++){
        _basis_bb.set(i,in._basis_bb.get_data(i));
    }

    _basis_model.reset();
    for(i=0;i<in._basis_model.get_dim();i++){
        _basis_model.set(i,in._basis_model.get_data(i));
    }

    _basis_lengths.reset();
    for(i=0;i<in._basis_lengths.get_dim();i++){
        _basis_lengths.set(i,in._basis_lengths.get_data(i));
    }

    _basis_vectors.reset();
    _basis_vectors.set_cols(in._basis_vectors.get_cols());
    for(i=0;i<in._basis_vectors.get_rows();i++){
        for(j=0;j<in._basis_vectors.get_cols();j++){
            _basis_vectors.set(i,j,in._basis_vectors.get_data(i,j));
        }
    }

    _true_min.reset();
    for(i=0;i<in._true_min.get_dim();i++){
        _true_min.set(i,in._true_min.get_data(i));
    }

    _true_max.reset();
    for(i=0;i<in._true_max.get_dim();i++){
        _true_max.set(i,in._true_max.get_data(i));
    }

    _transform.reset();
    for(i=0;i<in._transform.get_dim();i++){
        _transform.set(i,in._transform.get_data(i));
    }

    _transform_associates.reset();
    for(i=0;i<in._transform_associates.get_rows();i++){
        for(j=0;j<in._transform_associates.get_cols(i);j++){
            _transform_associates.set(i,j,in._transform_associates.get_data(i,j));
        }
    }

    _min_found.reset();
    for(i=0;i<in._min_found.get_dim();i++){
        _min_found.set(i,in._min_found.get_data(i));
    }

    _max_found.reset();
    for(i=0;i<in._max_found.get_dim();i++){
        _max_found.set(i,in._max_found.get_data(i));
    }

    _projected_min.reset();
    for(i=0;i<in._projected_min.get_dim();i++){
        _projected_min.set(i,in._projected_min.get_data(i));
    }

    _projected_max.reset();
    for(i=0;i<in._projected_max.get_dim();i++){
        _projected_max.set(i,in._projected_max.get_data(i));
    }

    _ricochet_particles.reset();
    for(i=0;i<in._ricochet_particles.get_dim();i++){
        _ricochet_particles.set(i,in._ricochet_particles.get_data(i));
    }

}

void node::merge(const node &in){
    if(this==&in){
        return;
    }

    int i;
    for(i=0;i<in._associates.get_dim();i++){
        if(_associates.contains(in._associates.get_data(i))==0){
            _associates.add(in._associates.get_data(i));
        }
    }

    for(i=0;i<in._boundary_points.get_dim();i++){
        if(_boundary_points.contains(in._boundary_points.get_data(i))==0){
            _boundary_points.add(in._boundary_points.get_data(i));
        }
    }

    for(i=0;i<in._off_center_origins.get_dim();i++){
        _off_center_origins.add(in._off_center_origins.get_data(i));
    }

    for(i=0;i<in._off_center_candidates.get_dim();i++){
        _off_center_candidates.add(in._off_center_candidates.get_data(i));
    }

    int old_n_particles=_ricochet_particles.get_dim();
    int j;
    for(i=0;i<in._ricochet_particles.get_dim();i++){
        _ricochet_particles.add(in._ricochet_particles.get_data(i));
        _ricochet_origins.add(in._ricochet_origins.get_data(i));
    }

    for(i=0;i<in._ricochet_strikes.get_dim();i++){
        _ricochet_strikes.add(in._ricochet_strikes.get_data(i));
    }

    if(in._chimin<_chimin){
        _chimin=in._chimin;
    }

    if(_chisquared->get_fn(in._centerdex)<_chisquared->get_fn(_centerdex)){
        _centerdex=in._centerdex;
    }

    array_1d<double> buffer;
    for(i=0;i<in._swarm.get_rows();i++){
        for(j=0;j<in._swarm.get_cols();j++){
            buffer.set(j,in._swarm.get_data(i,j));
        }
        _swarm.add_row(buffer);
    }

    _ricochet_growth*=in._ricochet_growth;
    _mcmc_growth*=in._mcmc_growth;
    _swarm_growth*=in._swarm_growth;
    _simplex_growth*=in._simplex_growth;
    _compass_growth*=in._compass_growth;

    recalibrate_max_min();
    _baseline_volume=volume();
    _active=1;
}

int node::get_convergence_ct(){
    return _convergence_ct;
}

double node::get_ricochet_growth(){
    return _ricochet_growth;
}

void node::increment_ricochet_growth(){
    double vv=volume();
    if(vv>_baseline_volume){
        _ricochet_growth*=vv/_baseline_volume;
    }
    _baseline_volume=vv;
}

double node::get_swarm_growth(){
    return _swarm_growth;
}

void node::increment_swarm_growth(){
   double vv=volume();
   if(vv>_baseline_volume){
       _swarm_growth*=vv/_baseline_volume;
   }
   _baseline_volume=vv;
}

double node::get_mcmc_growth(){
    return _mcmc_growth;
}

void node::increment_mcmc_growth(){
    double vv=volume();
    if(vv>_baseline_volume){
        _mcmc_growth*=vv/_baseline_volume;
    }
    _baseline_volume=vv;
}

double node::get_compass_growth(){
    return _compass_growth;
}

void node::increment_compass_growth(){
    double vv=volume();
    if(vv>_baseline_volume){
        _compass_growth*=vv/_baseline_volume;
    }
    _baseline_volume=vv;
}

double node::get_simplex_growth(){
    return _simplex_growth;
}

void node::increment_simplex_growth(){
    double vv=volume();
    if(vv>_baseline_volume){
        _simplex_growth*=vv/_baseline_volume;
    }
    _baseline_volume=vv;
}

int node::get_good_shots(){
    return _good_shots;
}

int node::get_bad_shots(){
    return _bad_shots;
}

int node::get_successful_ricochets(){
    return _successful_ricochets;
}

int node::get_total_kicks(){
    return _total_kicks;
}

int node::get_total_trimmed(){
    return _total_trimmed;
}

int node::get_failed_kicks(){
    return _failed_kicks;
}

int node::get_successful_kicks(){
    return _successful_kicks;
}

int node::get_activity(){
    return _active;
}

int node::get_ct_ricochet(){
    return _ct_ricochet;
}

void node::set_id_dex(int ii){
    _id_dex=ii;
}

double node::get_transform(int dex){
    return _transform.get_data(dex);
}

double node::get_norm(int dex){
    if(dex>=_max_found.get_dim() || dex>=_min_found.get_dim()){
        return 1.0;
    }

    double ans;
    ans=_max_found.get_data(dex)-_min_found.get_data(dex);
    if(ans<1.0e-20){
        return 1.0;
    }

    return ans;
}

double node::get_projected_norm(int dex){
    if(dex>=_projected_max.get_dim() || dex>=_projected_min.get_dim()){
        return 1.0;
    }

    double ans;
    ans=_projected_max.get_data(dex)-_projected_min.get_data(dex);
    if(ans<1.0e-20){
        return 1.0;
    }

    return ans;
}

void node::set_center(int ix){
    _centerdex=ix;
    _first_centerdex=ix;
    _min_changed=1;
    if(_chisquared!=NULL){
        _chimin=_chisquared->get_fn(ix);
    }
}

void node::set_log(asymm_array_2d<int> *ll){
    _log=ll;
}

void node::add_to_log(int kind, int pt){
    if(_log==NULL){
        return;
    }

    if(kind<_log->get_rows() && _log[0](kind)->contains(pt)==1){
        return;
    }

    _log->add(kind,pt);
}

void node::set_chisquared(chisq_wrapper *cc){
    _chisquared=cc;
    if(_centerdex>=0){
        _chimin=_chisquared->get_fn(_centerdex);
    }

    int i,j;
    _basis_vectors.reset();
    _basis_vectors.set_cols(_chisquared->get_dim());
    for(i=0;i<_chisquared->get_dim();i++){
        for(j=0;j<_chisquared->get_dim();j++){
            if(i==j){
                _basis_vectors.set(i,j,1.0);
            }
            else{
                _basis_vectors.set(i,j,0.0);
            }
        }
    }

    for(i=0;i<_chisquared->get_dim();i++){
        _transform.set(i,1.0);
    }

}

int node::get_n_boundary(){
    return _boundary_points.get_dim();
}

int node::get_boundary(int ii){
    return _boundary_points.get_data(ii);
}

int node::get_n_associates(){
    return _associates.get_dim();
}

int node::get_associate(int ii){
    return _associates.get_data(ii);
}

int node::get_center(){
    return _centerdex;
}

double node::volume(){
    is_it_safe("volume");
    if(_true_min.get_dim()!=_chisquared->get_dim() || _true_max.get_dim()!=_chisquared->get_dim()){
        return 0.0;
    }

    double ans;
    int i;
    ans=1.0;
    for(i=0;i<_chisquared->get_dim();i++){
        ans*=(_true_max.get_data(i)-_true_min.get_data(i));
    }

    return ans;
}

double node::projected_volume(){
    is_it_safe("projected_volume");
    if(_projected_min.get_dim()!=_chisquared->get_dim() || _projected_max.get_dim()!=_chisquared->get_dim()){
        return 0.0;
    }

    double ans;
    int i;
    ans=1.0;
    for(i=0;i<_chisquared->get_dim();i++){
        ans*=(_projected_max.get_data(i)-_projected_min.get_data(i));
    }

    return ans;
}

void node::set_basis(int ii, int jj, double vv){
    if(_chisquared==NULL){
        printf("WARNING cannot set node basis before assigning _chisquared\n");
        exit(1);
    }
    _basis_vectors.set(ii,jj,vv);
}

void node::is_it_safe(char *word){
    if(_chisquared==NULL){
        printf("WARNING in node::%s\n",word);
        printf("_chisquared is null\n");
        exit(1);
    }

    if(_basis_vectors.get_cols()!=_chisquared->get_dim() ||
       _basis_vectors.get_rows()!=_chisquared->get_dim()){

        printf("WARNING in node::%s\n",word);

        printf("_basis_vectors %d by %d\n",_basis_vectors.get_rows(),
        _basis_vectors.get_cols());

        printf("but dim should be %d\n",_chisquared->get_dim());

        exit(1);

    }

    if(_ricochet_particles.get_dim()!=_ricochet_origins.get_dim()){

        printf("WARNING in node::%s\n",word);
        printf("ricochet particles %d\n",_ricochet_particles.get_dim());
        printf("ricochet origins %d\n",_ricochet_origins.get_dim());

        exit(1);

    }

    if(_ricochet_strikes.get_dim()!=_ricochet_particles.get_dim()){
        printf("WARNING %d strikes but %d particles\n",
        _ricochet_strikes.get_dim(),
        _ricochet_particles.get_dim());

        exit(1);
    }

}

int node::is_this_an_associate_gross(int dex){
    array_1d<double> trial;
    trial.set_name("node_is_this_associate_gross_trial");
    double mu;
    int iFound;
    int i;
    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,0.5*(get_pt(dex,i)+get_pt(_first_centerdex,i)));
    }
    evaluate(trial,&mu,&iFound);
    if(mu<=_chisquared->target()){
         return 1;
    }
    return 0;
}

int node::is_this_an_associate(int dex){

    double dd,ddmin;
    int i,ibest;
    ibest=-1;
    for(i=0;i<_associates.get_dim();i++){
        if(dex==_associates.get_data(i))return 1;
        dd=node_distance(dex,_associates.get_data(i));
        if(ibest<0 || dd<ddmin){
            ibest=_associates.get_data(i);
            ddmin=dd;
        }
    }

    dd=node_distance(dex,_centerdex);
    if(dd<ddmin || ibest<0){
        ibest=_centerdex;
    }

    if(ibest<0){
        return 0;
    }

    array_1d<double> trial;
    trial.set_name("node_is_this_trial");
    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,0.5*(get_pt(ibest,i)+get_pt(dex,i)));
    }
    double mu;
    int newDex,j;
    evaluate(trial,&mu,&newDex);

    if(mu<_chisquared->target()){
        _associates.add(dex);
        return 1;
    }

    return 0;

}

double node::get_pt(int dex, int ii){
    return _chisquared->get_pt(dex, ii)/_transform.get_data(ii);
}

void node::transform_pt_to_node(array_1d<double> &pt_in, array_1d<double> &pt_out){
    int i;
    for(i=0;i<_chisquared->get_dim();i++){
        pt_out.set(i,pt_in.get_data(i)/_transform.get_data(i));
    }
}

void node::transform_pt_to_truth(array_1d<double> &pt_in, array_1d<double> &pt_out){
    int i;
    for(i=0;i<_chisquared->get_dim();i++){
        pt_out.set(i,pt_in.get_data(i)*_transform.get_data(i));
    }
}

void node::set_transform(){
    int i,ipt;
    int ix;
    array_1d<double> dd,dd_sorted;
    array_1d<int> dd_dexes;

    dd.set_name("node_set_transform_dd");
    dd_sorted.set_name("node_set_transform_dd_sorted");
    dd_dexes.set_name("node_set_transform_dd_dexes");

    double mu;

    for(ix=0;ix<_chisquared->get_dim();ix++){
        dd.reset_preserving_room();
        dd_sorted.reset_preserving_room();
        dd_dexes.reset_preserving_room();
        for(i=0;i<_transform_associates.get_cols(ix);i++){
            ipt=_transform_associates.get_data(ix,i);
            mu=0.5*fabs(_chisquared->get_pt(_centerdex,ix)-_chisquared->get_pt(ipt,ix));
            if(mu>1.0e-10*(_true_max.get_data(ix)-_true_min.get_data(ix))){
                dd.add(mu);
                dd_dexes.add(i);
            }
        }

        if(dd.get_dim()==1){
            _transform.set(ix,dd.get_data(0));
        }
        else if(dd.get_dim()>1){
            sort_and_check(dd,dd_sorted,dd_dexes);
            _transform.set(ix,dd_sorted.get_data(dd.get_dim()/2));
        }
    }
    find_bases(0);
    set_geo_center();
    recalibrate_max_min();

}

void node::get_true_pt(int i1, array_1d<double> &pt_out){
    int i;
    for(i=0;i<_chisquared->get_dim();i++){
        pt_out.set(i,_chisquared->get_pt(i1,i));
    }
}

void node::evaluate(array_1d<double> &pt_in, double *value, int *dex){
    is_it_safe("evaluate");

    array_1d<double> pt_true;
    pt_true.set_name("evaluate_pt_true");

    transform_pt_to_truth(pt_in,pt_true);
    _chisquared->evaluate(pt_true,value,dex);

    int i,j;
    array_1d<double> projected;
    projected.set_name("node_evaluate_projected");

    if(dex[0]>=0){
        if(value[0]<_chimin){
            _chimin=value[0];
            _centerdex=dex[0];
            _min_changed=1;
            if(_first_centerdex<0 || _id_dex==0){
                _first_centerdex=dex[0];
            }
        }

        if(value[0]<=_chisquared->target()){
            j=1;
            for(i=0;i<_associates.get_dim() && j==1;i++){
                if(_associates.get_data(i)==dex[0])j=0;
            }

            if(j==1){
                _associates.add(dex[0]);
            }

            for(i=0;i<pt_true.get_dim();i++){
                if(i>=_true_min.get_dim() || pt_true.get_data(i)<_true_min.get_data(i)){
                    _true_min.set(i,pt_true.get_data(i));
                }

                if(i>=_true_max.get_dim() || pt_true.get_data(i)>_true_max.get_data(i)){
                    _true_max.set(i,pt_true.get_data(i));
                }
            }

            for(i=0;i<pt_in.get_dim();i++){
                if(i>=_min_found.get_dim() || pt_in.get_data(i)<_min_found.get_data(i)){
                    _min_found.set(i,pt_in.get_data(i));
                }

                if(i>=_max_found.get_dim() || pt_in.get_data(i)>_max_found.get_data(i)){
                    _max_found.set(i,pt_in.get_data(i));
                }
            }

            project_to_bases(pt_in,projected);
            for(i=0;i<projected.get_dim();i++){
                if(i>=_projected_min.get_dim() || projected.get_data(i)<_projected_min.get_data(i)){
                    _projected_min.set(i,projected.get_data(i));
                }

                if(i>=_projected_max.get_dim() || projected.get_data(i)>_projected_max.get_data(i)){
                    _projected_max.set(i,projected.get_data(i));
                }
            }

            if(value[0]<0.5*(_chisquared->target()+_chisquared->get_fn(_centerdex))){
                if(value[0]>0.75*_chisquared->get_fn(_centerdex)+0.25*_chisquared->target()){
                    _off_center_candidates.add(dex[0]);
                }
            }
        }

    }
}

double node::normalized_node_distance(int i1, int i2){
    double ans=0.0;
    int i;
    for(i=0;i<_chisquared->get_dim();i++){
        ans+=power((get_pt(i1,i)-get_pt(i2,i))/get_norm(i),2);
    }
    return sqrt(ans);
}

double node::node_distance(array_1d<double> &p1, array_1d<double> &p2){
    int i;
    double ans;
    ans=0.0;
    for(i=0;i<_chisquared->get_dim();i++){
        ans+=power(p1.get_data(i)-p2.get_data(i),2);
    }
    return sqrt(ans);

}

double node::node_distance(int i1, int i2){
    array_1d<double> p1, p2;
    int i;
    for(i=0;i<_chisquared->get_dim();i++){
        p1.set(i,get_pt(i1,i));
        p2.set(i,get_pt(i2,i));
    }
    return node_distance(p1,p2);
}

double node::node_distance(int i1, array_1d<double> &p2){
    array_1d<double> p1;
    int i;
    for(i=0;i<_chisquared->get_dim();i++){
        p1.set(i,get_pt(i1,i));
    }
    return node_distance(p2, p1);
}

double node::node_second_derivative_different(int center, int ix, int iy){
    is_it_safe("node_second_derivative_different");

    double xnorm,ynorm;
    xnorm=1.0;
    ynorm=1.0;

    int ifpp,ifpm,ifmp,ifmm;
    ifpp=-1;
    ifpm=-1;
    ifmp=-1;
    ifmm=-1;

    double mu;
    array_1d<double> trial;
    trial.set_name("node_second_derivative_different_trial");

    int i;
    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,get_pt(center,i));
    }

    double xp,xm,yp,ym,fpp,fmp,fpm,fmm;
    double dx,xcenter,ycenter;

    int proceed,xpBound,xmBound,ypBound,ymBound;

    dx=1.0e-2;
    while(ifpp==ifpm || ifpp==ifmp || ifpp==ifmm ||
          ifpm==ifmp || ifpm==ifmm || ifmp==ifmm){


         proceed=0;

         while(proceed==0){
             xp=get_pt(center,ix)+dx*xnorm;
             xm=get_pt(center,ix)-dx*xnorm;

             yp=get_pt(center,iy)+dx*ynorm;
             ym=get_pt(center,iy)-dx*ynorm;

             xpBound=_chisquared->in_bounds(ix,xp);
             xmBound=_chisquared->in_bounds(ix,xm);
             ypBound=_chisquared->in_bounds(iy,yp);
             ymBound=_chisquared->in_bounds(iy,ym);

             if(xpBound==0 && xmBound==0){
                 throw -1;
             }

             if(ypBound==0 && ymBound==0){
                 throw -1;
             }

             proceed=1;
             if(xpBound==0){
                 xcenter=get_pt(center,ix)-1.5*dx*xnorm;
                 proceed=0;
             }
             else if(xmBound==0){
                 xcenter=get_pt(center,ix)+1.5*dx*xnorm;
                 proceed=0;
             }
             else{
                 xcenter=get_pt(center,ix);
             }

             if(ypBound==0){
                 ycenter=get_pt(center,iy)-1.5*dx*ynorm;
                 proceed=0;
             }
             else if(ymBound==0){
                 ycenter=get_pt(center,iy)+1.5*dx*ynorm;
                 proceed=0;
             }
             else{
                 ycenter=get_pt(center,iy);
             }

             if(proceed==0){
                 trial.set(ix,xcenter);
                 trial.set(iy,ycenter);
                 evaluate(trial,&mu,&center);
                 if(center<0){
                     throw -1;
                 }
             }

         }

         trial.set(ix,xp);
         trial.set(iy,yp);
         evaluate(trial,&fpp,&ifpp);

         trial.set(iy,ym);
         evaluate(trial,&fpm,&ifpm);

         trial.set(ix,xm);
         evaluate(trial,&fmm,&ifmm);

         trial.set(iy,yp);
         evaluate(trial,&fmp,&ifmp);

         if(ifpp<0 && ifpm<0 && ifmp<0 && ifmm<0){
             ifpp=-1;
             ifpm=-2;
             ifmp=-3;
             ifpm=-4;
         }

         dx*=1.5;

    }

    double ans;
    ans=(fpp+fmm-fmp-fpm)/((xp-xm)*(yp-ym));

    return ans;
}

double node::node_second_derivative_same(int center, int ix){
    is_it_safe("node_second_derivative_same");

    double xnorm;

    xnorm=1.0;

    double dx;
    dx=1.0e-2;

    int ifpp,ifmm;
    double fpp,fmm,xp,xpp,xm,xmm;

    ifpp=-1;
    ifmm=-1;

    int i;
    array_1d<double> trial;
    trial.set_name("node_second_derivative_same_trial");
    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,get_pt(center,i));
    }

    int proceed;
    int xppBound,xmmBound;
    double xcenter,ycenter,mu;

    while(ifpp==center || ifpp==ifmm || ifmm==center){
        proceed=0;
        while(proceed==0){
            xp=get_pt(center,ix)+dx*xnorm;
            xm=get_pt(center,ix)-dx*xnorm;
            xpp=get_pt(center,ix)+2.0*dx*xnorm;
            xmm=get_pt(center,ix)-2.0*dx*xnorm;

            xppBound = _chisquared->in_bounds(ix,xpp);
            xmmBound = _chisquared->in_bounds(ix,xmm);

            if(xppBound==0 && xmmBound==0){
                throw -1;
            }

            proceed=1;
            if(xppBound==0){
                xcenter=get_pt(center,ix)-2.5*dx*xnorm;
                proceed=0;
            }
            else if(xmmBound==0){
                xcenter=get_pt(center,ix)+2.5*dx*xnorm;
                proceed=0;
            }

            if(proceed==0){
                trial.set(ix,xcenter);
                evaluate(trial,&mu,&center);
                if(center<0){
                    throw -1;
                }
            }

        }

        trial.set(ix,xpp);
        evaluate(trial,&fpp,&ifpp);

        trial.set(ix,xmm);
        evaluate(trial,&fmm,&ifmm);

        if(ifpp<0 && ifmm<0){
            ifpp=-1;
            ifmm=-2;
        }

        dx*=1.5;

    }

    double ans;
    ans=(fpp+fmm-2.0*_chisquared->get_fn(center))/power(xp-xm,2);
    return ans;
}

double node::node_second_derivative(int center, int ix, int iy){
    if(center<0){
        return 0.0;
    }

    if(ix==iy){
        return node_second_derivative_same(center,ix);
    }
    else{
        return node_second_derivative_different(center,ix,iy);
    }
}

void node::node_gradient(int dex, array_1d<double> &grad){
    int ibefore=_chisquared->get_called();
    int tried=0;

    tried=_node_1sided_projected_gradient(dex,grad);

    if(tried==0){
        _node_2sided_gradient(dex,grad);
    }

    _gradient_calls+=_chisquared->get_called()-ibefore;
}

void node::_node_1sided_gradient(int dex, array_1d<double> &grad){
    is_it_safe("node_gradient");

    int i,j,if1,if2;
    array_1d<double> trial;
    double norm,x1,x2,y1,y2;
    trial.set_name("node_node_gradient_trial");

    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,get_pt(dex,i));
    }

    double dx,dxstart;
    dxstart=1.0e-2;

    for(i=0;i<_chisquared->get_dim();i++){
        norm=1.0;

        dx=dxstart;
        if1=-1;
        if2=-1;
        while((if1<0 || if1==dex) && (if2<0 || if2==dex)){
            x1=get_pt(dex,i)+dx*norm;
            trial.set(i,x1);
            evaluate(trial,&y1,&if1);

            if(if1<0 && if1!=dex){
                x2=get_pt(dex,i)-dx*norm;
                trial.set(i,x2);
                evaluate(trial,&y2,&if2);
            }

            if((if1>=0 && if1!=dex) || (if2>=0 && if2!=dex)){
                if(if1>=0 && if1!=dex){
                    grad.set(i,(y1-_chisquared->get_fn(dex))/(x1-get_pt(dex,i)));
                }
                else{
                    grad.set(i,(y2-_chisquared->get_fn(dex))/(x2-get_pt(dex,i)));
                }
            }
            else{
                printf("%d %d %d\n",if1,if2,dex);
                dx*=1.5;
            }
        }

        trial.set(i,get_pt(dex,i));
    }

}

void node::_node_2sided_gradient(int dex, array_1d<double> &grad){
    is_it_safe("node_gradient");

    int i,j,if1,if2;
    array_1d<double> trial;
    double norm,x1,x2,y1,y2;
    trial.set_name("node_node_gradient_trial");

    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,get_pt(dex,i));
    }

    double dx,dxstart;
    dxstart=1.0e-2;

    for(i=0;i<_chisquared->get_dim();i++){
        norm=1.0;

        dx=dxstart;
        if1=-1;
        if2=-1;
        while(if1==if2){
            x1=get_pt(dex,i)+dx*norm;
            trial.set(i,x1);
            evaluate(trial,&y1,&if1);

            x2=get_pt(dex,i)-dx*norm;
            trial.set(i,x2);
            evaluate(trial,&y2,&if2);

            if((if1<0 && if2>=0) or (if2<0 && if1>=0)){
                if(if1<0){
                    x1=get_pt(dex,i);
                    y1=_chisquared->get_fn(dex);
                    if1=dex;
                }
                else{
                    x2=get_pt(dex,i);
                    y2=_chisquared->get_fn(dex);
                    if2=dex;
                }
            }

            if(if1!=if2 || (if1<0 && if2<0)){
                grad.set(i,(y1-y2)/(x1-x2));
                if(if1<0 && if2<0){
                    if1=-1;
                    if2=-2;
                }
            }
            else{
                dx*=1.5;
            }
        }

        trial.set(i,get_pt(dex,i));
    }

}


int node::_node_1sided_projected_gradient(int dex, array_1d<double> &grad){

    if(_basis_vectors.get_rows()!=_chisquared->get_dim()){
        return 0;
    }

    array_1d<double> dir,projected_dir,sign;
    dir.set_name("node_proj_grad_dir");
    projected_dir.set_name("node_proj_grad_proj_dir");
    sign.set_name("node_proj_grad_sign");

    int local_center,i;
    local_center=find_local_center();
    for(i=0;i<_chisquared->get_dim();i++){
       dir.set(i,get_pt(dex,i)-get_pt(local_center,i));
    }

    project_to_bases(dir,projected_dir);
    for(i=0;i<_chisquared->get_dim();i++){
        if(projected_dir.get_data(i)<0.0){
            sign.set(i,1.0);
        }
        else{
            sign.set(i,-1.0);
        }
    }

    double dx,dxstart;
    dxstart=0.001;

    array_1d<double> trial;
    trial.set_name("node_proj_grad_trial");

    double mu,norm,y2;
    int ix,iFound;

    for(i=0;i<_chisquared->get_dim();i++){
        grad.set(i,0.0);
    }

    for(ix=0;ix<_chisquared->get_dim();ix++){
        dx=dxstart;

        norm=1.0;
        if(ix<_projected_max.get_dim()){
             if(_projected_max.get_data(ix)-_projected_min.get_data(ix)>0.0){
                 norm=(_projected_max.get_data(ix)-_projected_min.get_data(ix));
             }
        }

        iFound=-1;
        while(iFound<0 || iFound==dex){
            for(i=0;i<_chisquared->get_dim();i++){
                mu=get_pt(dex,i)+dx*sign.get_data(ix)*norm*_basis_vectors.get_data(ix,i);
                trial.set(i,mu);
            }

            evaluate(trial,&y2,&iFound);

            if(iFound<0 || iFound==dex){
                dx*=2.0;
            }

            if(dx>10.0*dxstart){
                return 0;
            }
        }

        mu=(y2-_chisquared->get_fn(dex))/(norm*dx*sign.get_data(ix));
        for(i=0;i<_chisquared->get_dim();i++){
            grad.add_val(i,mu*_basis_vectors.get_data(ix,i));
        }
    }

    return 1;
}

int node::node_bisection_origin_dir(int ii, array_1d<double> &dd){
    int iFound;
    iFound=node_bisection_origin_dir(ii,dd,_chisquared->target(),0.01*(_chisquared->target()-_chimin));
    if(iFound>=0){
        add_to_boundary(iFound);
    }
    return iFound;
}

int node::node_bisection_origin_dir(int iOrigin, array_1d<double> &dir, double target, double tol){
    array_1d<double> lowball,highball;
    double flow, fhigh;

    lowball.set_name("node_bisection_origin_lowball");
    highball.set_name("node_bisection_origin_highball");

    if(_chisquared->get_fn(iOrigin)>target){
        return -1;
    }

    flow=_chisquared->get_fn(iOrigin);
    int i;
    for(i=0;i<_chisquared->get_dim();i++){
        lowball.set(i,get_pt(iOrigin,i));
        highball.set(i,get_pt(iOrigin,i));
    }

    double norm=1.0;

    fhigh=-2.0*exception_value;
    while(fhigh<target){
        for(i=0;i<_chisquared->get_dim();i++){
            highball.add_val(i,norm*dir.get_data(i));
        }
        evaluate(highball,&fhigh,&i);
        norm*=2.0;
    }

    return node_bisection(lowball,flow,highball,fhigh,1,target,tol);
}

int node::node_bisection(array_1d<double> &ll, double fl,
                    array_1d<double> &hh, double fh,
                    int doSlope){

        int iFound;
        iFound=node_bisection(ll,fl,hh,fh,doSlope,_chisquared->target(),0.01*(_chisquared->target()-_chimin));
        if(iFound>=0){
            add_to_boundary(iFound);
        }
        return iFound;
}

int node::node_bisection(int i_low, int i_high, double target, double tol){
    array_1d<double> lowball,highball;
    lowball.set_name("node_bisection_int_lowball");
    highball.set_name("node_bisection_int_highball");
    double flow,fhigh;

    int i;
    flow=_chisquared->get_fn(i_low);
    fhigh=_chisquared->get_fn(i_high);

    if(fhigh<target || flow>target){
        return -1;
    }

    for(i=0;i<_chisquared->get_dim();i++){
        lowball.set(i,get_pt(i_low,i));
        highball.set(i,get_pt(i_high,i));
    }

    return node_bisection(lowball, flow, highball, fhigh, 1, target, tol);

}

int node::node_bisection(array_1d<double> &lowball_in, double flow,
                    array_1d<double> &highball_in, double fhigh,
                    int doSlope, double target_value, double tolerance){

    is_it_safe("bisection");

    int i;
    array_1d<double> lowball,highball;
    lowball.set_name("node_bisection_lowball");
    highball.set_name("node_bisection_highball");

    for(i=0;i<_chisquared->get_dim();i++){
        lowball.set(i,lowball_in.get_data(i));
        highball.set(i,highball_in.get_data(i));
    }

    if(flow>target_value){
        printf("cannot do bisection, flow %e target %e\n",flow,target_value);
        flow=_chisquared->get_fn(_centerdex);
        for(i=0;i<_chisquared->get_dim();i++){
            lowball.set(i,get_pt(_centerdex,i));
        }
    }

    int ibefore=_chisquared->get_called();

    array_1d<double> e_dir;
    double component;
    e_dir.set_name("node_bisection_e_dir");

    if(flow>fhigh || fhigh<target_value){
        printf("WARNING in node bisection flow %e fhigh %e\n",flow,fhigh);
        for(i=0;i<_chisquared->get_dim();i++){
            e_dir.set(i,highball.get_data(i)-lowball.get_data(i));
        }
        e_dir.normalize();
        component=1.0;
        while(fhigh<flow || fhigh<target_value){
             for(i=0;i<_chisquared->get_dim();i++){
                 highball.add_val(i,component*e_dir.get_data(i));
             }
             component*=2.0;
             evaluate(highball,&fhigh,&i);
        }
    }

    double ftrial;
    array_1d<double> trial;
    trial.set_name("node_bisection_trial");

    int took_a_step=0,ct,iout;
    double wgt;

    ct=0;
    iout=-1;
    while(ct<100 && (took_a_step==0 || target_value-flow>tolerance)){

        if(doSlope==1){
            wgt=(fhigh-target_value)/(fhigh-flow);
            if(wgt<0.1)wgt=0.1;
            else if(wgt>0.9)wgt=0.9;
        }
        else{
            wgt=0.5;
        }

        for(i=0;i<_chisquared->get_dim();i++){
            trial.set(i,wgt*lowball.get_data(i)+(1.0-wgt)*highball.get_data(i));
        }

        evaluate(trial,&ftrial,&i);
        if(ftrial<target_value){
            flow=ftrial;
            if(i>=0)iout=i;
            for(i=0;i<_chisquared->get_dim();i++){
                lowball.set(i,trial.get_data(i));
            }
            took_a_step=1;
            doSlope=1;
        }
        else{
            fhigh=ftrial;
            for(i=0;i<_chisquared->get_dim();i++){
                highball.set(i,trial.get_data(i));
            }
        }
        ct++;
    }

    _total_bisections++;
    _bisection_calls+=_chisquared->get_called()-ibefore;
    return iout;

}

void node::perturb_bases(int idim, array_1d<double> &dx, array_2d<double> &bases_out){

    is_it_safe("perturb_bases");

    int i,j;
    bases_out.set_cols(_basis_vectors.get_cols());
    for(i=0;i<_chisquared->get_dim();i++){
        for(j=0;j<_chisquared->get_dim();j++){
            bases_out.set(i,j,_basis_vectors.get_data(i,j));
        }
    }

    for(i=0;i<_chisquared->get_dim();i++){
        bases_out.add_val(idim,i,dx.get_data(i));
    }
    bases_out(idim)->normalize();

    int ix,jx;
    double mu;

    for(ix=idim+1;ix!=idim;){
        if(ix>=_chisquared->get_dim()){
            ix=0;
        }

        for(jx=idim;jx!=ix;){
            if(ix>=_chisquared->get_dim()){
                jx=0;
            }

            mu=0.0;
            for(i=0;i<_chisquared->get_dim();i++){
                mu+=bases_out.get_data(ix,i)*bases_out.get_data(jx,i);
            }
            for(i=0;i<_chisquared->get_dim();i++){
                bases_out.subtract_val(ix,i,mu*bases_out.get_data(jx,i));
            }

            if(jx<_chisquared->get_dim()-1)jx++;
            else jx=0;
        }

        bases_out(ix)->normalize();

        if(ix<_chisquared->get_dim()-1)ix++;
        else ix=0;

    }

    validate_bases(bases_out,"node_perturb_bases");

}


void node::validate_bases(array_2d<double> &bases, char *whereami){
    int ix,i,jx;
    double mu;
    /////////////////testing
    for(ix=0;ix<_chisquared->get_dim();ix++){
        bases(ix)->normalize();
        mu=0.0;
        for(i=0;i<_chisquared->get_dim();i++){
            mu+=bases.get_data(ix,i)*bases.get_data(ix,i);
        }
        if(fabs(mu-1.0)>1.0e-6){
            printf("WARNING in %s, square norm %e\n",whereami,mu);
            exit(1);
        }

        for(jx=ix+1;jx<_chisquared->get_dim();jx++){
            mu=0.0;
            for(i=0;i<_chisquared->get_dim();i++){
                mu+=bases.get_data(ix,i)*bases.get_data(jx,i);
            }

            if(fabs(mu)>1.0e-6){
                printf("WARNING in %s, dot product %e\n",whereami,mu);
                exit(1);
            }
        }
    }

}

double node::basis_error(array_2d<double> &trial_bases, array_1d<double> &trial_model){

    if(_basis_associates.get_dim()<=0){
        printf("WARNING cannot calculate basis error there are only %d associates\n",
        _basis_associates.get_dim());

        exit(1);
    }

    is_it_safe("basis_error");

    trial_model.zero();
    if(_basis_ddsq.get_rows()>_basis_associates.get_dim()){
        _basis_ddsq.reset();
    }

    if(_basis_ddsq.get_cols()!=_chisquared->get_dim()){
        _basis_ddsq.set_cols(_chisquared->get_dim());
    }

    if(_basis_mm.get_dim()!=_chisquared->get_dim()*_chisquared->get_dim()){
        _basis_mm.set_dim(_chisquared->get_dim()*_chisquared->get_dim());
    }

    if(_basis_bb.get_dim()!=_chisquared->get_dim()){
        _basis_bb.set_dim(_chisquared->get_dim());
    }

    if(_basis_vv.get_dim()!=_chisquared->get_dim()){
        _basis_vv.set_dim(_chisquared->get_dim());
    }

    _basis_mm.zero();
    _basis_bb.zero();
    _basis_vv.zero();
    _basis_ddsq.zero();

    int i,j,ix;
    double mu;
    for(ix=0;ix<_basis_associates.get_dim();ix++){
        for(i=0;i<_chisquared->get_dim();i++){
            mu=0.0;
            for(j=0;j<_chisquared->get_dim();j++){
                mu+=(get_pt(_basis_associates.get_data(ix),j)-get_pt(_centerdex,j))*trial_bases.get_data(i,j);
            }
            _basis_ddsq.set(ix,i,mu*mu);
        }
    }

    for(i=0;i<_chisquared->get_dim();i++){
        for(j=0;j<_basis_associates.get_dim();j++){
            _basis_bb.add_val(i,_basis_ddsq.get_data(j,i)*(_chisquared->get_fn(_basis_associates.get_data(j))-_chimin));
        }
    }

    int k;
    for(i=0;i<_chisquared->get_dim();i++){
        for(j=i;j<_chisquared->get_dim();j++){
            ix=i*_chisquared->get_dim()+j;
            for(k=0;k<_basis_associates.get_dim();k++){
                _basis_mm.add_val(ix,_basis_ddsq.get_data(k,i)*_basis_ddsq.get_data(k,j));
            }
            if(j!=i){
                _basis_mm.set(j*_chisquared->get_dim()+i,_basis_mm.get_data(ix));
            }
        }
    }

    try{
        naive_gaussian_solver(_basis_mm,_basis_bb,trial_model,_chisquared->get_dim());
    }
    catch(int iex){
        printf("WARNING basis_error was no good\n");
        return 2.0*exception_value;
    }

    double error=0.0,chi_model;
    for(i=0;i<_basis_associates.get_dim();i++){
        chi_model=_chimin;
        for(j=0;j<_chisquared->get_dim();j++){
            chi_model+=trial_model.get_data(j)*_basis_ddsq.get_data(i,j);
        }
        error+=power(_chisquared->get_fn(_basis_associates.get_data(i))-chi_model,2);
    }

    return error/double(_basis_associates.get_dim());

}

void node::add_to_compass(int dex){
    _compass_points.add(dex);
    add_to_log(_log_compass, dex);
}

void node::populate_basis_associates(){

    int n_associates_0=_basis_associates.get_dim();

    array_1d<double> dir;
    dir.set_name("populate_basis_associates_dir");
    int i,dimsq,iFound,iStart;

    double target,tol;

    dimsq=_chisquared->get_dim()*_chisquared->get_dim();
    while(_basis_associates.get_dim()<n_associates_0+2*dimsq){
        for(i=0;i<_chisquared->get_dim();i++){
            dir.set(i,normal_deviate(_chisquared->get_dice(),0.0,1.0));
        }
        dir.normalize();
        target=0.5*(_chisquared->target()+_chimin);
        tol=0.01*(_chisquared->target()-_chimin);

        iFound=node_bisection_origin_dir(_centerdex,dir,target,tol);
        if(iFound>=0 && fabs(_chisquared->get_fn(iFound)-target)<5.0*tol){
            if(_basis_associates.contains(iFound)==0){
                _basis_associates.add(iFound);
            }

            target=0.25*_chimin+0.75*_chisquared->target();

            if(_chisquared->get_fn(iFound)<target){
                iStart=iFound;
            }
            else{
                iStart=_centerdex;
            }

            iFound=node_bisection_origin_dir(iStart,dir,target,tol);
            if(iFound>=0 && fabs(_chisquared->get_fn(iFound)-target)<5.0*tol){
                if(_basis_associates.contains(iFound)==0){
                    _basis_associates.add(iFound);
                }
            }


        }

    }
}


void node::compass_search(){
    compass_search(_centerdex);
}

void node::compass_search(int local_center){

    is_it_safe("compass_search");
    _compass_points.reset();
    _transform_associates.reset();

    int ibefore=_chisquared->get_called();

    int ix,i,j,iFound;
    double sgn,flow,fhigh,dx,ftrial,step,bisection_target;
    double blength;
    array_1d<double> lowball,highball,trial;

    lowball.set_name("node_compass_search_lowball");
    highball.set_name("node_compass_search_highball");

    dx=0.0;
    for(ix=0;ix<_chisquared->get_dim();ix++){
        blength=0.0;
        for(sgn=-1.0;sgn<1.5;sgn+=2.0){
            flow=2.0*exception_value;
            fhigh=-2.0*exception_value;

            if(sgn<0.0){
                flow=_chisquared->get_fn(local_center);
                for(i=0;i<_chisquared->get_dim();i++){
                    lowball.set(i,get_pt(local_center,i));
                }
            }
            else{
                for(i=0;i<_chisquared->get_dim();i++){
                    trial.set(i,get_pt(local_center,i)+sgn*dx*_basis_vectors.get_data(ix,i));
                }
                evaluate(trial,&ftrial,&iFound);

                if(ftrial<_chisquared->target()){
                    flow=ftrial;
                    for(i=0;i<_chisquared->get_dim();i++){
                        lowball.set(i,trial.get_data(i));
                    }
                }
                else{
                    fhigh=ftrial;
                    for(i=0;i<_chisquared->get_dim();i++){
                        highball.set(i,trial.get_data(i));
                    }
                }
            }

            if(flow>_chisquared->target()){
                flow=_chisquared->get_fn(local_center);
                for(i=0;i<_chisquared->get_dim();i++){
                    lowball.set(i,get_pt(local_center,i));
                }
            }

            step=1.0;
            while(fhigh<=_chisquared->target()){
                for(i=0;i<_chisquared->get_dim();i++){
                    highball.set(i,lowball.get_data(i)+sgn*step*_basis_vectors.get_data(ix,i));
                }
                evaluate(highball,&fhigh,&iFound);
                step*=2.0;
            }

            iFound=-1;
            if(flow>_chisquared->target() || flow>fhigh){
                flow=_chisquared->get_fn(_centerdex);
                for(i=0;i<_chisquared->get_dim();i++){
                    lowball.set(i,get_pt(_centerdex,i));
                }
            }
            iFound=node_bisection(lowball,flow,highball,fhigh,1);

            dx=0.0;
            for(i=0;i<_chisquared->get_dim();i++){
                dx+=_basis_vectors.get_data(ix,i)*(get_pt(local_center,i)-get_pt(iFound,i));
            }

            dx=fabs(dx);

            if(iFound>=0){
                add_to_compass(iFound);
                if(fabs(_chisquared->get_fn(iFound)-_chisquared->target())<0.1*(_chisquared->target()-_chisquared->chimin())){
                    _transform_associates.add(ix,iFound);
                }
            }

        }
        _basis_lengths.set(ix,blength);
    }

    printf("leaving compass %d\n\n",_chisquared->get_called()-ibefore);
    _compass_calls+=_chisquared->get_called()-ibefore;
}

void node::compass_diagonal(int local_center){
    is_it_safe("compass_diagonal");

    int ix,iy;
    array_1d<double> trial,lowball,highball,dir;
    trial.set_name("node_off_diag_trial");
    lowball.set_name("node_off_diag_lowball");
    highball.set_name("node_off_diag_highball");
    dir.set_name("node_off_diag_dir");

    double dx,dy,dmin,step,bisection_target;
    double flow,fhigh,ftrial,mu;
    int i,j,k,iFound;
    int nGuessHigh=0,nGuessLow=0,isHigh;
    int spentHigh=0,spentLow=0,spentNeither=0,startFromMin=0;
    int ibefore;

    double sqrt2o2,xweight,yweight;
    double dmin_np,dmin_nn;

    sqrt2o2=0.5*sqrt(2.0);

    for(ix=0;ix<_chisquared->get_dim();ix++){
        i=2*ix;
        j=i+1;

        dx=-1.0;
        if(_compass_points.get_dim()>i && _compass_points.get_dim()>j){
            dx=0.0;
            for(k=0;k<_chisquared->get_dim();k++){
                dx+=0.5*(get_pt(j,k)-get_pt(i,k))*_basis_vectors.get_data(ix,k);
            }
        }

        if(dx<0.0){
            dx=2.0*exception_value;
        }

        for(iy=ix+1;iy<_chisquared->get_dim();iy++){
            i=2*iy;
            j=i+1;

            dy=-1.0;
            if(_compass_points.get_dim()>i && _compass_points.get_dim()>j){
                dy=0.0;
                for(k=0;k<_chisquared->get_dim();k++){
                    dy+=0.5*(get_pt(j,k)-get_pt(i,k))*_basis_vectors.get_data(iy,k);
                }
            }

            if(dy<0.0){
                dy=2.0*exception_value;
            }

            if(dx<dy){
                dmin=dx;
            }
            else{
                dmin=dy;
            }

            dmin_np=dmin;
            dmin_nn=dmin;

            for(xweight=-1.0*sqrt2o2;xweight<sqrt2o2*1.1;xweight+=2.0*sqrt2o2){
                for(yweight=-1.0*sqrt2o2;yweight<sqrt2o2*1.1;yweight+=2.0*sqrt2o2){
                    ibefore=_chisquared->get_called();
                    isHigh=-1;
                    if(xweight>0.0){
                        if(yweight<0.0)dmin=dmin_np;
                        if(yweight>0.0)dmin=dmin_nn;
                    }
                    for(i=0;i<_chisquared->get_dim();i++){
                        dir.set(i,xweight*_basis_vectors.get_data(ix,i)+yweight*_basis_vectors.get_data(iy,i));
                    }

                    flow=2.0*exception_value;
                    fhigh=-2.0*exception_value;
                    if(dmin<exception_value){
                        for(i=0;i<_chisquared->get_dim();i++){
                            trial.set(i,get_pt(local_center,i)+dmin*dir.get_data(i));
                        }
                        evaluate(trial,&ftrial,&iFound);

                        if(ftrial<_chisquared->target()){
                            nGuessLow++;
                            isHigh=0;
                            flow=ftrial;
                            for(i=0;i<_chisquared->get_dim();i++){
                                lowball.set(i,trial.get_data(i));
                            }
                        }
                        else{
                            nGuessHigh++;
                            isHigh=1;
                            fhigh=ftrial;
                            for(i=0;i<_chisquared->get_dim();i++){
                                highball.set(i,trial.get_data(i));
                            }
                        }
                    }

                    if(flow>_chisquared->target()){
                        startFromMin++;
                        flow=_chimin;
                        for(i=0;i<_chisquared->get_dim();i++){
                            lowball.set(i,get_pt(local_center,i));
                        }
                    }

                    step=1.0;
                    while(fhigh<_chisquared->target()){
                        for(i=0;i<_chisquared->get_dim();i++){
                            highball.set(i,lowball.get_data(i)+step*dir.get_data(i));
                        }
                        evaluate(highball,&fhigh,&iFound);
                        step*=2.0;
                    }

                    if(flow>_chisquared->target() || flow>fhigh){
                        flow=_chisquared->get_fn(local_center);
                        for(i=0;i<_chisquared->get_dim();i++){
                            lowball.set(i,get_pt(local_center,i));
                        }

                    }
                    iFound=node_bisection(lowball,flow,highball,fhigh,1);

                    if(iFound>=0){
                        dmin=0.0;
                        mu=0.0;
                        for(i=0;i<_chisquared->get_dim();i++){
                            mu+=(get_pt(local_center,i)-get_pt(iFound,i))*_basis_vectors.get_data(ix,i);
                        }
                        dmin+=mu*mu;
                        mu=0.0;
                        for(i=0;i<_chisquared->get_dim();i++){
                            mu+=(get_pt(local_center,i)-get_pt(iFound,i))*_basis_vectors.get_data(iy,i);
                        }
                        dmin+=mu*mu;
                        dmin=sqrt(dmin);
                        if(xweight<0.0 && yweight<0.0)dmin_nn=dmin;
                        if(xweight<0.0 && yweight>0.0)dmin_np=dmin;

                        add_to_compass(iFound);

                        if(fabs(_chisquared->get_fn(iFound)-_chisquared->target())<0.1*(_chisquared->target()-_chisquared->chimin())){
                            _transform_associates.add(ix,iFound);
                            _transform_associates.add(iy,iFound);
                        }

                    }

                    if(isHigh==1){
                        spentHigh+=_chisquared->get_called()-ibefore;
                    }
                    else if(isHigh==0){
                        spentLow+=_chisquared->get_called()-ibefore;
                    }
                    else{
                        spentNeither+=_chisquared->get_called()-ibefore;
                    }

                }
            }
        }
    }

    printf("nGuessHigh %d nGuessLow %d\n",nGuessHigh,nGuessLow);
    printf("spentHigh %d spentLow %d spentNeither %d\n",spentHigh,spentLow,spentNeither);
    printf("startFromMin %d\n",startFromMin);

}


void node::set_geo_center(){
    is_it_safe("set_geo_center");

    array_1d<double> geometric_dir,trial,highball,lowball;
    int iFound,i,j;
    double mu,fhigh,flow,bisection_target,tol;

    geometric_dir.set_name("set_geo_center_dir");
    trial.set_name("set_geo_center_trial");

    iFound=-1;

    trial.set_dim(_chisquared->get_dim());
    trial.zero();

    for(i=0;i<_chisquared->get_dim();i++){
        mu=0.5*(_projected_max.get_data(i)+_projected_min.get_data(i));
        for(j=0;j<_chisquared->get_dim();j++){
            trial.add_val(j,mu*_basis_vectors.get_data(i,j));
        }
    }

    evaluate(trial,&mu,&iFound);

    if(iFound>=0){
        printf("\nraw geometric center chisq %e\n",_chisquared->get_fn(iFound));
        for(i=0;i<_chisquared->get_dim();i++){
            printf("    %e\n",get_pt(iFound,i));
        }
        if(mu>_chisquared->target()){
            for(i=0;i<_chisquared->get_dim();i++){
                lowball.set(i,get_pt(_centerdex,i));
                highball.set(i,get_pt(iFound,i));
            }
            flow=_chisquared->get_fn(_centerdex);
            fhigh=_chisquared->get_fn(iFound);
            bisection_target=0.5*(_chisquared->get_fn(_centerdex)+_chisquared->target());
            tol=0.01*(_chisquared->target()-_chisquared->get_fn(_centerdex));

            iFound=node_bisection(lowball,flow,highball,fhigh,1,bisection_target,tol);

        }

        if(iFound>=0 && iFound!=_centerdex){
            for(i=0;i<_chisquared->get_dim();i++){
                printf("%e %e\n",get_pt(_centerdex,i),get_pt(iFound,i));
            }
            _geo_centerdex=iFound;
        }

    }
}


void node::compass_search_geometric_center(){
    set_geo_center();
    if(_geo_centerdex>=0 && _geo_centerdex!=_centerdex && _geo_centerdex!=_last_geo_center){
        compass_search(_geo_centerdex);
    }
    _volume_of_last_geom=volume();
    _last_geo_center=_geo_centerdex;
}


void node::findCovarianceMatrix(int iCenter, array_2d<double> &covar){
    is_it_safe("findCovarianceMatrix");


    array_1d<double> trial,center,norm;
    trial.set_name("node_findCovar_trial");
    norm.set_name("node_findCovar_norm");
    center.set_name("node_findCovar_center");

    double fcenter;
    int i;

    fcenter=_chisquared->get_fn(iCenter);

    for(i=0;i<_chisquared->get_dim();i++){
        center.set(i,get_pt(iCenter,i));
        norm.set(i,1.0);
    }

    array_2d<double> fpp,fpm,fmp,fmm;
    fpp.set_name("node_findCovar_fpp");
    fpm.set_name("node_findCovar_fpm");
    fmp.set_name("node_findCovar_fmp");
    fmm.set_name("node_findCovar_fmm");

    fpp.set_cols(_chisquared->get_dim());
    fpm.set_cols(_chisquared->get_dim());
    fmp.set_cols(_chisquared->get_dim());
    fmm.set_cols(_chisquared->get_dim());

    array_1d<double> f2p,f2m;
    f2p.set_name("node_findCovar_f2p");
    f2m.set_name("node_findCovar_f2m");

    int ifpp,ifpm,ifmp,ifmm,if2p,if2m;

    double mu;
    array_1d<double> dx;
    dx.set_name("node_findCovar_dx");
    for(i=0;i<_chisquared->get_dim();i++){
        dx.set(i,1.0e-2);
    }

    int ix,iy,keepGoing,ctAbort,ctAbortMax,calledMax;
    int ibefore=_chisquared->get_called();

    ctAbort=0;
    ctAbortMax=100;
    calledMax = 10*_chisquared->get_dim()*_chisquared->get_dim();

    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,center.get_data(i));
    }

    for(ix=0;ix<_chisquared->get_dim();ix++){
        for(i=0;i<_chisquared->get_dim();i++){
            trial.set(i,center.get_data(i));
        }

        if(_chisquared->get_called()-ibefore>calledMax ||
           ctAbort>=ctAbortMax){
                printf("Could not find CoVar; aborting\n");
                printf("ctAbort %d\n",ctAbort);
                printf("called %d\n",_chisquared->get_called()-ibefore);
                throw -1;
        }

        if(iCenter<0){
            printf("Center is invalid; aborting\n");
            throw -1;
        }

        keepGoing=1;
        trial.set(ix,center.get_data(ix)+2.0*dx.get_data(ix)*norm.get_data(ix));
        evaluate(trial,&mu,&if2p);

        if(if2p>=0 && if2p!=iCenter){
            f2p.set(ix,mu);
        }
        else if(if2p==iCenter){
            dx.multiply_val(ix,1.5);
            keepGoing=0;
            ix--;
            ctAbort++;
        }
        else if(if2p<0){
            center.subtract_val(ix,2.5*dx.get_data(ix)*norm.get_data(ix));
            evaluate(center,&fcenter,&iCenter);
            keepGoing=0;
            ix--;
            ctAbort++;
        }

        if(keepGoing==1){
            trial.set(ix,center.get_data(ix)-2.0*dx.get_data(ix)*norm.get_data(ix));
            evaluate(trial,&mu,&if2m);

            if(if2m>=0 && if2m!=iCenter){
                f2m.set(ix,mu);
            }
            else if(if2m==iCenter){
                dx.multiply_val(ix,1.5);
                keepGoing=0;
                ix--;
                ctAbort++;
            }
            else if(if2m<0){
                center.add_val(ix,2.5*dx.get_data(ix)*norm.get_data(ix));
                evaluate(center,&fcenter,&iCenter);
                keepGoing=0;
                ix--;
                ctAbort++;
            }
        }

        for(iy=ix-1;iy>=0 && keepGoing==1;iy--){
            for(i=0;i<_chisquared->get_dim();i++){
                trial.set(i,center.get_data(i));
            }

            if(_chisquared->get_called()-ibefore>calledMax ||
               ctAbort>=ctAbortMax){
                printf("Could not find CoVar; aborting\n");
                printf("ctAbort %d\n",ctAbort);
                printf("called %d\n",_chisquared->get_called()-ibefore);
                throw -1;
            }

            if(iCenter<0){
                printf("center is invalid; aborting\n");
                throw -1;
            }

            trial.set(ix,center.get_data(ix)+dx.get_data(ix)*norm.get_data(ix));
            trial.set(iy,center.get_data(iy)+dx.get_data(iy)*norm.get_data(iy));
            evaluate(trial,&mu,&ifpp);
            if(ifpp>=0 && ifpp!=iCenter){
                fpp.set(ix,iy,mu);
            }
            else if(ifpp==iCenter){
                dx.multiply_val(ix,1.5);
                dx.multiply_val(iy,1.5);
                keepGoing=0;
                ix--;
                ctAbort++;
            }
            else if(ifpp<0){
                center.subtract_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                center.subtract_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                evaluate(center,&fcenter,&iCenter);
                keepGoing=0;
                ix--;
                ctAbort++;
            }

            if(keepGoing==1){
               trial.set(iy,center.get_data(iy)-dx.get_data(iy)*norm.get_data(iy));
               evaluate(trial,&mu,&ifpm);
               if(ifpm>=0 && ifpm!=iCenter){
                   fpm.set(ix,iy,mu);
               }
               else if(ifpm==iCenter){
                   dx.multiply_val(ix,1.5);
                   dx.multiply_val(iy,1.5);
                   keepGoing=0;
                   ix--;
                   ctAbort++;
               }
               else if(ifpm<0){
                   center.subtract_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                   center.add_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                   evaluate(center,&fcenter,&iCenter);
                   keepGoing=0;
                   ix--;
                   ctAbort++;
               }
            }

            if(keepGoing==1){
                trial.set(ix,center.get_data(ix)-dx.get_data(ix)*norm.get_data(ix));
                evaluate(trial,&mu,&ifmm);
                if(ifmm>=0 && ifmm!=iCenter){
                    fmm.set(ix,iy,mu);
                }
                else if(ifmm==iCenter){
                    dx.multiply_val(ix,1.5);
                    dx.multiply_val(iy,1.5);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                }
                else if(ifmm<0){
                    center.add_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                    center.add_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                    evaluate(center,&fcenter,&iCenter);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                }
            }

            if(keepGoing==1){
                trial.set(iy,center.get_data(iy)+dx.get_data(iy)*norm.get_data(iy));
                evaluate(trial,&mu,&ifmp);
                if(ifmp>=0 && ifmp!=iCenter){
                    fmp.set(ix,iy,mu);
                }
                else if(ifmp==iCenter){
                    dx.multiply_val(ix,1.5);
                    dx.multiply_val(iy,1.5);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                }
                else if(ifmp<0){
                    center.add_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                    center.subtract_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                    evaluate(center,&fcenter,&iCenter);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                }
            }

        }
    }

    covar.set_cols(_chisquared->get_dim());
    for(ix=0;ix<_chisquared->get_dim();ix++){
        covar.set(ix,ix,0.25*(f2p.get_data(ix)+f2m.get_data(ix)-2.0*fcenter)/power(dx.get_data(ix)*norm.get_data(ix),2));
    }

    double num,denom;
    for(ix=0;ix<_chisquared->get_dim();ix++){
        for(iy=ix-1;iy>=0;iy--){
            num=0.25*(fpp.get_data(ix,iy)+fmm.get_data(ix,iy)-fmp.get_data(ix,iy)-fpm.get_data(ix,iy));
            denom=dx.get_data(ix)*norm.get_data(ix)*dx.get_data(iy)*norm.get_data(iy);
            covar.set(ix,iy,num/denom);
            covar.set(iy,ix,num/denom);
        }
    }


}

int node::findAcceptableCenter(){
    is_it_safe("findAcceptableCenter");

    array_1d<double> step;
    step.set_name("node_findAcceptableCenter_step");
    double norm,dx;

    dx=1.0e-2;
    int i;
    for(i=0;i<_chisquared->get_dim();i++){
        norm=1.0;
        step.set(i,dx*norm);
    }

    int ans=_centerdex;

    int ix,acceptable,ct,ctmax,locallyAcceptable,iFound;
    double mu,xtrial;
    array_1d<double> center;
    array_1d<double> movement;
    center.set_name("node_findAcceptableCenter_center");
    movement.set_name("node_findAcceptableCenter_movement");

    for(i=0;i<_chisquared->get_dim();i++){
        center.set(i,get_pt(ans,i));
    }

    ct=0;
    ctmax=100;
    acceptable=0;
    while(acceptable==0 && ct<ctmax){
        ct++;
        for(ix=0;ix<_chisquared->get_dim();ix++){
            movement.set(ix,0.0);
            xtrial=center.get_data(ix)+2.0*step.get_data(ix);
            locallyAcceptable=_chisquared->in_bounds(ix, xtrial);
            if(locallyAcceptable==0){
                movement.set(ix,-2.5);
            }
            else{
                xtrial=center.get_data(ix)-2.0*step.get_data(ix);
                locallyAcceptable=_chisquared->in_bounds(ix,xtrial);
                if(locallyAcceptable==0){
                    movement.set(ix,2.5);
                }
            }
        }

        norm=movement.get_square_norm();
        if(norm<1.0){
            acceptable=1;
        }
        else{

            for(i=0;i<_chisquared->get_dim();i++){
                center.add_val(i,movement.get_data(i)*step.get_data(i));
            }
            evaluate(center,&mu,&iFound);

            if(iFound>=0){
                ans=iFound;
            }
            else{
                printf("need to abort; new trial center was not in bounds\n");
                ct=ctmax+1;
            }
        }

    }


    if(ct>=ctmax){
        printf("WARNING did not find acceptable center; had to abort because of ct\n");
    }
    return ans;
}

void node::guess_bases(array_2d<double> &bases){
    is_it_safe("guess_bases");

    array_2d<double> covar;
    covar.set_name("node_guess_bases_covar");
    covar.set_cols(_chisquared->get_dim());
    bases.set_cols(_chisquared->get_dim());
    int ix,iy;
    double covarmax=-1.0;

    findCovarianceMatrix(_centerdex,covar);

    for(ix=0;ix<_chisquared->get_dim();ix++){
        for(iy=ix;iy<_chisquared->get_dim();iy++){
            if(fabs(covar.get_data(ix,iy))>covarmax){
                covarmax=fabs(covar.get_data(ix,iy));
            }
        }
    }

    for(ix=0;ix<_chisquared->get_dim();ix++){
        for(iy=0;iy<_chisquared->get_dim();iy++){
            covar.divide_val(ix,iy,covarmax);
        }
    }

    printf("assembled covariance matrix\n");

    array_2d<double> evecs;
    evecs.set_name("node_guess_bases_evecs");
    array_1d<double> evals;
    evals.set_name("node_guess_bases_evals");

    evecs.set_cols(_chisquared->get_dim());

    try{
        eval_symm(covar,evecs,evals,0.1);
    }
    catch(int iex){
        printf("Guess failed on of eigen vectors\n");
        throw -1;
    }

    for(ix=0;ix<_chisquared->get_dim();ix++){
        for(iy=0;iy<_chisquared->get_dim();iy++){
            bases.set(ix,iy,evecs.get_data(ix,iy));
        }
        bases(ix)->normalize();
    }

    validate_bases(bases,"node_guess_bases");

    printf("validated guessed bases\n");
}

void node::project_to_bases(array_1d<double> &in, array_1d<double> &out){
    is_it_safe("project_to_bases");
    int i;

    if(_basis_vectors.get_rows()!=_chisquared->get_dim()){
        for(i=0;i<_chisquared->get_dim();i++){
            out.set(i,in.get_data(i));
        }
        return;
    }

    int j;
    for(i=0;i<_chisquared->get_dim();i++){
        out.set(i,0.0);
        for(j=0;j<_chisquared->get_dim();j++){
            out.add_val(i,in.get_data(j)*_basis_vectors.get_data(i,j));
        }
    }
}

void node::recalibrate_max_min(){
    is_it_safe("recalibrate_max_min");

    _projected_max.reset();
    _projected_min.reset();
    _min_found.reset();
    _max_found.reset();
    _true_min.reset();
    _true_max.reset();

    array_1d<double> projected;
    projected.set_name("node_recalibrate_projected");
    int i,j;
    double mu;
    double tol=0.01;

    array_1d<double> trial;
    trial.set_name("recalibrate_trial");
    int ix;

    for(i=0;i<_associates.get_dim();i++){
        if(_chisquared->get_fn(_associates.get_data(i))<=_chisquared->target()+tol){

            for(ix=0;ix<_chisquared->get_dim();ix++){
                trial.set(ix, get_pt(_associates.get_data(i),ix));
            }
            project_to_bases(trial,projected);

            for(j=0;j<_chisquared->get_dim();j++){
                if(j>=_projected_min.get_dim() || projected.get_data(j)<_projected_min.get_data(j)){
                    _projected_min.set(j,projected.get_data(j));
                }

                if(j>=_projected_max.get_dim() || projected.get_data(j)>_projected_max.get_data(j)){
                    _projected_max.set(j,projected.get_data(j));
                }

                mu=get_pt(_associates.get_data(i),j);
                if(j>=_min_found.get_dim() || mu<_min_found.get_data(j)){
                    _min_found.set(j,mu);
                }

                if(j>=_max_found.get_dim() || mu>_max_found.get_data(j)){
                    _max_found.set(j,mu);
                }

                mu=_chisquared->get_pt(_associates.get_data(i),j);
                if(j>=_true_min.get_dim() || mu<_true_min.get_data(j)){
                    _true_min.set(j,mu);
                }

                if(j>=_true_max.get_dim() || mu>_true_max.get_data(j)){
                    _true_max.set(j,mu);
                }
            }
        }
        else{
            _associates.remove(i);
            i--;
        }
    }

    for(i=0;i<_boundary_points.get_dim();i++){
        if(_chisquared->get_fn(_boundary_points.get_data(i))>_chisquared->target()+tol){
            _boundary_points.remove(i);
            i--;
        }
    }

    _v0=volume();
}

void node::find_bases(){
    find_bases(1);
}

void node::find_bases(int do_populate){
    is_it_safe("find_bases");

    if(_centerdex_basis<0 || node_distance(_centerdex_basis,_centerdex)>0.01){
        _basis_associates.reset();
    }

    _chimin_bases=_chimin;
    _ellipse_center=_centerdex;

    int i,j;
    for(i=0;i<_basis_associates.get_dim();i++){
        if(_chisquared->get_fn(_basis_associates.get_data(i))-_chimin < 1.0e-6){
            _basis_associates.remove(i);
            i--;
        }
    }

    int iFound;
    double flow,fhigh,target,mu,tol;
    array_1d<double> lowball,highball,dir;
    lowball.set_name("node_find_bases_lowball");
    highball.set_name("node_find_bases_highball");
    dir.set_name("node_find_bases_dir");

    int dimsq=_chisquared->get_dim()*_chisquared->get_dim();

    if(do_populate==1 || _basis_associates.get_dim()<_chisquared->get_dim()){
        populate_basis_associates();
    }

    if(_basis_associates.get_dim()==0){
        printf("WARNING _basis associates is empty\n");
        exit(1);
    }

    array_2d<double> trial_bases;
    array_1d<double> trial_model,dx;

    trial_bases.set_name("node_find_bases_trial_bases");
    trial_model.set_name("node_find_bases_trial_model");
    dx.set_name("node_find_bases_dx");

    int ct,idim,aborted,max_abort,total_aborted,changed_bases;
    double error0,error,errorBest,stdev,stdevlim,error1;

    stdev=0.1/sqrt(double(_chisquared->get_dim()));
    stdevlim=1.0e-5/sqrt(double(_chisquared->get_dim()));
    max_abort=_chisquared->get_dim()*100;

    error=basis_error(_basis_vectors,_basis_model);
    error0=error;
    error1=error;
    errorBest=error;
    aborted=0;
    total_aborted=0;
    changed_bases=0;
    ct=0;

    printf("center %e min %e\n",_chisquared->get_fn(_centerdex),_chisquared->chimin());
    printf("error0 %e %d -- vol %e proj %e\n",
    error0,_basis_associates.get_dim(),volume(),projected_volume());

    if(_centerdex!=_centerdex_basis){
        printf("guessing basis from second derivative\n");
        try{
            guess_bases(trial_bases);
            error=basis_error(trial_bases,trial_model);
            printf("guess got error %e\n",error);
            if(error<error0){
                changed_bases=1;
                for(i=0;i<_chisquared->get_dim();i++){
                    _basis_model.set(i,trial_model.get_data(i));
                    for(j=0;j<_chisquared->get_dim();j++){
                        _basis_vectors.set(i,j,trial_bases.get_data(i,j));
                    }
                }
                error1=error;
                errorBest=error;
            }
        }
        catch(int iex){
            printf("never mind; guess did not work\n");
        }
    }


    while(stdev>stdevlim && aborted<max_abort && ct<10000){
        ct++;
        idim=-1;
        while(idim>=_chisquared->get_dim() || idim<0){
            idim=_chisquared->random_int()%_chisquared->get_dim();
        }

        for(i=0;i<_chisquared->get_dim();i++){
            dx.set(i,normal_deviate(_chisquared->get_dice(),0.0,stdev));
        }

        perturb_bases(idim,dx,trial_bases);
        error=basis_error(trial_bases,trial_model);

        if(error<errorBest){
            if(error1-error>1.0e-5*error){
                aborted=0;
            }
            else{
                aborted++;
                total_aborted++;
            }

            changed_bases=1;
            for(i=0;i<_chisquared->get_dim();i++){
                _basis_model.set(i,trial_model.get_data(i));
                for(j=0;j<_chisquared->get_dim();j++){
                    _basis_vectors.set(i,j,trial_bases.get_data(i,j));
                }
            }
            errorBest=error;
        }
        else{
            aborted++;
            total_aborted++;
        }

        if(ct%(max_abort/2)==0){
            if(total_aborted<(3*ct)/4)stdev*=1.5;
            else if(total_aborted>(3*ct)/4)stdev*=0.5;
        }

        if(ct%1000==0){
            error1=errorBest;
            printf("    ct %d error %e from %e min %e\n",ct,errorBest,error0,_chimin);
        }
    }

    array_1d<double> projected;
    projected.set_name("node_find_bases_projected");

    recalibrate_max_min();

    _found_bases++;
    _min_changed=0;
    _centerdex_basis=_centerdex;
    printf("done finding bases -- vol %e proj %e changed %d\n",
    volume(),projected_volume(),changed_bases);
    if(errorBest<_min_basis_error){
        _min_basis_error=errorBest;
        _min_basis_error_changed=1;
    }
    _v0=volume();
}


void node::off_center_compass(int iStart){

    if(iStart<0){
        return;
    }

    printf("\noff center compass begins\n");

    int ix;
    int i,j,k;
    double dd,ddmin;

    array_2d<double> dir;
    dir.set_name("off_center_compass_directions");
    array_1d<double> trial_dir;
    trial_dir.set_name("off_center_compass_trial_dir");
    dir.set_cols(_chisquared->get_dim());

    double component,norm;

    while(dir.get_rows()<_chisquared->get_dim()){
        for(i=0;i<_chisquared->get_dim();i++){
            trial_dir.set(i,normal_deviate(_chisquared->get_dice(),0.0,1.0));
        }
        for(i=0;i<dir.get_rows();i++){
            component=0.0;
            for(j=0;j<_chisquared->get_dim();j++){
                component+=trial_dir.get_data(j)*dir.get_data(i,j);
            }
            for(j=0;j<_chisquared->get_dim();j++){
                trial_dir.subtract_val(j,component*dir.get_data(i,j));
            }
        }
        norm=trial_dir.normalize();
        if(norm>1.0e-20){
            dir.add_row(trial_dir);
        }
    }

    _off_center_origins.add(iStart);
    int ibefore=_chisquared->get_called();

    array_1d<double> lowball,highball,trial;
    double flow,fhigh,ftrial,dx;
    int iFound;

    lowball.set_name("node_off_center_compass_lowball");
    highball.set_name("node_off_center_compass_highball");
    trial.set_name("node_off_center_compass_trial");

    double sgn;
    for(ix=0;ix<_chisquared->get_dim();ix++){
        dx=1.0;
        for(sgn=-1.0;sgn<1.1;sgn+=2.0){
            flow=2.0*exception_value;
            fhigh=-2.0*exception_value;

            if(_chisquared->get_fn(iStart)>_chisquared->target()){
                for(i=0;i<_chisquared->get_dim();i++){
                    trial.set(i,0.5*(get_pt(_centerdex,i)+get_pt(iStart,i)));
                }
                evaluate(trial,&ftrial,&iFound);
                if(iFound>=0 && _chisquared->get_fn(iFound)<_chisquared->target()){
                    iStart=iFound;
                }
                else{
                    return;
                }
            }

            if(sgn>0.0){
                for(i=0;i<_chisquared->get_dim();i++){
                   trial.set(i,get_pt(iStart,i)+dx*sgn*dir.get_data(ix,i));
                }
                evaluate(trial,&ftrial,&iFound);

                if(ftrial<_chisquared->target()){
                    flow=ftrial;
                    for(i=0;i<_chisquared->get_dim();i++){
                        lowball.set(i,trial.get_data(i));
                    }
                }
                else if(ftrial>_chisquared->target()){
                    fhigh=ftrial;
                    for(i=0;i<_chisquared->get_dim();i++){
                        highball.set(i,trial.get_data(i));
                    }
                }
            }

            if(flow>_chisquared->target()){
                flow=_chisquared->get_fn(iStart);
                for(i=0;i<_chisquared->get_dim();i++){
                    lowball.set(i,get_pt(iStart,i));
                }
            }

            while(fhigh<_chisquared->target()){
                for(i=0;i<_chisquared->get_dim();i++){
                    highball.set(i,lowball.get_data(i)+dx*sgn*dir.get_data(ix,i));
                }
                evaluate(highball,&fhigh,&iFound);
                dx+=1.0;
                dx*=2.0;
            }

            iFound=node_bisection(lowball,flow,highball,fhigh,1);
            printf("iFound %d iStart %d ff %e\n",iFound,iStart,_chisquared->get_fn(iFound));
            for(i=0;i<_chisquared->get_dim();i++){
                printf("%.3e ",sgn*dir.get_data(ix,i));
            }
            printf("\n");

            if(iFound>=0){
                _off_center_compass_points.add(iFound);
            }

            if(sgn<0.0 && iFound>=0 && iFound!=iStart){
                dx=0.0;
                for(i=0;i<_chisquared->get_dim();i++){
                    dx+=(get_pt(iStart,i)-get_pt(iFound,i))*dir.get_data(ix,i);
                }
                if(dx<0.0){
                    dx*=-1.0;
                }
            }
            else{
                dx=1.0;
            }
        }
    }

    printf("    done with off-center compass %d -- %d %d -- volume %e\n\n",
    _chisquared->get_called()-ibefore,iStart,_centerdex,volume());

}

double node::apply_quadratic_model(int i1){
    array_1d<double> pt;
    int i;
    pt.set_name("apply_quadratic_model_pt");
    for(i=0;i<_chisquared->get_dim();i++){
        pt.set(i,get_pt(i1,i));
    }
    return apply_quadratic_model(pt);
}

double node::apply_quadratic_model(array_1d<double> &pt){
    is_it_safe("apply_quadratic_model");
    if(_ellipse_center<0){
        printf("WARNING calling apply quadratic model before we have a model\n");
        exit(1);
    }
    int ix,j;
    double ans,mu;
    ans=_chisquared->get_fn(_ellipse_center);
    for(ix=0;ix<_basis_vectors.get_rows();ix++){
        mu=0.0;
        for(j=0;j<_chisquared->get_dim();j++){
            mu+=(pt.get_data(j)-get_pt(_ellipse_center,j))*_basis_vectors.get_data(ix,j);
        }
        ans+=_basis_model.get_data(ix)*mu*mu;
    }
    return ans;

}

int node::find_local_center(){
    if(_geo_centerdex>=0 && _chisquared->get_fn(_geo_centerdex)<=_chisquared->target()){
        return _geo_centerdex;
    }
    else{
        return _centerdex;
    }
}

int node::_are_connected(int i1, int i2){
    int iFound,i;
    double mu;
    array_1d<double> trial;
    trial.set_name("node_are_connected_trial");

    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,0.5*(get_pt(i1,i)+get_pt(i2,i)));
    }

    evaluate(trial,&mu,&iFound);
    if(mu<=_chisquared->target()){
         return 1;
    }
    else{
        return 0;
    }

}


double node::evaluate_dir(int iStart, array_1d<double> &dir_in){

    double flow,fhigh;
    array_1d<double> lowball, highball;
    lowball.set_name("evaluate_dir_lowball");
    highball.set_name("evaluate_dir_highball");

    int i;
    for(i=0;i<_chisquared->get_dim();i++){
        lowball.set(i,get_pt(iStart,i));
        highball.set(i,get_pt(iStart,i));
    }

    flow=ricochet_model(lowball);
    double target;
    if(flow<_chisquared->target()){
        target=_chisquared->target();
    }
    else{
        target=flow+0.1*(_chisquared->target()-_chisquared->chimin());
    }

    int ct;
    double component=1.0;

    double f_max;
    for(i=0;i<_chisquared->get_dim();i++){
        if(fabs(dir_in.get_data(i))>1.0e-20){
            if(i==0 || fabs(dir_in.get_data(i))>f_max){
                if(_max_found.get_data(i)-_min_found.get_data(i)>1.0e-20){
                    f_max=fabs(dir_in.get_data(i));
                    component=0.1*(_max_found.get_data(i)-_min_found.get_data(i));
                }
            }
        }
    }

    fhigh=-2.0*exception_value;
    for(ct=0;fhigh<target && ct<100;ct++){
        for(i=0;i<_chisquared->get_dim();i++){
            highball.add_val(i,component*dir_in.get_data(i));
        }
        fhigh=ricochet_model(highball);
        component*=2.0;
    }

    if(fhigh<target){
        //return a bad value, since we cannot seem to get fhigh to work
        return -2.0*exception_value;
    }

    array_1d<double> trial,best_pt;
    double ftrial,dfbest;
    trial.set_name("evaluate_dir_trial");
    best_pt.set_name("evaluate_dir_best_pt");
    if(fabs(flow-target)<fabs(fhigh-ftrial)){
        dfbest=fabs(flow-target);
        for(i=0;i<_chisquared->get_dim();i++){
            best_pt.set(i,lowball.get_data(i));
        }
    }
    else{
        dfbest=fabs(fhigh-target);
        for(i=0;i<_chisquared->get_dim();i++){
            best_pt.set(i,highball.get_data(i));
        }
    }

    double tol=1.0e-10;
    double distance=node_distance(lowball,highball);
    while(distance>tol){
        for(i=0;i<_chisquared->get_dim();i++){
            trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
        }
        ftrial=ricochet_model(trial);
        if(ftrial<target){
            for(i=0;i<_chisquared->get_dim();i++){
                lowball.set(i,trial.get_data(i));
            }
        }
        else{
            for(i=0;i<_chisquared->get_dim();i++){
                highball.set(i,trial.get_data(i));
            }
        }

        if(fabs(ftrial-target)<dfbest){
            dfbest=fabs(ftrial-target);
            for(i=0;i<_chisquared->get_dim();i++){
                best_pt.set(i,trial.get_data(i));
            }
        }

        distance=node_distance(lowball,highball);
    }

    double sigma;
    //array_1d<int> neigh;
    ftrial=ricochet_model(best_pt,&sigma);
    //return node_distance(neigh.get_data(0),best_pt);
    return sigma;
}

double node::ricochet_model(array_1d<double> &pt, kd_tree &tree,
                            array_1d<int> &neigh_out){
    double sig;
    return _ricochet_model(pt, tree, &sig, 0, neigh_out);
}


double node::ricochet_model(array_1d<double> &pt, kd_tree &tree){
    double sig;
    array_1d<int> nn;
    return _ricochet_model(pt, tree, &sig, 0, nn);
}

double node::ricochet_model(array_1d<double> &pt, kd_tree &tree,
                            double *sig){

    array_1d<int> nn;
    return _ricochet_model(pt, tree, sig, 1, nn);
}

double node::ricochet_model(array_1d<double> &pt){
    array_1d<int> neigh;
    double sig;
    return _ricochet_model(pt,_chisquared->get_tree()[0],&sig,0,neigh);
}

double node::ricochet_model(array_1d<double> &pt, array_1d<int> &neigh){
    double sig;
    return _ricochet_model(pt,_chisquared->get_tree()[0],&sig,0,neigh);
}

double node::ricochet_model(array_1d<double> &pt, double *sig){
    array_1d<int> neigh;
    return _ricochet_model(pt,_chisquared->get_tree()[0],sig,1,neigh);
}

double node::_ricochet_model(array_1d<double> &pt, kd_tree &tree,
                              double *sig, int doSig,
                              array_1d<int> &neigh_out){
    is_it_safe("ricochet_model");

    int i,j,k;
    int npts=2*_chisquared->get_dim();
    double ell;
    array_2d<double> covar,covarin;
    array_1d<int> neigh;
    array_1d<double> dd,min_cache,max_cache,fn_array;

    fn_array.set_name("node_ricochet_model_fn_array");

    min_cache.set_name("node_ricochet_model_min_cache");
    max_cache.set_name("node_ricochet_model_max_cache");
    for(i=0;i<_chisquared->get_dim();i++){
        min_cache.set(i,tree.get_min(i));
        max_cache.set(i,tree.get_max(i));
        tree.set_min(i,_min_found.get_data(i));
        tree.set_max(i,_max_found.get_data(i));
    }

    covar.set_name("node_ricochet_model_covar");
    covarin.set_name("node_ricochet_model_covarin");
    neigh.set_name("node_ricochet_model_neigh");
    dd.set_name("node_ricochet_model_dd");

    array_1d<int> raw_neigh;
    array_1d<double> raw_dd;

    array_1d<double> pt_true;
    pt_true.set_name("node_ricochet_model_pt_true");
    transform_pt_to_truth(pt,pt_true);

    tree.nn_srch(pt_true,2*npts,raw_neigh,raw_dd);

    if(raw_dd.get_data(0)<1.0e-20){
        raw_neigh.remove(0);
        raw_dd.remove(0);
    }

    neigh.add(raw_neigh.get_data(0));
    dd.add(raw_dd.get_data(0));
    double test_dd,test_ddmin;
    int use_it;
    for(i=1;i<raw_neigh.get_dim() && neigh.get_dim()<npts;i++){
        use_it=0;
        if(raw_neigh.get_dim()-i<=npts-neigh.get_dim()){
            use_it=1;
        }
        else{
            test_ddmin=2.0*exception_value;
            for(j=0;j<neigh.get_dim();j++){
                test_dd=node_distance(raw_neigh.get_data(i),neigh.get_data(j));
                if(test_dd<test_ddmin){
                    test_ddmin=test_dd;
                }
            }
            if(test_ddmin>0.001){
                use_it=1;
            }
        }

        if(use_it==1){
            neigh.add(raw_neigh.get_data(i));
            dd.add(raw_dd.get_data(j));
        }
    }

    if(neigh.get_dim()!=npts){
        printf("WARNING in ricochet model neigh %d shldbe %d\n",
        neigh.get_dim(),npts);

        exit(1);
    }
    array_1d<double> mutual_dd,mutual_dd_sorted;
    array_1d<int> mutual_dexes;

    mutual_dd.set_name("node_ricochet_model_mutual_dd");
    mutual_dd_sorted.set_name("node_ricochet_mutual_dd_sorted");
    mutual_dexes.set_name("node_ricochet_mutual_dexes");

    double max_dfn=-1.0;
    for(i=0;i<npts;i++){
        ell=_chisquared->get_fn(neigh.get_data(i))-apply_quadratic_model(neigh.get_data(i));
        if(_chisquared->get_fn(neigh.get_data(i))<exception_value && ell>max_dfn){
            max_dfn=ell;
        }
    }

    if(max_dfn<0.0){
        max_dfn*=(-1.0);
    }

    double mu;
    for(i=0;i<npts;i++){
        ell=_chisquared->get_fn(neigh.get_data(i));
        if(ell<exception_value){
            fn_array.set(i,ell);
        }
        else{
            mu=max_dfn+fabs(apply_quadratic_model(neigh.get_data(i)));
            if(mu>_chisquared->target()){
                fn_array.set(i,mu);
            }
            else{
                fn_array.set(i,max_dfn+_chisquared->target());
            }
        }
    }

    k=0;
    for(i=0;i<npts;i++){
        for(j=i+1;j<npts;j++){
            mutual_dd.set(k,node_distance(neigh.get_data(i),neigh.get_data(j)));
            mutual_dexes.set(k,k);
            k++;
        }
    }

    sort_and_check(mutual_dd,mutual_dd_sorted,mutual_dexes);

    ell=mutual_dd_sorted.get_data(k/2);
    covar.set_dim(npts,npts);
    covarin.set_dim(npts,npts);

    double nugget;
    nugget=1.0e-4;
    for(i=0;i<npts;i++){
        covar.set(i,i,1.0+nugget);
        for(j=i+1;j<npts;j++){
            mu=node_distance(neigh.get_data(i),neigh.get_data(j));
            covar.set(i,j,exp(-0.5*power(mu/ell,2)));
            covar.set(j,i,covar.get_data(i,j));
        }
    }
    invert_lapack(covar,covarin,0);

    array_1d<double> qq;
    qq.set_name("node_ricochet_model_qq");
    for(i=0;i<npts;i++){
        mu=node_distance(neigh.get_data(i),pt);
        qq.set(i,exp(-0.5*power(mu/ell,2)));
    }

    double fbar;
    fbar=apply_quadratic_model(pt);
    /*fbar=0.0;
    for(i=0;i<npts;i++){
        fbar+=_chisquared->get_fn(neigh.get_data(i));
    }
    fbar=fbar/double(npts);*/

    array_1d<double> qbar;
    qbar.set_name("node_ricochet_model_qbar");
    for(i=0;i<npts;i++){
        qbar.set(i,apply_quadratic_model(neigh.get_data(i)));
        //qbar.set(i,fbar);
    }

    mu=fbar;
    for(i=0;i<npts;i++){
        for(j=0;j<npts;j++){
            mu+=qq.get_data(i)*covarin.get_data(i,j)*(fn_array.get_data(j)-qbar.get_data(j));
        }
    }

    double covar_norm;
    int ix,jx;
    if(doSig==1){
        covar_norm=0.0;
        sig[0]=0.0;
        for(i=0;i<npts;i++){
            ix=neigh.get_data(i);
            for(j=0;j<npts;j++){
                jx=neigh.get_data(j);
                covar_norm+=(fn_array.get_data(i)-qbar.get_data(i))*covarin.get_data(i,j)*(fn_array.get_data(j)-qbar.get_data(j));
                sig[0]-=qq.get_data(i)*covarin.get_data(i,j)*qq.get_data(j);
            }
        }
        sig[0]+=1.0+nugget;
        sig[0]*=covar_norm;
    }
    for(i=0;i<npts;i++){
        neigh_out.set(i,neigh.get_data(i));
    }

    for(i=0;i<_chisquared->get_dim();i++){
        tree.set_min(i,min_cache.get_data(i));
        tree.set_max(i,max_cache.get_data(i));
    }

    return mu;
}

void node::add_to_boundary(int dex){
    _boundary_points.add(dex);
}

void node::cull_ricochet(){

    if(_since_culled<_allowed_ricochet_strikes){
        _since_culled++;
        return;
    }

    int i,max_strikes,max_dex;
    max_strikes=-1;
    max_dex=-1;

    array_1d<int> strikes,strikes_sorted,strikes_dexes;
    strikes.set_name("cull_ricochet_strikes");
    strikes_sorted.set_name("cull_ricochet_strikes_sorted");
    strikes_dexes.set_name("cull_ricochet_strikes_dexes");

    for(i=0;i<_ricochet_particles.get_dim();i++){
        if(_ricochet_strikes.get_data(i)>=_allowed_ricochet_strikes){
            strikes.add(-1*_ricochet_strikes.get_data(i));
            strikes_dexes.add(i);
        }
    }

    if(strikes.get_dim()==0){
        return;
    }

    sort_and_check(strikes,strikes_sorted,strikes_dexes);
    int n_started=_ricochet_particles.get_dim();
    int target_particle,j;
    for(i=0;i<strikes_sorted.get_dim() && i<n_started-_chisquared->get_dim()/2;i++){
        target_particle=strikes_dexes.get_data(i);
        remove_particle(strikes_dexes.get_data(i));
        for(j=i+1;j<strikes_sorted.get_dim();j++){
            if(strikes_dexes.get_data(j)>target_particle){
                strikes_dexes.subtract_val(j,1);
            }
        }
    }

    while(_ricochet_particles.get_dim()<_chisquared->get_dim()){
        originate_particle_simplex();
    }

    _since_culled=0;

}

void node::remove_particle(int ip){
    _ricochet_particles.remove(ip);
    _ricochet_origins.remove(ip);
    _ricochet_strikes.remove(ip);
}


void node::set_particle(int ip, int ii){
    if(ii<0){
        return;
    }

    if(_ricochet_particles.get_dim()>ip){
        _ricochet_origins.set(ip,_ricochet_particles.get_data(ip));
    }
    else{
        _ricochet_origins.set(ip,ii);
    }

    _ricochet_particles.set(ip,ii);

    int j;

    increment_ricochet_growth();

    double v1=volume();
    if(v1>1.05*_v0 || ip>=_ricochet_strikes.get_dim()){
        _ricochet_strikes.set(ip,0);

        _v0=v1;
        return;
    }
    else{
        _ricochet_strikes.add_val(ip,1);
    }

    if(_mcmc_acceptances+_mcmc_rejections>_ricochet_particles.get_dim()*10){
        if(_mcmc_acceptances>_mcmc_rejections/3){
            _mcmc_step*=1.25;
        }
        else if(_mcmc_acceptances<_mcmc_rejections/5){
            _mcmc_step*=0.75;
        }
        _mcmc_acceptances=0;
        _mcmc_rejections=0;
    }

    int n_steps=20;

    int ix,i_found;
    array_1d<int> local_associates;
    local_associates.set_name("set_particle_associates");

    double tol=0.05*(_chisquared->target()-_chisquared->chimin());

    for(ix=0;ix<_boundary_points.get_dim();ix++){
        if(_chisquared->get_fn(_associates.get_data(ix))<_chisquared->target()+tol){
            local_associates.add(_boundary_points.get_data(ix));
        }
    }
    if(local_associates.get_dim()==0){
        _v0=v1;
        return;
    }

    mcmc_walk(_ricochet_particles.get_data(ip), &i_found, n_steps, local_associates);
    if(i_found>=0 && i_found!=_ricochet_particles.get_data(ip)){
        _ricochet_origins.set(ip,_ricochet_particles.get_data(ip));
        _ricochet_particles.set(ip,i_found);
    }
    increment_mcmc_growth();

    _v0=volume();
}

void node::initialize_ricochet(){
    is_it_safe("initialize_ricochet");

    if(_found_bases==0){
        find_bases();
        if(_compass_calls==0){
            compass_search();
            compass_diagonal(_centerdex);
        }
        set_transform();
    }

    int i;

    int nParticles=_chisquared->get_dim();
    if(_ricochet_particles.get_dim()>nParticles){
        nParticles=_ricochet_particles.get_dim();
    }

    _ricochet_particles.reset();

    array_1d<int> dexes;
    array_1d<double> dmu;

    int j,ix,iChosen;
    double dist,dist_best,dist_local_best;
    dexes.set_name("node_initialize_ricochet_dexes");
    dmu.set_name("node_initialize_ricochet_dmu");

    array_1d<double> dir;
    dir.set_name("node_initialize_ricochet_dir");

    dir.set_dim(_chisquared->get_dim());

    int iFound,local_center,i_origin;

    local_center=find_local_center();

    int n_0=_ricochet_particles.get_dim();

    while(_ricochet_particles.get_dim()<n_0+nParticles/2){
        originate_particle_shooting();
    }

    while(_ricochet_particles.get_dim()<n_0+nParticles){
        originate_particle_simplex();
    }


    _min_basis_error_changed=0;
    _since_expansion=0;
    _convergence_ct=0;
    _chimin_ricochet=_chimin;
}

int node::mcmc_kick(int iStart, int *iFound, array_1d<double> &dir_out, int max_steps){
    //in this case iStart will be the actual point's index

    double time_before=double(time(NULL));
    double strad,trial_strad;
    double mu,distance_wgt;

    int local_center=find_local_center();

    distance_wgt=_chisquared->target()-_chisquared->chimin();

    int i;
    array_1d<double> min_cache,max_cache;
    for(i=0;i<_chisquared->get_dim();i++){
        min_cache.set(i,_chisquared->get_tree()->get_min(i));
        max_cache.set(i,_chisquared->get_tree()->get_max(i));
        _chisquared->get_tree()->set_min(i,_min_found.get_data(i));
        _chisquared->get_tree()->set_max(i,_max_found.get_data(i));
    }

    array_1d<int> neigh;
    neigh.set_name("node_mcmc_kick_neigh");

    /*mu=ricochet_model(get_pt(iStart)[0],
                      _chisquared->get_tree()[0],
                      neigh);
    strad=mu-distance_wgt*node_distance(iStart,neigh.get_data(0));*/

    double sigma;
    array_1d<double> model_pt;
    model_pt.set_name("node_mcmc_kick_model_pt");
    for(i=0;i<_chisquared->get_dim();i++){
        model_pt.set(i,get_pt(iStart,i));
    }
    ricochet_model(model_pt,&sigma);
    //strad=mu-distance_wgt*sqrt(fabs(sigma));
    strad=fabs(mu-_chisquared->target())-sqrt(fabs(sigma));

    int ix;

    array_1d<double>step,trial,pt;
    step.set_name("node_mcmc_kick_step");
    trial.set_name("node_mcmc_kick_step");
    pt.set_name("node_mcmc_kick_pt");

    int j;
    for(i=0;i<_chisquared->get_dim();i++){
        pt.set(i,get_pt(iStart,i));
    }

    double radius,norm;
    int accept,n_accepted;

    array_1d<double> coulomb_step;
    double rr;
    int ib;
    coulomb_step.set_name("node_mcmc_kick_coulomb_step");
    coulomb_step.set_dim(_chisquared->get_dim());

    n_accepted=0;
    norm=0.1;
    for(ix=0;ix<max_steps;ix++){
        accept=0;
        for(i=0;i<_chisquared->get_dim();i++){
            step.set(i,normal_deviate(_chisquared->get_dice(),0.0,1.0));
        }
        step.normalize();

        if(ix>20){
            coulomb_step.zero();
            for(ib=0;ib<_boundary_points.get_dim();ib++){
                rr=node_distance(_boundary_points.get_data(ib),pt);
                if(rr>1.0e-20){
                    for(i=0;i<_chisquared->get_dim();i++){
                        coulomb_step.add_val(i,(pt.get_data(i)-get_pt(_boundary_points.get_data(ib),i))/(rr*rr));
                    }
                }
            }
            coulomb_step.normalize();

            for(i=0;i<_chisquared->get_dim();i++){
                step.add_val(i,10.0*coulomb_step.get_data(i));
            }
            step.normalize();
        }

        radius=normal_deviate(_chisquared->get_dice(),0.0,1.0);
        for(i=0;i<_chisquared->get_dim();i++){
            trial.set(i,pt.get_data(i)+norm*radius*step.get_data(i));
        }
        mu=ricochet_model(trial,&sigma);
        //trial_strad=mu-distance_wgt*node_distance(neigh.get_data(0),trial);
        //trial_strad=mu-distance_wgt*sqrt(fabs(sigma));
        trial_strad=fabs(mu-_chisquared->target())-sqrt(fabs(sigma));
        if(trial_strad<strad){
            accept=1;
        }
        else{
            norm=_chisquared->random_double();
            if(exp(-0.5*(trial_strad-strad))>norm){
                accept=1;
            }
        }

        if(accept==1){
            n_accepted++;
            for(i=0;i<_chisquared->get_dim();i++){
                pt.set(i,trial.get_data(i));
            }
            strad=trial_strad;
        }

    }

    int iEnd;
    evaluate(pt,&mu,&iEnd);
    double first_found=mu;
    double first_predicted=strad;
    if(iEnd==iStart){
        //printf("mcmc_kick just ended where it started\n");
        for(i=0;i<_chisquared->get_dim();i++){
            _chisquared->get_tree()->set_min(i,min_cache.get_data(i));
            _chisquared->get_tree()->set_max(i,max_cache.get_data(i));
        }
        return 0;
    }
    array_1d<double> dir;
    dir.set_name("node_mcmc_kick_dir");
    double component;
    if(mu<=_chisquared->target()){
        for(i=0;i<_chisquared->get_dim();i++){
            dir.set(i,pt.get_data(i)-get_pt(iStart,i));
        }
        dir.normalize();
        component=1.0;
        while(mu<=_chisquared->target()){
            for(i=0;i<_chisquared->get_dim();i++){
                pt.add_val(i,component*dir.get_data(i));
            }
            evaluate(pt,&mu,&iEnd);
            component*=2.0;
        }
    }

    array_1d<double> e_dir;
    double e_dir_norm,dx,fstart;
    e_dir.set_name("mcmc_kick_e_dir");

    array_1d<double> start;
    start.set_name("node_mcmc_kick_start");
    if(_chisquared->get_fn(iStart)<_chisquared->target()){
        for(i=0;i<_chisquared->get_dim();i++){
           start.set(i,get_pt(iStart,i));
        }
        fstart=_chisquared->get_fn(iStart);
    }
    else{
        for(i=0;i<_chisquared->get_dim();i++){
            e_dir.set(i,get_pt(iStart,i)-get_pt(local_center,i));
        }
        e_dir_norm=e_dir.normalize();
        fstart=_chisquared->get_fn(iStart);
        for(dx=0.1*e_dir_norm;dx<e_dir_norm*1.05 && fstart>_chisquared->target();dx+=0.1*e_dir_norm){
            for(i=0;i<_chisquared->get_dim();i++){
               start.set(i,get_pt(iStart,i)+dx*e_dir.get_data(i));
            }
            evaluate(start,&fstart,&i);
        }
    }


    int iBisection;
    iBisection=node_bisection(start,fstart,pt,mu,0);

    if(iBisection==iStart || iBisection<0){
        //printf("mcmc_kick bisection found starting point\n");
        for(i=0;i<_chisquared->get_dim();i++){
            _chisquared->get_tree()->set_min(i,min_cache.get_data(i));
            _chisquared->get_tree()->set_max(i,max_cache.get_data(i));
        }
        return 0;
    }

    for(i=0;i<_chisquared->get_dim();i++){
        dir_out.set(i,get_pt(iBisection,i)-get_pt(iStart,i));
    }

    iFound[0]=iBisection;
    //printf("mcmc_kick n_accepted %d distance %e -- %e\n",
    //n_accepted,node_distance(iStart,iFound[0]),double(time(NULL))-time_before);
    //printf("first found %e -- %e %d\n",first_found,first_predicted,iFound[0]);

    for(i=0;i<_chisquared->get_dim();i++){
        _chisquared->get_tree()->set_min(i,min_cache.get_data(i));
        _chisquared->get_tree()->set_max(i,max_cache.get_data(i));
    }

    return 1;

}


int node::choose_off_center_point(){
    double target,min_chisq;
    target=0.5*(_chisquared->get_fn(_centerdex)+_chisquared->target());
    min_chisq=0.5*(target+_chisquared->get_fn(_centerdex));

    double dd,ddmax,dd_nn,min_allowable_dd;
    int i,j,use_it,iFound=-1;
    int ipt;

    min_allowable_dd=1.0e-3;

    ddmax=-2.0*exception_value;
    for(i=0;i<_off_center_candidates.get_dim();i++){
        ipt=_off_center_candidates.get_data(i);
        if(_chisquared->get_fn(ipt)<=target && _chisquared->get_fn(ipt)>min_chisq){
            use_it=1;
            dd_nn=2.0*exception_value;
            for(j=0;j<_off_center_origins.get_dim() && use_it==1;j++){
                dd=normalized_node_distance(ipt,_off_center_origins.get_data(j));
                if(dd<min_allowable_dd){
                    use_it=0;
                }
                else if(dd<dd_nn){
                    dd_nn=dd;
                }
            }

            if(use_it==1){
                dd=normalized_node_distance(_centerdex,ipt);
                if(dd<dd_nn){
                    dd_nn=dd;
                }

                if(dd_nn>min_allowable_dd && dd_nn>ddmax){
                    ddmax=dd_nn;
                    iFound=ipt;
                }
            }
        }
        else{
            _off_center_candidates.remove(i);
            i--;
        }
    }

    if(iFound>=0){
        for(i=0;i<_off_center_candidates.get_dim();i++){
            if(_off_center_candidates.get_data(i)==iFound){
                _off_center_candidates.remove(i);
                i--;
            }
        }
        return iFound;
    }

    array_1d<double> lowball,highball,random_dir;

    double flow,fhigh;
    lowball.set_name("node_originate_lowball");
    highball.set_name("node_originate_highball");
    random_dir.set_name("node_originate_random_dir");


    for(i=0;i<_chisquared->get_dim();i++){
        random_dir.set(i,normal_deviate(_chisquared->get_dice(),0.0,1.0));
    }
    random_dir.normalize();

    flow=_chisquared->get_fn(_centerdex);
    for(i=0;i<_chisquared->get_dim();i++){
        lowball.set(i,get_pt(_centerdex,i));
        highball.set(i,get_pt(_centerdex,i));
    }

    fhigh=-2.0*exception_value;
    while(fhigh<=target){
        for(i=0;i<_chisquared->get_dim();i++){
            highball.add_val(i,random_dir.get_data(i));
        }
        evaluate(highball,&fhigh,&iFound);
    }

    iFound=node_bisection(lowball,flow,highball,fhigh,1,target,0.01);
    return iFound;
}

void node::originate_particle_simplex(){
    double chimin0=_chisquared->chimin();

    printf("\norig simplex %d %d -- %e %e\n",
    _chisquared->get_called(),_ricochet_particles.get_dim(),
    _chisquared->target(),_chisquared->get_fn(_centerdex));

    array_1d<double> min,max;
    min.set_name("orig_particle_simplex_min");
    max.set_name("orig_particle_simplex_max");
    array_1d<int> local_associates;
    local_associates.set_name("orig_particle_simplex_local_associates");

    int i,pts_start;
    pts_start=_chisquared->get_pts();
    for(i=0;i<_boundary_points.get_dim();i++){
        local_associates.add(_boundary_points.get_data(i));
    }

    for(i=0;i<_chisquared->get_dim();i++){
        min.set(i,_chisquared->get_min(i));
        max.set(i,_chisquared->get_max(i));
    }

    dchi_boundary_simplex dchi_fn(_chisquared, local_associates);

    simplex_minimizer ffmin;
    ffmin.set_minmax(min,max);
    ffmin.set_chisquared(&dchi_fn);

    array_2d<double> seed;
    seed.set_name("orig_particle_simplex_seed");
    array_1d<double> trial,trial_node;
    trial.set_name("orig_particle_simplex_trial");
    trial_node.set_name("orig_particle_simplex_trial_node");
    double mu,dx,midx;
    double start_min=2.0*exception_value;
    int iFound,i_min;

    while(seed.get_rows()<_chisquared->get_dim()){
        for(i=0;i<_chisquared->get_dim();i++){
            dx=(max.get_data(i)-min.get_data(i));
            midx=0.5*(max.get_data(i)+min.get_data(i));
            trial.set(i,midx+dx*(_chisquared->random_double()-0.5));
        }
        transform_pt_to_node(trial,trial_node);
        evaluate(trial_node,&mu,&iFound);
        if(mu<exception_value && iFound>=0){
            seed.add_row(_chisquared->get_pt(iFound)[0]);
            mu=dchi_fn(_chisquared->get_pt(iFound)[0]);
            if(mu<start_min){
                start_min=mu;
                i_min=iFound;
            }
        }
    }

    double dd,dd_best;
    int i_center;
    i_center=_centerdex;
    dd_best=node_distance(i_min,_centerdex);
    for(i=0;i<_boundary_points.get_dim();i++){
        dd=node_distance(i_min,_boundary_points.get_data(i));
        if(dd<dd_best){
            i_center=_boundary_points.get_data(i);
            dd_best=dd;
        }
    }

    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,0.5*(_chisquared->get_pt(i_center,i)+_chisquared->get_pt(i_min,i)));
    }
    seed.add_row(trial);

    int i_actual;

    trial.reset();
    ffmin.use_gradient();
    ffmin.set_dice(_chisquared->get_dice());
    ffmin.find_minimum(seed,trial);
    transform_pt_to_node(trial,trial_node);
    evaluate(trial_node,&mu,&i_actual);

    int i_walk;
    int n_steps=50;
    double tol=0.1*(_chisquared->target()-_chisquared->chimin());

    double desired=_chisquared->chimin()+2.0*(_chisquared->target()-_chisquared->chimin());

    iFound=i_actual;

    int use_it;

    for(i=_chisquared->get_pts()-1;i>pts_start;i--){
        if(_chisquared->get_fn(i)<_chisquared->target()){
            _associates.add(i);
        }

        if(fabs(_chisquared->get_fn(i)-_chisquared->target())<tol){
            _boundary_points.add(i);
        }

        if(_chisquared->get_fn(i_actual)>desired && _chisquared->get_fn(i)<desired){
            use_it=0;

            if(_chisquared->get_fn(iFound)>desired){
                use_it=1;
            }
            else if(_chisquared->get_fn(iFound)>_chisquared->target()+tol &&
                    _chisquared->get_fn(i)<_chisquared->target()+tol){

                use_it=1;
            }

            if(use_it==1){
                iFound=i;
            }
        }
    }

    if(iFound>=0){
        add_to_log(_log_dchi_simplex, iFound);
    }

    int i_p=-1,i_o=-1;
    array_1d<int> acceptable;
    acceptable.set_name("orig_simplex_acceptable");

    printf("found %e desired %e\n",
    _chisquared->get_fn(iFound),
    desired);

    if(iFound>=0 && _centerdex>=0 && iFound!=_centerdex &&
       _chisquared->get_fn(iFound)<desired){
        for(i=pts_start+1;i<_chisquared->get_pts();i++){
            if(_chisquared->get_fn(i)<_chisquared->target()){
                acceptable.add(i);
            }
        }

        if(acceptable.get_dim()>0){
            for(i=0;i<acceptable.get_dim();i++){
                dd=node_distance(_centerdex,acceptable.get_data(i));
                if(i_p<0 || dd>dd_best){
                    dd_best=dd;
                    i_p=acceptable.get_data(i);
                }
            }
        }

        if(acceptable.get_dim()>1){
            for(i=0;i<acceptable.get_dim();i++){
                if(acceptable.get_data(i)!=i_p){
                    dd=node_distance(acceptable.get_data(i),i_p);
                    if(i_o<0 || dd<dd_best){
                        dd_best=dd;
                        i_o=acceptable.get_data(i);
                    }
                }
            }
        }

        if(i_p>=0 && i_o>=0){
            _ricochet_particles.add(i_p);
            _ricochet_origins.add(i_o);
            _ricochet_strikes.add(0);
        }
        else{
            if(i_p<0){
                i_p=iFound;
            }
            _ricochet_particles.add(i_p);
            _ricochet_origins.add(_centerdex);
            _ricochet_strikes.add(0);

            local_associates.reset();
            for(i=0;i<_boundary_points.get_dim();i++){
                if(_chisquared->get_fn(_associates.get_data(i))<_chisquared->target()+tol){
                    local_associates.add(_boundary_points.get_data(i));
                }
            }
            mcmc_walk(_ricochet_particles.get_data(_ricochet_particles.get_dim()-1),
                      &i_walk, n_steps, local_associates);

            if(i_walk!=i_p){
                _ricochet_origins.set(_ricochet_origins.get_dim()-1,i_p);
                _ricochet_particles.set(_ricochet_particles.get_dim()-1,i_walk);
            }
        }
    }
    else{
        originate_particle_shooting();
    }

    if(_chisquared->chimin()<chimin0){
        _centerdex=_chisquared->mindex();
        _chimin=_chisquared->chimin();
    }

    recalibrate_max_min();

    printf("simplex found %e %e from %e -- %e %e %e\n",_chisquared->get_fn(iFound),
    dchi_fn(_chisquared->get_pt(iFound)[0]),start_min,
    _chisquared->chimin(),_chimin,_chisquared->get_fn(_centerdex));
    /*printf("    min pt\n");
    for(i=0;i<_chisquared->get_dim();i++){
        printf("    %.3e\n",_chisquared->get_pt(_centerdex,i));
    }*/

    increment_simplex_growth();

}

void node::originate_particle_shooting(){
    int iFound, i_origin;
    iFound=_originate_particle_shooting(&i_origin);
    if(iFound>=0 && i_origin>=0){
        _ricochet_particles.add(iFound);
        _ricochet_origins.add(i_origin);
        _ricochet_strikes.add(0);
    }
}

int node::_originate_particle_shooting(int *origin_out){
    array_1d<double> dir;
    dir.set_name("originate_particle_shooting_dir");
    int iOrigin=choose_off_center_point();
    _off_center_origins.add(iOrigin);

    int i,j;

    double norm,rr;
    array_1d<double> sub_dir,pointing_dir;
    sub_dir.set_name("node_originate_sub_dir");
    pointing_dir.set_name("node_originate_pointing_dir");
    for(i=0;i<_chisquared->get_dim();i++){
        pointing_dir.set(i,0.0);
    }
    for(i=0;i<_boundary_points.get_dim();i++){
        for(j=0;j<_chisquared->get_dim();j++){
            sub_dir.set(j,get_pt(iOrigin,j)-get_pt(_boundary_points.get_data(i),j));
        }
        norm=sub_dir.normalize();
        rr=node_distance(iOrigin,_boundary_points.get_data(i));
        if(norm>1.0e-10){
            for(j=0;j<_chisquared->get_dim();j++){
                pointing_dir.add_val(j,sub_dir.get_data(j)/rr);
            }
        }
    }

    norm=pointing_dir.normalize();
    double dotproduct=-1.0;

    for(i=0;i<_chisquared->get_dim();i++){
        dir.set(i,normal_deviate(_chisquared->get_dice(),0.0,1.0));
    }
    dir.normalize();

    dotproduct=0.0;
    for(i=0;i<_chisquared->get_dim();i++){
        dotproduct+=dir.get_data(i)*pointing_dir.get_data(i);
    }

    while(dotproduct<0.0 && norm>1.0e-10){
        for(i=0;i<_chisquared->get_dim();i++){
            dir.set(i,normal_deviate(_chisquared->get_dice(),0.0,1.0));
        }
        dir.normalize();

        dotproduct=0.0;
        for(i=0;i<_chisquared->get_dim();i++){
            dotproduct+=dir.get_data(i)*pointing_dir.get_data(i);
        }
    }

    int iFound;
    iFound=node_bisection_origin_dir(iOrigin,dir);

    origin_out[0]=iOrigin;
    increment_compass_growth();
    return iFound;

}


void node::search(){

    printf("calling node search\n");


    if(_chisquared->get_fn(_centerdex)>_chisquared->target()){
        simplex_search();
        if(_chisquared->get_fn(_centerdex)>_chisquared->target()){
            _active=0;
            return;
        }
    }

    if(_ct_simplex<=_ct_ricochet &&
        _failed_simplexes<3 &&
        _chisquared->could_it_go_lower(_chimin)>0 &&
        _do_simplex==1){

        simplex_search();
    }

    cull_ricochet();

    double projectedVolume0=projected_volume();
    double volume0=volume();
    double target0=_chisquared->target();
    int ibefore=_chisquared->get_called();

    ricochet();

    double tol;
    tol=0.1*(_chisquared->target()-_chisquared->get_fn(_centerdex_basis));

    if(_chisquared->get_fn(_centerdex)<_chisquared->get_fn(_centerdex_basis)-tol){
        find_bases();
        compass_search();
        compass_diagonal(_centerdex);
        set_transform();
        increment_compass_growth();
    }

    if(fabs(_chisquared->target()-target0)>0.01){
        recalibrate_max_min();
    }

    double volume1=volume();
    double projectedVolume1=projected_volume();

    if(fabs(volume0-volume1)>(0.01*_chisquared->get_dim())*volume0){
        _since_expansion=0;
        _convergence_ct=0;
    }
    else{
        _convergence_ct++;
    }

    if(_convergence_ct>_allowed_ricochet_strikes){
        printf("deactivating because we did not expand\n");
        _active=0;
    }

    if(_ricochet_particles.get_dim()==0){
        initialize_ricochet();
    }

    array_1d<double> dir;
    dir.set_name("node_search_dir");
    int iParticle,i,j,n_particles;

    if(_active==0){
        target0=_chisquared->target();
        volume0=volume();
        projectedVolume0=projected_volume();

        find_bases();
        set_geo_center();
        increment_compass_growth();
        n_particles=_ricochet_particles.get_dim();
        _ricochet_particles.reset();
        _ricochet_strikes.reset();
        _ricochet_origins.reset();

        while(_ricochet_particles.get_dim()<n_particles){
            originate_particle_simplex();
        }

        volume1=volume();
        projectedVolume1=projected_volume();

        _active=1;
        _convergence_ct=0;
    }

    if(_ricochet_particles.get_dim()==0){
        _active=0;
    }

    printf("failed kicks %d successful kicks %d\n",
    _failed_kicks,_successful_kicks);

}


void node::simplex_search(){
    printf("    node simplex search %e\n",_chisquared->get_fn(_centerdex));

    is_it_safe("simplex_search");

    int ibefore=_chisquared->get_called();
    int is_acceptable;
    int i,j,iFound;
    double target,tolerance;

    array_1d<double> trial,simplex_min,simplex_max;
    array_1d<int> found_dexes;
    trial.set_name("node_simplex_search_trial");
    simplex_min.set_name("node_simplex_min");
    simplex_max.set_name("node_simplex_max");
    found_dexes.set_name("found_dexes");

    is_acceptable=0;
    while(is_acceptable==0){
        if(_true_min.get_dim()==_chisquared->get_dim() && _true_max.get_dim()==_chisquared->get_dim()){
            for(i=0;i<_chisquared->get_dim();i++){
                simplex_min.set(i,_true_min.get_data(i));
                simplex_max.set(i,_true_max.get_data(i));
            }
        }
        is_acceptable=1;
        if(simplex_min.get_dim()!=_chisquared->get_dim() || simplex_max.get_dim()!=_chisquared->get_dim()){
            is_acceptable=0;
            printf("not acceptable because of dim %d %d\n",
            simplex_min.get_dim(),simplex_max.get_dim());
        }

        for(i=0;i<_chisquared->get_dim() && is_acceptable==1;i++){
            if(simplex_max.get_data(i)-simplex_min.get_data(i)<1.0e-30){
                is_acceptable=0;
                printf("not acceptable because %d %.3e %.3e %.3e -- %.3e\n",
                i,simplex_min.get_data(i),simplex_max.get_data(i),
                simplex_max.get_data(i)-simplex_min.get_data(i),
                _true_max.get_data(i)-_true_min.get_data(i));
            }
        }

        if(is_acceptable==0){
            for(i=0;i<_chisquared->get_dim();i++){
                trial.set(i,normal_deviate(_chisquared->get_dice(),0.0,1.0));
            }
            trial.normalize();
            target=_chisquared->get_fn(_centerdex)+(_chisquared->target()-_chisquared->chimin());
            tolerance = 0.01*(_chisquared->target()-_chisquared->chimin());
            iFound=node_bisection_origin_dir(_centerdex, trial, target, tolerance);
        }
    }

    int center0 = _centerdex;

    array_1d<int> dexes;
    array_1d<double> minpt,dmu,dmusorted;
    array_2d<double> seed;
    minpt.set_name("node_simplex_search_minpt");
    seed.set_name("node_simplex_search_minpt");
    dexes.set_name("node_simplex_search_dexes");
    dmu.set_name("node_simplex_search_dmu");
    dmusorted.set_name("node_simplex_search_dmusorted");


    seed.set_cols(_chisquared->get_dim());
    double mu;

    array_1d<double> pt_node;
    pt_node.set_name("node_simplex_pt_node");

    while(seed.get_rows()<_chisquared->get_dim()+1){
        for(i=0;i<_chisquared->get_dim();i++){
            trial.set(i,normal_deviate(_chisquared->get_dice(),0.0,1.0));
        }
        trial.normalize();
        mu=_chisquared->random_double();
        for(i=0;i<_chisquared->get_dim();i++){
            trial.multiply_val(i,mu*(simplex_max.get_data(i)-simplex_min.get_data(i)));
            trial.add_val(i,0.5*(simplex_min.get_data(i)+simplex_max.get_data(i)));
        }

        transform_pt_to_node(trial, pt_node);

        evaluate(pt_node,&mu,&iFound);
        if(iFound>=0){
            seed.add_row(trial);
        }
    }

    if(seed.get_rows()!=_chisquared->get_dim()+1){
        printf("WARNING in node_simplex_search seed has %d rows want %d\n",
        seed.get_rows(),_chisquared->get_dim()+1);

        exit(1);
    }
    simplex_minimizer ffmin;
    ffmin.set_minmax(simplex_min,simplex_max);
    ffmin.set_chisquared(_chisquared);
    ffmin.set_dice(_chisquared->get_dice());
    ffmin.find_minimum(seed,minpt);

    transform_pt_to_node(minpt, pt_node);

    evaluate(pt_node,&mu,&i);

    _ct_simplex+=_chisquared->get_called()-ibefore;
    if(_centerdex==center0){
        _failed_simplexes++;
    }
    else{
        _failed_simplexes=0;
    }

}


int node::_ricochet(int iparticle){
    if(_ricochet_particles.get_data(iparticle)==_ricochet_origins.get_data(iparticle)){
        return _ricochet_particles.get_data(iparticle);
    }

    if(_chisquared->get_fn(_ricochet_particles.get_data(iparticle))>=_chisquared->target()){
        _shift_ricochet(iparticle);
    }

    array_1d<double> gradient,dir_0;
    gradient.set_name("_ricochet_gradient");
    dir_0.set_name("_ricochet_dir_0");

    double reflection_coeff,component,gnorm,dirnorm;

    array_1d<double> dir;
    dir.set_name("_ricochet_dir");
    int i,j,iFound;

    node_gradient(_ricochet_particles.get_data(iparticle),gradient);

    reflection_coeff=2.0;

    gnorm=gradient.normalize();
    component=0.0;
    for(i=0;i<_chisquared->get_dim();i++){
        dir_0.set(i,get_pt(_ricochet_particles.get_data(iparticle),i)-get_pt(_ricochet_origins.get_data(iparticle),i));
    }
    dir_0.normalize();

    for(i=0;i<_chisquared->get_dim();i++){
        component+=dir_0.get_data(i)*gradient.get_data(i);
    }

    for(i=0;i<_chisquared->get_dim();i++){
        dir.set(i,dir_0.get_data(i)-reflection_coeff*component*gradient.get_data(i));
    }

    dirnorm=dir.normalize();

    iFound=node_bisection_origin_dir(_ricochet_particles.get_data(iparticle), dir);

    if(iFound>=0){
        add_to_log(_log_ricochet, iFound);
    }

    return iFound;

}


void node::_shift_ricochet(int ix){

    if(_chisquared->get_fn(_ricochet_particles.get_data(ix))<_chisquared->target()){
        return;
    }

    array_1d<double> trial,dir;
    trial.set_name("shift_ricochet_trial");
    dir.set_name("shift_ricochet_dir");

    int i,iFound;
    double mu;
    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,0.5*(get_pt(_ricochet_particles.get_data(ix),i)+
                         get_pt(_ricochet_origins.get_data(ix),i)));
    }
    evaluate(trial,&mu,&iFound);
    if(iFound>=0 && mu<_chisquared->target()){
        for(i=0;i<_chisquared->get_dim();i++){
            dir.set(i,get_pt(_ricochet_particles.get_data(ix),i)-trial.get_data(i));
        }
        iFound=node_bisection_origin_dir(iFound,dir);
        if(iFound>=0 && _chisquared->get_fn(iFound)<_chisquared->target()){
            _ricochet_particles.set(ix,iFound);
            return;
        }
    }

    for(i=0;i<_chisquared->get_dim();i++){
        dir.set(i,get_pt(_ricochet_particles.get_data(ix),i)-get_pt(_centerdex,i));
    }

    dir.normalize();
    iFound=node_bisection_origin_dir(_centerdex, dir);
    if(iFound>=0 && _chisquared->get_fn(iFound)<_chisquared->target()){
        _ricochet_particles.set(ix,iFound);
    }

}

void node::ricochet(){
    printf("starting ricochet\n");
    is_it_safe("ricochet");

    if(_ricochet_particles.get_dim()==0){
        initialize_ricochet();
    }

   swarm_search();

   double volume0=volume();
   double projectedVolume0=projected_volume();

   int ibefore=_chisquared->get_called();
   int rcalls_before=_bisection_calls;
   int rbefore=_total_bisections;

   double dd_max=-1.0;
   double dd_min=2.0*exception_value;

   printf("    starting ricochet with volume %e and pts %d\n",volume(),
   _ricochet_particles.get_dim());

   int ix,iFound,i,j;

   double dx,x1,x2,y1,y2,distanceMin;
   double gnorm,dirnorm;
   int local_pts0;

    int kicked,local_center,is_connected,highball_call_0;
    double reflection_coeff;
    local_center=find_local_center();

    int abort_kicks;

   array_1d<double> scratch;
   scratch.set_name("node_ricochet_scratch");

   int i_origin;
   double dd;

   double v0,v1;

   distanceMin=1.0e-2;
   for(ix=0;ix<_ricochet_particles.get_dim();ix++){
       _v0=volume();
       v0=volume();
       local_pts0=_chisquared->get_pts();

       while(_chisquared->get_fn(_ricochet_particles.get_data(ix))>_chisquared->target()){
           _shift_ricochet(ix);
       }

       i_origin=_ricochet_particles.get_data(ix);

       _proper_ricochets++;
       iFound=_ricochet(ix);

       if(iFound>=0){
           set_particle(ix,iFound);
       }
       else{
           iFound=_ricochet_particles.get_data(ix);
           set_particle(ix,iFound);

       }

       dd=node_distance(i_origin,_ricochet_particles.get_data(ix));
       if(dd<dd_min)dd_min=dd;
       if(dd>dd_max)dd_max=dd;

   }

   printf("    done with actual ricochet: volume %e dd %e %e\n",
   volume(),dd_min,dd_max);

   double ricochet_dd,ddmax;
   int iChosen;

    _total_ricochets+=_ricochet_particles.get_dim();

   int r_called=_chisquared->get_called()-ibefore;

   _ricochet_calls+=_chisquared->get_called()-ibefore;
   _ricochet_bisection_calls+=_bisection_calls-rcalls_before;
   _ricochet_bisections+=_total_bisections-rbefore;


    _ct_ricochet+=_chisquared->get_called()-ibefore;

    double volume1=volume();
    double projectedVolume1=projected_volume();

    printf("    ending ricochet with volume %e from %e -- %d -- %d\n\n",
    volume(),volume0,r_called,_ricochet_particles.get_dim());
}

void node::mcmc_walk(int i_start, int *i_found, int n_steps,
                     array_1d<int> &local_associates){

    array_1d<double> step,trial;
    step.set_name("mcmc_walk_step");
    trial.set_name("mcmc_walk_trial");

    double temp_term,roll;
    double delta,dmu,dd,ddmin;
    int i_pt=i_start;

    int i;
    ddmin=2.0*exception_value;
    for(i=0;i<local_associates.get_dim();i++){
        dd=normalized_node_distance(local_associates.get_data(i),i_pt);
        if(dd<ddmin){
            ddmin=dd;
        }
    }

    dmu=fabs(_chisquared->get_fn(i_pt)-_chisquared->target());
    delta=_chisquared->target()-_chisquared->chimin();
    double chi_new,chi_old,exp_term;
    if(_chisquared->target()<_chisquared->get_fn(i_pt)){
        exp_term=exp(-0.1*dmu/delta);
    }
    else{
        exp_term=1.0;
    }
    chi_old=dmu-exp_term*delta*ddmin*2.0;

    int i_step,local_acceptances,accept_it;
    double mu;
    int iFound,j;

    local_acceptances=0;
    for(i_step=0;i_step<n_steps;i_step++){

        if(i_step<n_steps/2){
            temp_term=exp(-1.0*(n_steps/2-i_step));
        }
        else{
            temp_term=1.0;
        }

        delta=_chisquared->target()-_chisquared->chimin();

        for(i=0;i<_chisquared->get_dim();i++){
            step.set(i,0.0);
        }

        for(i=0;i<_chisquared->get_dim();i++){
            mu=get_projected_norm(i)*normal_deviate(_chisquared->get_dice(),0.0,1.0);
            for(j=0;j<_chisquared->get_dim();j++){
                step.add_val(j,_mcmc_step*mu*_basis_vectors.get_data(i,j));
            }
        }

        for(i=0;i<_chisquared->get_dim();i++){
            trial.set(i,get_pt(i_pt,i)+step.get_data(i));
        }
        evaluate(trial,&mu,&iFound);
        accept_it=0;
        if(iFound>=0){
            ddmin=2.0*exception_value;
            for(i=0;i<local_associates.get_dim();i++){
                dd=normalized_node_distance(local_associates.get_data(i),iFound);
                if(dd<ddmin){
                    ddmin=dd;
                }
            }
            dmu=fabs(mu-_chisquared->target());

            if(_chisquared->target()<mu){
                exp_term=exp(-0.1*dmu/delta);
            }
            else{
                exp_term=1.0;
            }

            chi_new=dmu-exp_term*delta*ddmin*2.0;

            if(chi_new<chi_old){
                accept_it=1;
            }
            else{
                roll=_chisquared->random_double()*temp_term;
                if(exp(-0.5*(chi_new-chi_old))>roll){
                    accept_it=1;
                }
            }
        }

        if(temp_term>0.95){
            if(accept_it==0){
                _mcmc_rejections++;
            }
            else{
                _mcmc_acceptances++;
            }
        }

        if(accept_it==1){
            local_acceptances++;
            i_pt=iFound;
            chi_old=chi_new;
        }

    }

    array_1d<double> lowball,highball;
    double flow, fhigh;
    highball.set_name("mcmc_walk_highball");
    lowball.set_name("mcmc_walk_lowball");
    int i_nearest;
    double dd_nearest;

    if(i_pt>=0 && _chisquared->get_fn(i_pt)>_chisquared->target()){
        i_nearest=-1;
        for(i=0;i<_associates.get_dim();i++){
            if(_chisquared->get_fn(_associates.get_data(i))<_chisquared->target()){
                mu=node_distance(i_pt,_associates.get_data(i));
                if(i_nearest<0 || mu<dd_nearest){
                    i_nearest=_associates.get_data(i);
                    dd_nearest=mu;
                }
            }
        }

        if(i_nearest>=0){
            for(i=0;i<_chisquared->get_dim();i++){
                lowball.set(i,get_pt(i_nearest,i));
                highball.set(i,get_pt(i_pt,i));
            }
            flow=_chisquared->get_fn(i_nearest);
            fhigh=_chisquared->get_fn(i_pt);
            i_pt=node_bisection(lowball,flow,highball,fhigh,1);
        }

    }


    i_found[0]=i_pt;

    if(i_pt>=0){
        add_to_log(_log_mcmc, i_pt);
    }

    printf("    local_acceptances %d %e\n",local_acceptances,volume());
}


int node::get_highball_calls(){
    return _highball_calls;
}

int node::get_gradient_calls(){
    return _gradient_calls;
}

int node::get_total_bisections(){
    return _total_bisections;
}

int node::get_bisection_calls(){
    return _bisection_calls;
}

int node::get_proper_ricochets(){
    return _proper_ricochets;
}

int node::get_total_ricochets(){
    return _total_ricochets;
}

int node::get_ricochet_calls(){
    return _ricochet_calls;
}

int node::get_ricochet_bisection_calls(){
    return _ricochet_bisection_calls;
}

int node::get_ricochet_bisections(){
    return _ricochet_bisections;
}

int node::get_n_particles(){
    return _ricochet_particles.get_dim();
}

int node::get_swarm_outside(){
    return _swarm_outsiders;
}

int node::get_swarm_expand(){
    return _swarm_expanders;
}

void node::swarm_shoot(int i_start){

    array_1d<double> gradient;
    gradient.set_name("swarm_shoot_gradient");
    node_gradient(i_start, gradient);
    gradient.normalize();
    int iFound;
    iFound=node_bisection_origin_dir(i_start,gradient);
    if(iFound<0){
        printf("WARNING swarm shot failed\n");
        exit(1);
    }

    if(iFound>=0 && i_start>=0){
        add_to_log(_log_swarm, iFound);
        _ricochet_particles.add(iFound);
        _ricochet_origins.add(i_start);
        _ricochet_strikes.add(0);
    }

}

void node::swarm_evaluate(array_1d<double> &pt, double *mu){
    int i,j;
    double v0=volume();
    array_1d<double> buffer;
    buffer.set_name("swarm_evaluate_buffer");
    for(i=0;i<_chisquared->get_dim();i++){
        buffer.set(i,_swarm_center.get_data(i));
    }
    for(i=0;i<_chisquared->get_dim();i++){
        for(j=0;j<_chisquared->get_dim();j++){
            buffer.add_val(j,_basis_vectors.get_data(i,j)*pt.get_data(i)*_swarm_norm.get_data(i));
        }
    }
    int iFound;
    evaluate(buffer,mu,&iFound);

    if(mu[0]<_chisquared->target()){
        add_to_log(_log_swarm, iFound);
    }

    double v1=volume();
    if(v1>v0*1.05){
        swarm_shoot(iFound);
        _swarm_expanders++;
    }

    int outside=0;
    if(mu[0]<_chisquared->target()){
        for(i=0;i<_chisquared->get_dim() && outside==0;i++){
            if(buffer.get_data(i)<_swarm_min.get_data(i)){
                outside=1;
            }
            if(buffer.get_data(i)>_swarm_max.get_data(i)){
                outside=1;
            }
        }
        if(outside==1){
            _swarm_outsiders++;
        }
    }

    double dmu=fabs(mu[0]-_chisquared->target());
    double chisq=mu[0];
    double deltachisq=_chisquared->target()-_chisquared->chimin();
    double dd,ddmin;

    if(iFound>=0){
        ddmin=2.0*exception_value;
        for(i=0;i<_swarm_associates.get_dim();i++){
            dd=normalized_node_distance(iFound,_swarm_associates.get_data(i));
            if(dd<ddmin){
                ddmin=dd;
            }
        }
    }
    else{
        ddmin=0.0;
    }

    double exp_term;
    if(_chisquared->target()<mu[0]){
        exp_term=exp(-0.1*dmu/deltachisq);
    }
    else{
        exp_term=1.0;
    }

    mu[0]=dmu-exp_term*deltachisq*ddmin*2.0;

}

void node::swarm_search(){

    int i;
    double tol=0.05*(_chisquared->target()-_chisquared->chimin());
    _swarm_associates.reset();
    for(i=0;i<_boundary_points.get_dim();i++){
        if(_chisquared->get_fn(_boundary_points.get_data(i))<_chisquared->target()+tol){
            _swarm_associates.add(_boundary_points.get_data(i));
        }
    }
    if(_swarm_associates.get_dim()==0){
        increment_swarm_growth();
        return;
    }

    if(_swarm_acceptances+_swarm_rejections>10*_swarm.get_rows()){
        if(_swarm_acceptances<(_swarm_rejections)/5){
            _swarm_step*=0.75;
        }
        else if(_swarm_acceptances>(_swarm_rejections)/3){
            _swarm_step*=1.25;
        }
        _swarm_acceptances=0;
        _swarm_rejections=0;
    }

    array_1d<double> pp,vv,max,min;

    pp.set_name("swarm_pp");
    vv.set_name("swarm_vv");
    max.set_name("swarm_norm");
    min.set_name("swarm_norm");

    int j;
    for(i=0;i<_chisquared->get_dim();i++){
        _swarm_center.set(i,0.0);
    }

    for(i=0;i<_ricochet_particles.get_dim();i++){
        for(j=0;j<_chisquared->get_dim();j++){
            vv.set(j,get_pt(_ricochet_particles.get_data(i),j));
            _swarm_center.add_val(j,get_pt(_ricochet_particles.get_data(i),j));
        }
        project_to_bases(vv,pp);
        for(j=0;j<_chisquared->get_dim();j++){
            if(j>=min.get_dim() || pp.get_data(j)<min.get_data(j)){
                min.set(j,pp.get_data(j));
            }
            if(j>=max.get_dim() || pp.get_data(j)>max.get_data(j)){
                max.set(j,pp.get_data(j));
            }
        }
    }

    for(i=0;i<_chisquared->get_dim();i++){
        _swarm_max.set(i,max.get_data(i));
        _swarm_min.set(i,min.get_data(i));
        _swarm_norm.set(i,max.get_data(i)-min.get_data(i));
        _swarm_center.divide_val(i,double(_ricochet_particles.get_dim()));
    }

    double mu;
    int i_found;
    evaluate(_swarm_center,&mu,&i_found);
    _f_swarm_center=mu;

    array_1d<double> step,trial;
    step.set_name("swarm_step");
    trial.set_name("swarm_trial");

    array_1d<double> chi_old;
    chi_old.set_name("swarm_chi_old");

    int ipt;

    if(_swarm.get_rows()==0){
        _swarm.set_cols(_chisquared->get_dim());
        for(ipt=0;ipt<_chisquared->get_dim();ipt++){
            for(i=0;i<_chisquared->get_dim();i++){
                trial.set(i,2.0*(_chisquared->random_double()-0.5));
                _swarm.set(ipt,i,trial.get_data(i));
            }
            swarm_evaluate(trial,&mu);
            chi_old.set(ipt,mu);

            if(mu>exception_value){
                ipt--;
            }
        }
    }
    else{
        for(ipt=0;ipt<_swarm.get_rows();ipt++){
            swarm_evaluate(_swarm(ipt)[0],&mu);
            chi_old.set(ipt,mu);
        }
    }

    int n_steps=20;
    int i_step;
    double rr;
    int accept_it;
    double temp_term;
    for(i_step=0;i_step<n_steps;i_step++){

        if(i_step<n_steps/2){
            temp_term=exp(-0.5*(n_steps/2-i_step));
        }
        else{
            temp_term=1.0;
        }

        for(ipt=0;ipt<_swarm.get_rows();ipt++){

            rr=normal_deviate(_chisquared->get_dice(),_swarm_step,0.2*_swarm_step);

            for(i=0;i<_chisquared->get_dim();i++){
                step.set(i,normal_deviate(_chisquared->get_dice(),0.0,1.0));
            }

            step.normalize();

            for(i=0;i<_chisquared->get_dim();i++){
                trial.set(i,_swarm.get_data(ipt,i)+rr*step.get_data(i));
            }

            swarm_evaluate(trial,&mu);

            accept_it=0;

            if(mu<chi_old.get_data(ipt)){
                accept_it=1;
            }
            else{
                rr=_chisquared->random_double()*temp_term;
                if(exp(-0.5*(mu-chi_old.get_data(ipt)))>rr){
                    accept_it=1;
                }
            }

            if(accept_it==1){
                chi_old.set(ipt,mu);
                for(i=0;i<_chisquared->get_dim();i++){
                    _swarm.set(ipt,i,trial.get_data(i));
                }
            }

            if(temp_term>0.95){
                if(accept_it==1){
                    _swarm_acceptances++;
                }
                else{
                    _swarm_rejections++;
                }
            }
        }
    }

    increment_swarm_growth();

}


///////////////arrayOfNodes code below//////////

arrayOfNodes::arrayOfNodes(){
    _data=NULL;
    _ct=0;
    _room=0;
}

arrayOfNodes::~arrayOfNodes(){
    if(_data!=NULL){
        delete [] _data;
    }
}

void arrayOfNodes::add(int cc, chisq_wrapper *gg, asymm_array_2d<int> *a){

    node *buffer;
    int i,j;

    if(_ct==_room){
        if(_ct>0){
            buffer=new node[_ct];
            for(i=0;i<_ct;i++){
                buffer[i].copy(_data[i]);
            }

            delete [] _data;
        }

        _room+=5;
        _data=new node[_room];

        if(_ct>0){
            for(i=0;i<_ct;i++){
                _data[i].copy(buffer[i]);
            }
            delete [] buffer;
        }
    }

    _data[_ct].set_chisquared(gg);
    _data[_ct].set_center(cc);
    _data[_ct].set_id_dex(_ct);
    _data[_ct].set_log(a);
    printf("time to initialize ricochet\n");
    _data[_ct].initialize_ricochet();
    printf("done with that\n");
    _ct++;

}

void arrayOfNodes::add(int i, chisq_wrapper *g){
    add(i,g,NULL);
}

void arrayOfNodes::add(chisq_wrapper *g, int i){
    add(i,g,NULL);
}

int arrayOfNodes::get_dim(){
    return _ct;
}

void arrayOfNodes::cull(){
    if(_ct<2){
        return;
    }

    int i,j,kill_it;
    for(i=0;i<_ct;i++){
        if(_data[i].get_activity()==1){
            for(j=i+1;j<_ct;j++){
                if(_data[j].get_activity()==1){
                    kill_it=_data[i].is_this_an_associate(_data[j].get_center());
                    if(kill_it==1){
                        _data[i].merge(_data[j]);
                        _data[j].deactivate();
                    }
                }
            }
        }
    }
}

void arrayOfNodes::remove(int ii){
    int i;
    if(ii>=_ct) return;

    for(i=ii+1;i<_ct;i++){
        _data[i-1].copy(_data[i]);
    }
    _ct--;
}

node* arrayOfNodes::operator()(int dex){
    if(dex<0 || dex>=_ct){
        printf("WARNING asked arrayOfNodes for dex %d but only have %d (operator)\n",
        dex,_ct);

        exit(1);
    }
    return &_data[dex];
}
