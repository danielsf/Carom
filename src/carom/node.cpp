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

void node::initialize(){
    _chisquared=NULL;
    _chimin=2.0*exception_value;
    _centerdex=-1;
    _centerdex_basis=-1;
    _min_changed=0;
    _active=1;
    _found_bases=0;
    _ct_ricochet=0;
    _ct_simplex=0;
    _calls_to_ricochet=0;
    _allowed_ricochet_strikes=3;
    _since_expansion=0;
    _min_basis_error=exception_value;
    _min_basis_error_changed=0;
    _failed_simplexes=0;
    _failed_kicks=0;
    _successful_kicks=0;
    _volume_of_last_geom=0.0;
    
    _compass_points.set_name("node_compass_points");
    _ricochet_candidates.set_name("node_ricochet_candidates");
    _ricochet_candidate_velocities.set_name("node_ricochet_candidate_velocities");
    _off_center_compass_points.set_name("node_off_center_compass_points");
    _off_center_origins.set_name("node_off_center_origins");
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
    _distance_traveled.set_name("node_distance_traveled");
    _ricochet_particles.set_name("node_ricochet_particles");
    _ricochet_origins.set_name("node_ricochet_origins");
    _ricochet_velocities.set_name("node_ricochet_velocities");
    _ricochet_strikes.set_name("node_ricochet_strikes");
    _ricochet_discoveries.set_name("node_ricochet_discoveries");
    _ricochet_discovery_time.set_name("node_ricochet_discovery_time");
    _ricochet_discovery_dexes.set_name("node_ricochet_discovery_dexes");
    _ricochet_grad_norm.set_name("node_ricochet_grad_norm");
    _ricochet_dir_norm.set_name("node_ricochet_dir_norm");
    _ricochet_distances.set_name("node_ricochet_distances");
    _ricochet_mu.set_name("node_ricochet_mu");
    _ricochet_strike_log.set_name("node_ricochet_strike_log");
    _associates.set_name("node_associates");
    _boundary_points.set_name("node_boundary_points");
}

void node::copy(const node &in){
    if(this==&in){
        return;
    }
    _centerdex=in._centerdex;
    _centerdex_basis=in._centerdex_basis;
    _chimin=in._chimin;
    _chimin_bases=in._chimin_bases;
    _min_changed=in._min_changed;
    _active=in._active;
    _failed_simplexes=in._failed_simplexes;
    _found_bases=in._found_bases;
    _ct_ricochet=in._ct_ricochet;
    _ct_simplex=in._ct_simplex;
    _calls_to_ricochet=in._calls_to_ricochet;
    _allowed_ricochet_strikes=in._allowed_ricochet_strikes;
    _ellipse_center=in._ellipse_center;
    _since_expansion=in._since_expansion;
    _min_basis_error=in._min_basis_error;
    _min_basis_error_changed=in._min_basis_error_changed;
    _volume_of_last_geom=in._volume_of_last_geom;
    
    int i,j;
    
    _chisquared=in._chisquared;

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
    
    _ricochet_candidates.reset();
    for(i=0;i<in._ricochet_candidates.get_dim();i++){
        _ricochet_candidates.set(i,in._ricochet_candidates.get_data(i));
    }
    
    _ricochet_candidate_velocities.reset();
    _ricochet_candidate_velocities.set_cols(_chisquared->get_dim());
    for(i=0;i<in._ricochet_candidate_velocities.get_rows();i++){
        for(j=0;j<in._ricochet_candidate_velocities.get_cols();j++){
            _ricochet_candidate_velocities.set(i,j,in._ricochet_candidate_velocities.get_data(i,j));
        }
    }
    
    _off_center_compass_points.reset();
    for(i=0;i<in._off_center_compass_points.get_dim();i++){
        _off_center_compass_points.set(i,in._off_center_compass_points.get_data(i));
    }
    
    _off_center_origins.reset();
    for(i=0;i<in._off_center_origins.get_dim();i++){
        _off_center_origins.set(i,in._off_center_origins.get_data(i));
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
    
    _distance_traveled.reset();
    for(i=0;i<in._distance_traveled.get_dim();i++){
        _distance_traveled.set(i,in._distance_traveled.get_data(i));
    }
    
    _ricochet_particles.reset();
    _ricochet_velocities.reset();
    _ricochet_velocities.set_cols(in._ricochet_velocities.get_cols());
    for(i=0;i<in._ricochet_velocities.get_rows();i++){
        for(j=0;j<in._ricochet_velocities.get_cols();i++){
            _ricochet_velocities.set(i,j,in._ricochet_velocities.get_data(i,j));
        }
    }
    for(i=0;i<in._ricochet_particles.get_dim();i++){
        _ricochet_particles.set(i,in._ricochet_particles.get_data(i));
    }
    
    _ricochet_strikes.reset();
    for(i=0;i<in._ricochet_strikes.get_dim();i++){
        _ricochet_strikes.set(i,in._ricochet_strikes.get_data(i));
    }
    
    _ricochet_discoveries.reset();
    for(i=0;i<in._ricochet_discoveries.get_rows();i++){
        for(j=0;j<in._ricochet_discoveries.get_cols(i);j++){
            _ricochet_discoveries.set(i,j,in._ricochet_discoveries.get_data(i,j));
        }
    }
    
    _ricochet_discovery_time.reset();
    for(i=0;i<in._ricochet_discovery_time.get_rows();i++){
        for(j=0;j<in._ricochet_discovery_time.get_cols(i);j++){
            _ricochet_discovery_time.set(i,j,in._ricochet_discovery_time.get_data(i,j));
        }
    }
    
    _ricochet_grad_norm.reset();
    for(i=0;i<in._ricochet_grad_norm.get_rows();i++){
        for(j=0;j<in._ricochet_dir_norm.get_cols(i);j++){
            _ricochet_grad_norm.set(i,j,in._ricochet_grad_norm.get_data(i,j));
        }
    }
    
    _ricochet_dir_norm.reset();
    for(i=0;i<in._ricochet_dir_norm.get_rows();i++){
        for(j=0;j<in._ricochet_dir_norm.get_cols(i);j++){
            _ricochet_dir_norm.set(i,j,in._ricochet_dir_norm.get_data(i,j));
        }
    }
    
    _ricochet_mu.reset();
    for(i=0;i<in._ricochet_mu.get_rows();i++){
        for(j=0;j<in._ricochet_mu.get_cols(i);j++){
            _ricochet_mu.set(i,j,in._ricochet_mu.get_data(i,j));
        }
    }
    
    _ricochet_strike_log.reset();
    for(i=0;i<in._ricochet_strike_log.get_rows();i++){
        for(j=0;j<in._ricochet_strike_log.get_cols(i);j++){
            _ricochet_strike_log.set(i,j,in._ricochet_strike_log.get_data(i,j));
        }
    }
    
    _ricochet_discovery_dexes.reset();
    for(i=0;i<in._ricochet_discovery_dexes.get_dim();i++){
        _ricochet_discovery_dexes.set(i,in._ricochet_discovery_dexes.get_data(i));
    }
    
    _ricochet_distances.reset();
    for(i=0;i<in._ricochet_distances.get_rows();i++){
        for(j=0;j<in._ricochet_distances.get_cols(i);j++){
            _ricochet_distances.set(i,j,in._ricochet_distances.get_data(i,j));
        }
    }
    
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

void node::set_center(int ix){
    _centerdex=ix;
    _min_changed=1;
    if(_chisquared!=NULL){
        _chimin=_chisquared->get_fn(ix);
    }
}

double node::distance_traveled(int ix){
    if(ix<0 || ix>=_chisquared->get_dim()){
        printf("WARNING asked for distance_traveled %d but %d\n",
        ix,_chisquared->get_dim());
    }
    
    if(ix>=_distance_traveled.get_dim()){
        return 0.0;
    }
    
    return _distance_traveled.get_data(ix);
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
}

int node::get_center(){
    return _centerdex;
}

double node::volume(){
    is_it_safe("volume");
    if(_min_found.get_dim()!=_chisquared->get_dim() || _max_found.get_dim()!=_chisquared->get_dim()){
        return 0.0;
    }
    
    double ans;
    int i;
    ans=1.0;
    for(i=0;i<_chisquared->get_dim();i++){
        ans*=(_max_found.get_data(i)-_min_found.get_data(i));
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
    
    if(_ricochet_particles.get_dim()!=_ricochet_velocities.get_rows()){
    
        printf("WARNING in node::%s\n",word);
        printf("ricochet particles %d\n",_ricochet_particles.get_dim());
        printf("ricochet velocities %d\n",_ricochet_velocities.get_rows());
        printf("ricochet strikes %d\n",_ricochet_strikes.get_dim());
        
        exit(1);
    
    }
    
    if(_ricochet_particles.get_dim()!=_ricochet_strikes.get_dim()){
        printf("WARNING in node::%s\n",word);
        printf("ricochet particles %d\n",_ricochet_particles.get_dim());
        printf("ricochet strikes %d\n",_ricochet_strikes.get_dim());
        exit(1);
    }
    
    if(_ricochet_candidates.get_dim()!=_ricochet_candidate_velocities.get_rows()){
        printf("WARNING in node::%s\n",word);
        printf("ricochet candidates %d\n",_ricochet_candidates.get_dim());
        printf("velocities %d\n",_ricochet_candidate_velocities.get_rows());
        exit(1);
    }
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
        trial.set(i,0.5*(_chisquared->get_pt(ibest,i)+_chisquared->get_pt(dex,i)));
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

void node::evaluate(array_1d<double> &pt, double *value, int *dex){
    is_it_safe("evaluate");
    
    _chisquared->evaluate(pt,value,dex);
    
    int i,j;
    array_1d<double> projected;
    projected.set_name("node_evaluate_projected");
    
    if(dex[0]>=0){
        if(value[0]<_chimin){
            _chimin=value[0];
            _centerdex=dex[0];
            _min_changed=1;
        }
        
        if(value[0]<=_chisquared->target()){
            j=1;
            for(i=0;i<_associates.get_dim() && j==1;i++){
                if(_associates.get_data(i)==dex[0])j=0;
            }
            
            if(j==1){
                _associates.add(dex[0]);
            }
            
            for(i=0;i<pt.get_dim();i++){
                if(i>=_min_found.get_dim() || pt.get_data(i)<_min_found.get_data(i)){
                    _min_found.set(i,pt.get_data(i));
                }
                
                if(i>=_max_found.get_dim() || pt.get_data(i)>=_max_found.get_data(i)){
                    _max_found.set(i,pt.get_data(i));
                }
            }
            
            project_to_bases(pt,projected);
            for(i=0;i<projected.get_dim();i++){
                if(i>=_projected_min.get_dim() || projected.get_data(i)<_projected_min.get_data(i)){
                    _projected_min.set(i,projected.get_data(i));
                }
                
                if(i>=_projected_max.get_dim() || projected.get_data(i)>_projected_max.get_data(i)){
                    _projected_max.set(i,projected.get_data(i));
                }
            }
            
            
        }
        
    }
}

double node::node_distance(array_1d<double> &p1, array_1d<double> &p2){
    int i;
    double norm,ans;
    ans=0.0;
    for(i=0;i<_chisquared->get_dim();i++){
        norm=_max_found.get_data(i)-_min_found.get_data(i);
        if(norm>0.0){
            ans+=power((p1.get_data(i)-p2.get_data(i))/norm,2);
        }
    }
    return sqrt(ans);

}

double node::node_distance(int i1, int i2){
    return node_distance(_chisquared->get_pt(i1)[0], _chisquared->get_pt(i2)[0]);
}

double node::node_distance(int i1, array_1d<double> &p2){
    return node_distance(p2, _chisquared->get_pt(i1)[0]);
}

double node::node_second_derivative_different(int center, int ix, int iy){
    is_it_safe("node_second_derivative_different");
    
    double xnorm,ynorm;
    xnorm=-1.0;
    ynorm=-1.0;
    if(_max_found.get_dim()>ix && _min_found.get_dim()>ix){
        xnorm=_max_found.get_data(ix)-_min_found.get_data(ix);
    }
    if(!(xnorm>0.0)){
        xnorm=_chisquared->get_max(ix)-_chisquared->get_min(ix);
    }
    
    if(_max_found.get_dim()>iy && _min_found.get_dim()>iy){
        ynorm=_max_found.get_data(iy)-_min_found.get_data(iy);
    }
    if(!(ynorm>0.0)){
        ynorm=_chisquared->get_max(iy)-_chisquared->get_min(iy);
    }
    
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
        trial.set(i,_chisquared->get_pt(center,i));
    }
    
    double xp,xm,yp,ym,fpp,fmp,fpm,fmm;
    double dx,xcenter,ycenter;
    
    int proceed,xpBound,xmBound,ypBound,ymBound;
    
    dx=1.0e-2;
    while(ifpp==ifpm || ifpp==ifmp || ifpp==ifmm ||
          ifpm==ifmp || ifpm==ifmm || ifmp==ifmm){
         
         
         proceed=0;
         
         while(proceed==0){
             xp=_chisquared->get_pt(center,ix)+dx*xnorm;
             xm=_chisquared->get_pt(center,ix)-dx*xnorm;
         
             yp=_chisquared->get_pt(center,iy)+dx*ynorm;
             ym=_chisquared->get_pt(center,iy)-dx*ynorm;
             
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
                 xcenter=_chisquared->get_pt(center,ix)-1.5*dx*xnorm;
                 proceed=0;
             }
             else if(xmBound==0){
                 xcenter=_chisquared->get_pt(center,ix)+1.5*dx*xnorm;
                 proceed=0;
             }
             else{
                 xcenter=_chisquared->get_pt(center,ix);
             }
             
             if(ypBound==0){
                 ycenter=_chisquared->get_pt(center,iy)-1.5*dx*ynorm;
                 proceed=0;
             }
             else if(ymBound==0){
                 ycenter=_chisquared->get_pt(center,iy)+1.5*dx*ynorm;
                 proceed=0;
             }
             else{
                 ycenter=_chisquared->get_pt(center,iy);
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
    
    xnorm=-1.0;
    if(_max_found.get_dim()>ix && _min_found.get_dim()>ix){
        xnorm=_max_found.get_data(ix)-_min_found.get_data(ix);
    }
    if(!(xnorm>0.0)){
        xnorm=_chisquared->get_max(ix)-_chisquared->get_min(ix);
    }
    
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
        trial.set(i,_chisquared->get_pt(center,i));
    }
    
    int proceed;
    int xppBound,xmmBound;
    double xcenter,ycenter,mu;
    
    while(ifpp==center || ifpp==ifmm || ifmm==center){
        proceed=0;
        while(proceed==0){
            xp=_chisquared->get_pt(center,ix)+dx*xnorm;
            xm=_chisquared->get_pt(center,ix)-dx*xnorm;
            xpp=_chisquared->get_pt(center,ix)+2.0*dx*xnorm;
            xmm=_chisquared->get_pt(center,ix)-2.0*dx*xnorm;
            
            xppBound = _chisquared->in_bounds(ix,xpp);
            xmmBound = _chisquared->in_bounds(ix,xmm);
            
            if(xppBound==0 && xmmBound==0){
                throw -1;
            }
            
            proceed=1;
            if(xppBound==0){
                xcenter=_chisquared->get_pt(center,ix)-2.5*dx*xnorm;
                proceed=0;
            }
            else if(xmmBound==0){
                xcenter=_chisquared->get_pt(center,ix)+2.5*dx*xnorm;
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
    is_it_safe("node_gradient");

    int i,j,if1,if2;
    array_1d<double> trial;
    double norm,x1,x2,y1,y2;
    trial.set_name("node_node_gradient_trial");
    
    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,_chisquared->get_pt(dex,i));
    }
    
    double dx,dxstart;
    dxstart=1.0e-2;
    
    for(i=0;i<_chisquared->get_dim();i++){
        norm=-1.0;
        if(_max_found.get_dim()>i && _min_found.get_dim()>i){
            norm=_max_found.get_data(i)-_min_found.get_data(i);
        }
        if(!(norm>0.0)){
           norm = _chisquared->get_max(i)-_chisquared->get_min(i);
        }
        
        dx=dxstart;
        if1=-1;
        if2=-1;
        while(if1==if2){
            x1=_chisquared->get_pt(dex,i)+dx*norm;
            trial.set(i,x1);
            evaluate(trial,&y1,&if1);
            
            x2=_chisquared->get_pt(dex,i)-dx*norm;
            trial.set(i,x2);
            evaluate(trial,&y2,&if2);
            
            if((if1<0 && if2>=0) or (if2<0 && if1>=0)){
                if(if1<0){
                    x1=_chisquared->get_pt(dex,i);
                    y1=_chisquared->get_fn(dex);
                    if1=dex;
                }
                else{
                    x2=_chisquared->get_pt(dex,i);
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
        
        trial.set(i,_chisquared->get_pt(dex,i));
    }

}

int node::node_bisection(array_1d<double> &ll, double fl,
                    array_1d<double> &hh, double fh,
                    int doSlope){

        return node_bisection(ll,fl,hh,fh,doSlope,_chisquared->target(),0.01*(_chisquared->target()-_chimin));
}

int node::node_bisection(array_1d<double> &lowball, double flow,
                    array_1d<double> &highball, double fhigh,
                    int doSlope, double target_value, double tolerance){

    is_it_safe("bisection");
    
    if(flow>fhigh){
        printf("WARNING in node bisection flow %e fhigh %e\n",flow,fhigh);
        exit(1);
    }
    
    double ftrial;
    array_1d<double> trial;
    trial.set_name("node_bisection_trial");
 
    int took_a_step=0,ct,i,iout;
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
        }
        else{
            fhigh=ftrial;
            for(i=0;i<_chisquared->get_dim();i++){
                highball.set(i,trial.get_data(i));
            }
        }
        
        
        ct++;
    }

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
                mu+=(_chisquared->get_pt(_basis_associates.get_data(ix),j)-_chisquared->get_pt(_centerdex,j))*trial_bases.get_data(i,j);
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

void node::compass_search(){
    compass_search(_centerdex);
}

void node::compass_search(int local_center){
        
    is_it_safe("compass_search");
    _compass_points.reset();
    _chisquared->set_iWhere(iCompass);
    
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
                flow=_chimin;
                for(i=0;i<_chisquared->get_dim();i++){
                    lowball.set(i,_chisquared->get_pt(local_center,i));
                }
            }
            else{
                for(i=0;i<_chisquared->get_dim();i++){
                    trial.set(i,_chisquared->get_pt(local_center,i)+sgn*dx*_basis_vectors.get_data(ix,i));
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
                flow=_chimin;
                for(i=0;i<_chisquared->get_dim();i++){
                    lowball.set(i,_chisquared->get_pt(local_center,i));
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
                flow=_chisquared->get_fn(local_center);
                for(i=0;i<_chisquared->get_dim();i++){
                    lowball.set(i,_chisquared->get_pt(local_center,i));
                }
            }
            iFound=node_bisection(lowball,flow,highball,fhigh,1);
            
            dx=0.0;
            for(i=0;i<_chisquared->get_dim();i++){
                dx+=_basis_vectors.get_data(ix,i)*(_chisquared->get_pt(local_center,i)-_chisquared->get_pt(iFound,i));
            }   
            
            dx=fabs(dx);
            
            if(iFound>=0){
                add_to_boundary(iFound);
                _compass_points.add(iFound);
                
                if(_chisquared->get_fn(iFound)>0.5*(_chimin+_chisquared->target()) && local_center==_centerdex){
                    for(i=0;i<_chisquared->get_dim();i++){
                        lowball.set(i,_chisquared->get_pt(local_center,i));
                        highball.set(i,_chisquared->get_pt(iFound,i));
                    }
                    
                    fhigh=_chisquared->get_fn(iFound);
                    bisection_target=0.5*(_chimin+_chisquared->target());
                    iFound=node_bisection(lowball,_chimin,highball,fhigh,1,bisection_target,0.1*(bisection_target-_chimin));

                    if(iFound>=0){
                        _basis_associates.add(iFound);

                        for(i=0;i<_chisquared->get_dim();i++){
                            lowball.set(i,_chisquared->get_pt(local_center,i));
                            highball.set(i,_chisquared->get_pt(iFound,i));
                        }
                        fhigh=_chisquared->get_fn(iFound);
                        bisection_target=0.5*(_chimin+fhigh);
                        iFound=node_bisection(lowball,_chimin,highball,fhigh,1,bisection_target,0.1*(bisection_target-_chimin));
                        if(iFound>=0){
                            _basis_associates.add(iFound);
                        }

                    }
                    
                }
                
            }
            
        }
        _basis_lengths.set(ix,blength);
    }
    
    printf("before off_diag %d\n",_chisquared->get_called()-ibefore);
    compass_diagonal(local_center);
    printf("leaving compass %d\n\n",_chisquared->get_called()-ibefore);
}

void node::compass_diagonal(int local_center){
    is_it_safe("compass_diagonal");
    _chisquared->set_iWhere(iCompass);

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
                dx+=0.5*(_chisquared->get_pt(j,k)-_chisquared->get_pt(i,k))*_basis_vectors.get_data(ix,k);
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
                    dy+=0.5*(_chisquared->get_pt(j,k)-_chisquared->get_pt(i,k))*_basis_vectors.get_data(iy,k);
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
                            trial.set(i,_chisquared->get_pt(local_center,i)+dmin*dir.get_data(i));
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
                            lowball.set(i,_chisquared->get_pt(local_center,i));
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
                            lowball.set(i,_chisquared->get_pt(local_center,i));
                        }
                        
                    }
                    iFound=node_bisection(lowball,flow,highball,fhigh,1);
                    
                    if(iFound>=0){
                        dmin=0.0;
                        mu=0.0;
                        for(i=0;i<_chisquared->get_dim();i++){
                            mu+=(_chisquared->get_pt(local_center,i)-_chisquared->get_pt(iFound,i))*_basis_vectors.get_data(ix,i);
                        }
                        dmin+=mu*mu;
                        mu=0.0;
                        for(i=0;i<_chisquared->get_dim();i++){
                            mu+=(_chisquared->get_pt(local_center,i)-_chisquared->get_pt(iFound,i))*_basis_vectors.get_data(iy,i);
                        }
                        dmin+=mu*mu;
                        dmin=sqrt(dmin);
                        if(xweight<0.0 && yweight<0.0)dmin_nn=dmin;
                        if(xweight<0.0 && yweight>0.0)dmin_np=dmin;
                        
                        add_to_boundary(iFound);
                        _compass_points.add(iFound);
                        if(_chisquared->get_fn(iFound)>0.5*(_chimin+_chisquared->target()) && local_center==_centerdex){
                            for(i=0;i<_chisquared->get_dim();i++){
                                lowball.set(i,_chisquared->get_pt(local_center,i));
                                highball.set(i,_chisquared->get_pt(iFound,i));
                            }
                            
                            fhigh=_chisquared->get_fn(iFound);
                            bisection_target=0.5*(_chimin+_chisquared->target());
                            iFound=node_bisection(lowball,_chimin,highball,fhigh,1,bisection_target,0.1*(bisection_target-_chimin));
                            if(iFound>=0){
                                _basis_associates.add(iFound);
                            }
                        
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

void node::compass_search_geometric_center(){   
    is_it_safe("compass_geometric_center");
    
    _chisquared->set_iWhere(iCompass); 
    array_1d<double> geometric_dir,trial,highball,lowball;
    int iFound,i;
    double mu,fhigh,flow,bisection_target,tol;
    
    geometric_dir.set_name("compass_search_geometric_dir");
    trial.set_name("compass_search_trial");
    
    iFound=-1;
    
    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,0.5*(_min_found.get_data(i)+_max_found.get_data(i)));
    }
    
    evaluate(trial,&mu,&iFound);
    
    if(iFound>=0){
        if(mu>_chisquared->target()){
            for(i=0;i<_chisquared->get_dim();i++){
                lowball.set(i,_chisquared->get_pt(_centerdex,i));
                highball.set(i,_chisquared->get_pt(iFound,i));
            }
            flow=_chisquared->get_fn(_centerdex);
            fhigh=_chisquared->get_fn(iFound);
            bisection_target=0.5*(_chisquared->get_fn(_centerdex)+_chisquared->target());
            tol=0.01*(_chisquared->target()-_chisquared->get_fn(_centerdex));
            
            iFound=node_bisection(lowball,flow,highball,fhigh,1,bisection_target,tol);
            
        }
        
        if(iFound>=0 && iFound!=_centerdex){
            compass_search(iFound);
        }
    
    }
    
    _volume_of_last_geom=volume();
    
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
        center.set(i,_chisquared->get_pt(iCenter,i));
        norm.set(i,-1.0);
        if(_max_found.get_dim()>i && _min_found.get_dim()>i){
            norm.set(i,_max_found.get_data(i)-_min_found.get_data(i));
        }
        if(!(norm.get_data(i)>0.0)){
            norm.set(i,_chisquared->get_max(i)-_chisquared->get_min(i));
        }
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
        norm=-1.0;
        if(_max_found.get_dim()>i && _min_found.get_dim()>i){
            norm=_max_found.get_data(i)-_min_found.get_data(i);
        }
        if(!(norm>0.0)){
            norm=_chisquared->get_max(i)-_chisquared->get_min(i);
        }
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
        center.set(i,_chisquared->get_pt(ans,i));
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

void node::recalibrate_projected_max_min(){
    is_it_safe("recalibrate_projected_max_min");

    _projected_max.reset();
    _projected_min.reset();
    _since_expansion=0;
    
    array_1d<double> projected;
    projected.set_name("node_recalibrate_projected");
    int i,j;
    
    
    for(i=0;i<_associates.get_dim();i++){
        project_to_bases(_chisquared->get_pt(_associates.get_data(i))[0],projected);
        for(j=0;j<_chisquared->get_dim();j++){
            if(j>=_projected_min.get_dim() || projected.get_data(j)<_projected_min.get_data(j)){
                _projected_min.set(j,projected.get_data(j));
            }
            
            if(j>=_projected_max.get_dim() || projected.get_data(j)>_projected_max.get_data(j)){
                _projected_max.set(j,projected.get_data(j));
            }
        }
    }
}

void node::find_bases(){
    is_it_safe("find_bases");
    
    if(_centerdex_basis<0 || node_distance(_centerdex_basis,_centerdex)>0.01){
        _basis_associates.reset();
    }
    
    if(_basis_associates.get_dim()==0){
        compass_search();
    }
    
    if(_basis_associates.get_dim()==0){
        printf("WARNING cannot find bases; no associates\n");
        return;
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
    
    target=0.5*(_chisquared->get_fn(_centerdex)+_chisquared->target());
    tol=0.01*(_chisquared->target()-_chisquared->get_fn(_centerdex));
    
    while(_basis_associates.get_dim()<_chisquared->get_dim()*_chisquared->get_dim()){
        for(i=0;i<_chisquared->get_dim();i++){
            dir.set(i,normal_deviate(_chisquared->get_dice(),0.0,1.0));
        }
        dir.normalize();
        
        flow=_chisquared->get_fn(_centerdex);
        for(i=0;i<_chisquared->get_dim();i++){
            lowball.set(i,_chisquared->get_pt(_centerdex,i));
            highball.set(i,_chisquared->get_pt(_centerdex,i));
        }
        fhigh=-2.0*exception_value;
        
        while(fhigh<target){
            for(i=0;i<_chisquared->get_dim();i++){
                highball.add_val(i,dir.get_data(i));
            }
            evaluate(highball,&fhigh,&iFound);
        }
        
        iFound=node_bisection(lowball,flow,highball,fhigh,1,target,tol);
        if(iFound>=0 && _chisquared->get_fn(iFound)>_chisquared->get_fn(_centerdex)+tol){
            _basis_associates.add(iFound);
        }
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
    
    printf("error0 %e %d\n",error0,_basis_associates.get_dim());

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
    
    
    while(stdev>stdevlim && aborted<max_abort){
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
    
    if(changed_bases==1){
        compass_search();
    }
   
    array_1d<double> projected;
    projected.set_name("node_find_bases_projected");
    
    if(changed_bases==1){
        recalibrate_projected_max_min();
    }
   
   
    _found_bases++;
    _min_changed=0;
    _centerdex_basis=_centerdex;
    printf("done finding bases\n");
    if(errorBest<_min_basis_error){
        _min_basis_error=errorBest;
        _min_basis_error_changed=1;
    }
}

void node::off_center_compass(int iStart){
    
    if(iStart<0){
        return;
    }
    
    int i,j,goAhead;
    double dd,ddmin;
    
    _chisquared->set_iWhere(iCompass);
    
    ddmin=1.0e-3;
    goAhead=1;
    for(i=0;i<_off_center_origins.get_dim() && goAhead==1;i++){
        dd=node_distance(iStart,_off_center_origins.get_data(i));
        if(dd<ddmin){
            goAhead=0;
        }
    }
    
    if(node_distance(iStart,_centerdex)<ddmin){
        goAhead=0;
    }
    
    if(goAhead==0){
        return;
    }
    
    _off_center_origins.add(iStart); 
    int ibefore=_chisquared->get_called();
    
    array_1d<double> lowball,highball,trial;
    double flow,fhigh,ftrial,dx;
    int iFound;
    
    lowball.set_name("node_off_center_compass_lowball");
    highball.set_name("node_off_center_compass_highball");
    trial.set_name("node_off_center_compass_trial");
    
    int ix;
    double sgn;
    for(ix=0;ix<_chisquared->get_dim();ix++){
        dx=1.0;
        for(sgn=-1.0;sgn<1.1;sgn+=2.0){
            flow=2.0*exception_value;
            fhigh=-2.0*exception_value;
            
            if(_chisquared->get_fn(iStart)>_chisquared->target()){
                for(i=0;i<_chisquared->get_dim();i++){
                    trial.set(i,0.5*(_chisquared->get_pt(_centerdex,i)+_chisquared->get_pt(iStart,i)));
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
                   trial.set(i,_chisquared->get_pt(iStart,i)+dx*sgn*_basis_vectors.get_data(ix,i));
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
                    lowball.set(i,_chisquared->get_pt(iStart,i));
                }
            }
            
            while(fhigh<_chisquared->target()){
                for(i=0;i<_chisquared->get_dim();i++){
                    highball.set(i,lowball.get_data(i)+dx*sgn*_basis_vectors.get_data(ix,i));
                }
                evaluate(highball,&fhigh,&iFound);
                dx+=1.0;
                dx*=2.0;
            }
            
            iFound=node_bisection(lowball,flow,highball,fhigh,1);
            
            if(iFound>=0){
                add_to_boundary(iFound);
                _off_center_compass_points.add(iFound);
                if(iFound!=iStart){
                    _ricochet_candidates.add(iFound);
                    j=_ricochet_candidate_velocities.get_rows();
                    for(i=0;i<_chisquared->get_dim();i++){
                        _ricochet_candidate_velocities.set(j,i,_chisquared->get_pt(iFound,i)-_chisquared->get_pt(iStart,i));
                    }

                    if(_ricochet_candidates.get_dim()!=_ricochet_candidate_velocities.get_rows()){
                        printf("WARNING out of sync in off_diag\n");
                        exit(1);
                    }
                }
            }
            
            if(sgn<0.0 && iFound>=0 && iFound!=iStart){
                dx=0.0;
                for(i=0;i<_chisquared->get_dim();i++){
                    dx+=(_chisquared->get_pt(iStart,i)-_chisquared->get_pt(iFound,i))*_basis_vectors.get_data(ix,i);
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
    
    printf("    done with off-center compass %d -- %d %d -- volume %e\n",
    _chisquared->get_called()-ibefore,iStart,_centerdex,volume());

}

double node::apply_quadratic_model(array_1d<double> &pt){
    is_it_safe("apply_quadratic_model");
    int ix,j;
    double ans,mu;
    ans=_chimin;
    for(ix=0;ix<_basis_vectors.get_rows();ix++){
        mu=0.0;
        for(j=0;j<_chisquared->get_dim();j++){
            mu+=(pt.get_data(j)-_chisquared->get_pt(_ellipse_center,j))*_basis_vectors.get_data(ix,j);
        }
        ans+=_basis_model.get_data(ix)*mu*mu;
    }
    return ans;
    
}

double node::ricochet_model(array_1d<double> &pt, kd_tree &tree){
    is_it_safe("ricochet_model");

    int npts=_chisquared->get_dim();
    double ell;
    array_2d<double> covar,covarin;
    array_1d<int> neigh;
    array_1d<double> dd;
    
    covar.set_name("node_ricochet_model_covar");
    covarin.set_name("node_ricochet_model_covarin");
    neigh.set_name("node_ricochet_model_neigh");
    dd.set_name("node_ricochet_model_dd");
    
    tree.nn_srch(pt,npts+1,neigh,dd);
    
    if(dd.get_data(0)<1.0e-20){
        neigh.remove(0);
        dd.remove(0);
    }
    
    array_1d<double> mutual_dd,mutual_dd_sorted;
    array_1d<int> mutual_dexes;
    
    mutual_dd.set_name("node_ricochet_model_mutual_dd");
    mutual_dd_sorted.set_name("node_ricochet_mutual_dd_sorted");
    mutual_dexes.set_name("node_ricochet_mutual_dexes");
    
    int i,j,k;
    
    k=0;
    for(i=0;i<npts;i++){
        for(j=i+1;j<npts;j++){
            mutual_dd.set(k,_chisquared->distance(neigh.get_data(i),neigh.get_data(j)));
            mutual_dexes.set(k,k);
            k++;
        }
    }
    
    sort_and_check(mutual_dd,mutual_dd_sorted,mutual_dexes);
    
    ell=mutual_dd_sorted.get_data(k/2);
    
    
    
    covar.set_dim(npts,npts);
    covarin.set_dim(npts,npts);
    
    double mu,nugget;
    nugget=1.0e-4;
    for(i=0;i<npts;i++){
        covar.set(i,i,1.0+nugget);
        for(j=i+1;j<npts;j++){
            mu=tree.distance(neigh.get_data(i),neigh.get_data(j));
            covar.set(i,j,exp(-0.5*power(mu/ell,2)));
            covar.set(j,i,covar.get_data(i,j));
        }
    }
    invert_lapack(covar,covarin,0);
    
    array_1d<double> qq;
    qq.set_name("node_ricochet_model_qq");
    for(i=0;i<npts;i++){
        mu=tree.distance(pt,neigh.get_data(i));
        qq.set(i,exp(-0.5*power(mu/ell,2)));
    }
    
    double fbar;
    fbar=apply_quadratic_model(pt);
    
    
    array_1d<double> qbar;
    qbar.set_name("node_ricochet_model_qbar");
    for(i=0;i<npts;i++){
        qbar.set(i,apply_quadratic_model(tree.get_pt(neigh.get_data(i))[0]));
    }
    
    mu=fbar;
    for(i=0;i<npts;i++){
        for(j=0;j<npts;j++){
            mu+=qq.get_data(i)*covarin.get_data(i,j)*(_chisquared->get_fn(neigh.get_data(j))-qbar.get_data(j));
        }
    }
    
    return mu;
}

double node::ricochet_distance(int i1, int i2){
    is_it_safe("ricochet_distance");
    
    if(_basis_lengths.get_dim()!=_chisquared->get_dim()){
        printf("WARNING in node::ricochet_distance basis_lengths.dim %d need %d\n",
        _basis_lengths.get_dim(),_chisquared->get_dim());
        
        exit(1);
    }
    
    double mu,ell,dd=0.0;
    int ix,i;
    for(ix=0;ix<_chisquared->get_dim();ix++){
        mu=0.0;
        for(i=0;i<_chisquared->get_dim();i++){
            mu+=(_chisquared->get_pt(i1,i)-_chisquared->get_pt(i2,i))*_basis_vectors.get_data(ix,i);
        }
        
        if(_basis_lengths.get_data(ix)>1.0e-10){
            dd+=power(mu/_basis_lengths.get_data(ix),2);
        }
        else{
            dd+=mu*mu;
        }
        
    }
    return sqrt(dd);
}

void node::add_to_boundary(int dex){
    _boundary_points.add(dex);
}

void node::_filter_candidates(){
    int i;

    for(i=0;i<_ricochet_candidates.get_dim();i++){
        if(_chisquared->get_fn(_ricochet_candidates.get_data(i))<0.5*(_chisquared->target()+_chisquared->chimin())){
            _ricochet_candidates.remove(i);
            _ricochet_candidate_velocities.remove_row(i);
            i--;
        }
    }
      
}

void node::initialize_ricochet(){
    is_it_safe("initialize_ricochet");
    
    if(_ricochet_candidates.get_dim()==0){
        find_bases();
    }
    
    int i;
    if(_ricochet_particles.get_dim()>0){
        for(i=0;i<_ricochet_particles.get_dim();i++){
	    _ricochet_candidates.add(_ricochet_particles.get_data(i));
            _ricochet_candidate_velocities.add_row(_ricochet_velocities(i)[0]);
	}
    }
    else{
        _ricochet_candidate_velocities.set_cols(_chisquared->get_dim());
    }
    
    _ricochet_velocities.reset();
    _ricochet_particles.reset();
    _ricochet_strikes.reset();
    _ricochet_velocities.set_cols(_chisquared->get_dim());
    _distance_traveled.reset();
    
    int nParticles=2*_chisquared->get_dim();
    
    array_1d<int> dexes;
    array_1d<double> dmu;
    
    int j,ix,iChosen;
    double dist,dist_best,dist_local_best;
    dexes.set_name("node_initialize_ricochet_dexes");
    dmu.set_name("node_initialize_ricochet_dmu");

    _filter_candidates();
    
    for(i=0;i<nParticles;i++){
        _ricochet_discovery_dexes.set(i,i);
    }
    
    for(i=0;i<_chisquared->get_dim();i++){
        _distance_traveled.set(i,0.0);
    }

    array_1d<double> dir;
    dir.set_name("node_initialize_ricochet_dir");
    
    dir.set_dim(_chisquared->get_dim());
    
    for(i=0;i<nParticles;i++){
        _ricochet_strikes.set(i,0);
        _ricochet_particles.set(i,0);
        _ricochet_velocities.add_row(dir);
        originate_particle(i,dir);
        for(j=0;j<_chisquared->get_dim();j++){
            _ricochet_velocities.set(i,j,dir.get_data(j));
        }
    }

    _min_basis_error_changed=0;
    _since_expansion=0;

    FILE *output;
    output=fopen("ricochet_particles.sav","w");
    for(ix=0;ix<_ricochet_particles.get_dim();ix++){
        for(i=0;i<_chisquared->get_dim();i++){
            fprintf(output,"%e ",_chisquared->get_pt(_ricochet_particles.get_data(ix),i));
        }
        fprintf(output,"\n");
    }
    fclose(output);
    output=fopen("compass_points.sav","w");
    for(ix=0;ix<_compass_points.get_dim();ix++){
        for(i=0;i<_chisquared->get_dim();i++){
            fprintf(output,"%e ",_chisquared->get_pt(_compass_points.get_data(ix),i));
        }
        fprintf(output,"\n");
    }
    fclose(output);

}

int node::smart_step_kick(int ix, double ratio, array_1d<double> &dir){
    array_1d<double> radial;
    radial.set_name("node_smart_step_kick_radial");
    int i;
    double mu;
    for(i=0;i<_chisquared->get_dim();i++){
        radial.set(i,_chisquared->get_pt(_centerdex,i)-_chisquared->get_pt(_ricochet_particles.get_data(ix),i));
    }
    mu=radial.normalize();
    
    if(mu<1.0e-20){
        printf("radial norm in smart_step_kick is small\n");
        return 0;
    }

    array_1d<double> trial;
    trial.set_name("node_smart_step_trial");
    int iFound;
    
    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,ratio*_chisquared->get_pt(_ricochet_particles.get_data(ix),i)+(1.0-ratio)*_chisquared->get_pt(_centerdex,i));
    }
    evaluate(trial,&mu,&iFound);
    
    if(iFound<0){
        printf("smart step kick found bogus particle\n");
        return 0;
    }
    
    _ricochet_particles.set(ix,iFound);
    
    array_2d<double> unitSphere;
    unitSphere.set_name("node_smart_step_kick_unitSphere");
    int j,k;
    for(i=0;i<_ricochet_candidates.get_dim();i++){
        for(j=0;j<_chisquared->get_dim();j++){
            trial.set(j,_chisquared->get_pt(_ricochet_candidates.get_data(i),j)-_chisquared->get_pt(_ricochet_particles.get_data(ix),j));
        }
        mu=trial.normalize();
        if(mu>1.0e-20){
            unitSphere.add_row(trial);
        }
    }
    
    for(i=0;i<_boundary_points.get_dim();i++){
        for(k=0;k<_chisquared->get_dim();k++){
            trial.set(k,_chisquared->get_pt(_boundary_points.get_data(i),k)-_chisquared->get_pt(_ricochet_particles.get_data(ix),k));
        }
        mu=trial.normalize();
        if(mu>1.0e-20){
            unitSphere.add_row(trial);
        }
    }
    
    int nAccepted,nCandidates=200;
    double dot,dotmax,dotbest;
    
    dotbest=2.0*exception_value;
    nAccepted=0;
    while(nAccepted<nCandidates){
        for(i=0;i<_chisquared->get_dim();i++){
            trial.set(i,normal_deviate(_chisquared->get_dice(),0.0,1.0));
        }
        trial.normalize();
        dot=0.0;
        for(i=0;i<_chisquared->get_dim();i++){
            dot+=trial.get_data(i)*radial.get_data(i);
        }
        if(dot<-1.0e-20){
            nAccepted++;
            
            dotmax=-2.0*exception_value;
            for(i=0;i<unitSphere.get_rows();i++){
                dot=0.0;
                for(j=0;j<_chisquared->get_dim();j++){
                    dot+=trial.get_data(j)*unitSphere.get_data(i,j);
                }
                
                if(dot>dotmax){
                    dotmax=dot;
                }
            }
            
            if(dotmax<dotbest){
                dotbest=dotmax;
                for(i=0;i<_chisquared->get_dim();i++){
                    dir.set(i,trial.get_data(i));
                }
            }
            
        }
    }
    
    return 1;
    

}

int node::t_kick(int ix, array_1d<double> &dir){

    if(_ricochet_origins.get_dim()<=ix || 
       _ricochet_origins.get_data(ix)<0 ||
       _ricochet_origins.get_data(ix)==_ricochet_particles.get_data(ix)){
       
        return 0;
    }
    
    double tol=1.0e-20;
    
    int i1,i0;
    
    i1=_ricochet_particles.get_data(ix);
    i0=_ricochet_origins.get_data(ix);
    
    int i,j,iMid;
    array_1d<double> midpt;
    midpt.set_name("node_t_kick_midpt");
    
    for(i=0;i<_chisquared->get_dim();i++){
        midpt.set(i,0.5*(_chisquared->get_pt(i1,i)+_chisquared->get_pt(i0,i)));
    }
    
    double mu;
    evaluate(midpt,&mu,&iMid);
    
    if(iMid<0){
        return 0;
    }
    
    double flow,fhigh;
    array_1d<double> lowball;
    int iFound;
    
    if(mu>_chisquared->target()){
        flow=_chisquared->get_fn(_centerdex);
        for(i=0;i<_chisquared->get_dim();i++){
            lowball.set(i,_chisquared->get_pt(_centerdex,i));
        }
        
        iFound=node_bisection(lowball, flow, midpt, mu, 0);
        if(iFound>=0){
            _ricochet_particles.set(ix,iFound);
            for(i=0;i<_chisquared->get_dim();i++){
                _ricochet_velocities.set(ix,i,_chisquared->get_pt(iFound,i)-_chisquared->get_pt(_centerdex,i));
            }
            _ricochet_velocities(ix)->normalize();
        }
        
        return 0;
    }
    
    int iOrigin;
    array_1d<double> origin;
    origin.set_name("node_t_kick_origin");
    for(i=0;i<_chisquared->get_dim();i++){
        origin.set(i,0.5*(_chisquared->get_pt(_centerdex,i)+_chisquared->get_pt(iMid,i)));
    }
    
    evaluate(origin,&mu,&iOrigin);
    
    if(iOrigin<0){
        return 0;
    }
    
    if(mu>_chisquared->target()){
        flow=_chisquared->get_fn(_centerdex);
        for(i=0;i<_chisquared->get_dim();i++){
            lowball.set(i,_chisquared->get_pt(_centerdex,i));
        }
        
        iFound=node_bisection(lowball, flow, origin, mu, 0);
        if(iFound>=0){
            _ricochet_particles.set(ix,iFound);
            for(i=0;i<_chisquared->get_dim();i++){
                _ricochet_velocities.set(ix,i,_chisquared->get_pt(iFound,i)-_chisquared->get_pt(_centerdex,i));
            }
            _ricochet_velocities(ix)->normalize();
        }
        
        return 0;
    }
    
    array_1d<double> trial;
    trial.set_name("node_t_kick_trial");
    double dnorm;
    for(i=0;i<_chisquared->get_dim();i++){
        trial.set(i,_chisquared->get_pt(i1,i)-_chisquared->get_pt(i0,i));
    }
    
    dnorm=trial.normalize();
    
    if(dnorm<tol){
        return 0;
    }
    
    _ricochet_particles.set(ix,iOrigin);
    for(i=0;i<_chisquared->get_dim();i++){
        dir.set(i,trial.get_data(i));
    }
    
    
    return 1;

}

int node::step_kick(int ix, double ratio, array_1d<double> &dir){

    int i,nearestParticle;
    double x1,x2,ddbest,ddmin,dd;
    nearestParticle=-1;
    ddmin=1.0e-10;
    
    for(i=0;i<_ricochet_candidates.get_dim();i++){
        dd=node_distance(_ricochet_candidates.get_data(i), _chisquared->get_pt(_ricochet_particles.get_data(ix))[0]);
        if(dd>ddmin){
            if(nearestParticle<0 || dd<ddbest){
                ddbest=dd;
                nearestParticle=_ricochet_candidates.get_data(i);
            }
        }
    }
    
    for(i=0;i<_boundary_points.get_dim();i++){
        if(_boundary_points.get_data(i)!=_ricochet_particles.get_data(ix)){
            dd=node_distance(_boundary_points.get_data(i),_ricochet_particles.get_data(ix));
            if(dd>ddmin){
                if(nearestParticle<0 || dd<ddbest){
                    ddbest=dd;
                    nearestParticle=_boundary_points.get_data(i);
                }
            }
        }
            
        
    }
    
    array_1d<double> trial;
    trial.set_name("node_step_kick_trial");
    int iFound;
    double mu;
    
    for(i=0;i<_chisquared->get_dim();i++){
           x1=_chisquared->get_pt(_ricochet_particles.get_data(ix),i);
           trial.set(i,ratio*x1+(1.0-ratio)*_chisquared->get_pt(_centerdex,i));
    }
    evaluate(trial,&mu,&iFound);
    if(iFound<0){
        return 0;
    }
    _ricochet_particles.set(ix,iFound);       
           
     if(nearestParticle>=0){
         for(i=0;i<_chisquared->get_dim();i++){
             dir.set(i,_chisquared->get_pt(_ricochet_particles.get_data(ix),i)-_chisquared->get_pt(nearestParticle,i));
         }
     }
     else{
         for(i=0;i<_chisquared->get_dim();i++){
             dir.set(i,_chisquared->get_pt(_ricochet_particles.get_data(ix),i)-_chisquared->get_pt(_centerdex,i));
         }
     }
     dir.normalize();
     return 1;
     //maybe should try reflecting about the gradient...
}


void node::originate_particle(int ix, array_1d<double> &dir){

    //choose new origin
    
    array_1d<double> lowball, highball,random_dir;
    double flow, fhigh,target;
    int iFound,i,j;
    
    while(_ricochet_candidates.get_dim()<2*_chisquared->get_dim()){
        lowball.set_name("node_originate_lowball");
        highball.set_name("node_originate_highball");
        random_dir.set_name("node_originate_random_dir");
        
        target=0.5*(_chisquared->get_fn(_centerdex)+_chisquared->target());
        
        for(i=0;i<_chisquared->get_dim();i++){
            random_dir.set(i,normal_deviate(_chisquared->get_dice(),0.0,1.0));
        }
        random_dir.normalize();
        
        flow=_chisquared->get_fn(_centerdex);
        for(i=0;i<_chisquared->get_dim();i++){
            lowball.set(i,_chisquared->get_pt(_centerdex,i));
            highball.set(i,_chisquared->get_pt(_centerdex,i));
        }
        
        fhigh=-2.0*exception_value;
        while(fhigh<=target){
            for(i=0;i<_chisquared->get_dim();i++){
                highball.add_val(i,random_dir.get_data(i));
            }
            evaluate(highball,&fhigh,&iFound);
        }
        
        iFound=node_bisection(lowball,flow,highball,fhigh,1,target,0.01);
        
        if(iFound>=0){
            off_center_compass(iFound);
        }
        
        _filter_candidates();
        
    }

    int iChosen=-1,iCandidate=-1;;
    double mu,dmu,dmubest;
    for(i=0;i<_ricochet_candidates.get_dim();i++){
        mu=apply_quadratic_model(_chisquared->get_pt(_ricochet_candidates.get_data(i))[0]);

        dmu=fabs(mu-_chisquared->get_fn(_ricochet_candidates.get_data(i)));
        if(iChosen<0 || dmu>dmubest){
            iChosen=_ricochet_candidates.get_data(i);
            iCandidate=i;
            dmubest=dmu;
        }
    }
    
    _ricochet_particles.set(ix,iChosen);
    _ricochet_origins.set(ix,-1);
   
    printf("iCandidate %d -- %d %d\n",iCandidate,_ricochet_candidates.get_dim(),_ricochet_candidate_velocities.get_rows());
    for(i=0;i<_chisquared->get_dim();i++){
        dir.set(i,_ricochet_candidate_velocities.get_data(iCandidate,i));
    }
    dir.normalize();
        
    _ricochet_candidate_velocities.remove_row(iCandidate);
    _ricochet_candidates.remove(iCandidate);
    printf("done with iCandidate\n");

    _ricochet_grad_norm.add(_ricochet_discovery_dexes.get_data(ix),-1.0);
    _ricochet_dir_norm.add(_ricochet_discovery_dexes.get_data(ix),-1.0);
    _ricochet_discoveries.add(_ricochet_discovery_dexes.get_data(ix),iChosen);
    _ricochet_distances.add(_ricochet_discovery_dexes.get_data(ix),-1.0);
    _ricochet_discovery_time.add(_ricochet_discovery_dexes.get_data(ix),_chisquared->get_called());
    _ricochet_mu.add(_ricochet_discovery_dexes.get_data(ix),-2.0*exception_value);
    _ricochet_strike_log.add(_ricochet_discovery_dexes.get_data(ix),-2);


    

}

int node::kick_particle(int ix, array_1d<double> &dir){
    return step_kick(ix,1.0-0.1*_ricochet_strikes.get_data(ix),dir);
}

void node::search(){

    if(_chisquared->get_fn(_centerdex)>_chisquared->target()){
        simplex_search();
        if(_chisquared->get_fn(_centerdex)>_chisquared->target()){
            _active=0;
            return;
        }
    }

    double minExpansionFactor=1.001;
 
    double projectedVolume0=projected_volume(); 
    double volume0=volume();
    int ibefore=_chisquared->get_called();
    
    ricochet();

    if(_ct_simplex<_ct_ricochet && 
        _failed_simplexes<3 && 
        _chisquared->could_it_go_lower(_chimin)>0){
        
        simplex_search();
    }
    
    
    double tol;
    tol=0.1*(_chisquared->target()-_chisquared->get_fn(_centerdex_basis));
    
    if(_chisquared->get_fn(_centerdex)<_chisquared->get_fn(_centerdex_basis)-tol){
        find_bases();
        compass_search_geometric_center();
    }
    
    double volume1=volume();
    double projectedVolume1=projected_volume();
    
    if(volume1>volume0*minExpansionFactor || projectedVolume1>projectedVolume0*minExpansionFactor){
        _since_expansion=0;
    }
    else{
        _since_expansion+=_ricochet_particles.get_dim();
    }
    
    if(_since_expansion>6*_chisquared->get_dim()){
        printf("deactivating because we did not expand\n");
        _active=0;
    }
    
    if(volume()>2.0*_volume_of_last_geom){
        compass_search_geometric_center();
    }
    
    if(_ricochet_particles.get_dim()==0){
        printf("deactivating because no particles are worth it\n");
        _active=0;
    }

    double maxDeltaV=0.1;

    if(_active==0){
        volume0=volume();
        projectedVolume0=projected_volume();
        
        find_bases();
        compass_search_geometric_center();
        
        volume1=volume();
        projectedVolume1=projected_volume();
        
        if(fabs(volume0-volume1)>maxDeltaV*volume0){
           initialize_ricochet();
           _active=1;
       }
    }
    
    if(_ricochet_particles.get_dim()==0){
        _active=0;
    }
    
    printf("failed kicks %d successful kicks %d\n",
    _failed_kicks,_successful_kicks);

}


void node::simplex_search(){
    is_it_safe("simplex_search");
    
    if(_min_found.get_dim()!=_chisquared->get_dim() || _max_found.get_dim()!=_chisquared->get_dim()){
        printf("    leaving node simplex search because there are no bounds\n");
        return;
    }
    
    printf("    node simplex search\n");
    
    int center0 = _centerdex;
    int ibefore=_chisquared->get_called();
    
    array_1d<int> dexes;
    array_1d<double> minpt,dmu,dmusorted;
    array_2d<double> seed;
    minpt.set_name("node_simplex_search_minpt");
    seed.set_name("node_simplex_search_minpt");
    dexes.set_name("node_simplex_search_dexes");
    dmu.set_name("node_simplex_search_dmu");
    dmusorted.set_name("node_simplex_search_dmusorted");
    
    seed.set_cols(_chisquared->get_dim());
    int i,j;
    double mu;

    _chisquared->set_iWhere(iNodeSimplex);

    if(_ricochet_particles.get_dim()<=_chisquared->get_dim()+1){
        for(i=0;i<_ricochet_particles.get_dim();i++){
            seed.add_row(_chisquared->get_pt(_ricochet_particles.get_data(i))[0]);
        }
        
        if(seed.get_rows()<_chisquared->get_dim()+1){
            for(i=0;i<_compass_points.get_dim();i++){
                dexes.set(i,_compass_points.get_data(i));
                mu=_chisquared->get_fn(_compass_points.get_data(i));
                dmu.set(i,fabs(mu-apply_quadratic_model(_chisquared->get_pt(_compass_points.get_data(i))[0])));
            }
            sort_and_check(dmu,dmusorted,dexes);
            for(i=dexes.get_dim()-1;i>=0 && seed.get_rows()<_chisquared->get_dim()+1;i--){
                seed.add_row(_chisquared->get_pt(dexes.get_data(i))[0]);
            }
        }
        
    }
    else{
        for(i=0;i<_ricochet_particles.get_dim();i++){
            dexes.set(i,i);
            mu=_chisquared->get_fn(_ricochet_particles.get_data(i));
            dmu.set(i,fabs(mu-apply_quadratic_model(_chisquared->get_pt(_ricochet_particles.get_data(i))[0])));
        }
        sort_and_check(dmu,dmusorted,dexes);
        
        for(i=dexes.get_dim()-1;seed.get_rows()<_chisquared->get_dim()+1;i--){
            seed.add_row(_chisquared->get_pt(_ricochet_particles.get_data(dexes.get_data(i)))[0]);
        }
        
    }
    
    if(seed.get_rows()!=_chisquared->get_dim()+1){
        printf("WARNING in node_simplex_search seed has %d rows want %d\n",
        seed.get_rows(),_chisquared->get_dim()+1);
        
        exit(1);
    }
    
    simplex_minimizer ffmin;
    ffmin.set_minmax(_min_found,_max_found);
    ffmin.set_chisquared(_chisquared);
    ffmin.set_dice(_chisquared->get_dice());
    ffmin.find_minimum(seed,minpt);

    evaluate(minpt,&mu,&i);

    _ct_simplex+=_chisquared->get_called()-ibefore;
    if(_centerdex==center0){
        _failed_simplexes++;
    }
    else{
        _failed_simplexes=0;
    }
    
}

void node::ricochet(){
    is_it_safe("ricochet");
    
    if(_ricochet_velocities.get_rows()==0){
        initialize_ricochet();
    }
    
    if(_ricochet_velocities.get_rows()!=_ricochet_particles.get_dim()){
        printf("WARNING in node n_velocities %d n_particles %d\n",
        _ricochet_velocities.get_rows(),_ricochet_particles.get_dim());
        
        exit(1);
    }
    
    if(_ricochet_velocities.get_cols()!=_chisquared->get_dim()){
   
       printf("WARNING in node ricochet velocity dim %d shld be %d\n",
       _ricochet_velocities.get_cols(),_chisquared->get_dim());
       
       exit(1);
       
   }
   
   _chisquared->set_iWhere(iRicochet);
   _calls_to_ricochet++;
   
   int ibefore=_chisquared->get_called();
   
   printf("    starting ricochet with volume %e and pts %d\n",volume(),
   _ricochet_particles.get_dim());
   
   int ix,i,j,iFound,updated;
   double flow,fhigh,eflow,efhigh;
   array_1d<double> lowball,highball,elowball,ehighball,edir,kick;
   
   lowball.set_name("node_ricochet_lowball");
   highball.set_name("node_ricochet_highball");
   elowball.set_name("node_ricochet_elowball");
   ehighball.set_name("node_ricochet_ehighball"); 
   edir.set_name("node_ricochet_edir");
   kick.set_name("node_ricochet_kick");
   
   kd_tree kd_copy(_chisquared->get_tree()[0]);
    
   double dx,x1,x2,y1,y2,component,distanceMin;
   double gnorm,dirnorm;
   array_1d<double> gradient,trial,dir,distanceMoved,chiFound;
   array_1d<int> end_pts,boundsChanged;
   
   array_2d<double> start_pts;
   array_1d<double> ricochet_max,ricochet_min,min0,max0;
   
   ricochet_max.set_name("node_ricochet_max");
   ricochet_min.set_name("node_ricochet_min");
   min0.set_name("node_ricochet_min0");
   max0.set_name("node_ricochet_max0");
   boundsChanged.set_name("node_ricochet_boundsChanges");
   distanceMoved.set_name("node_ricochet_distanceMoved");
   chiFound.set_name("node_ricochet_chiFound");
   
   for(i=0;i<_compass_points.get_dim();i++){
       for(j=0;j<_chisquared->get_dim();j++){
           if(i==0 || _chisquared->get_pt(_compass_points.get_data(i),j)<ricochet_min.get_data(j)){
               ricochet_min.set(j,_chisquared->get_pt(_compass_points.get_data(i),j));
           }
           if(i==0 || _chisquared->get_pt(_compass_points.get_data(i),j)>ricochet_max.get_data(j)){
               ricochet_max.set(j,_chisquared->get_pt(_compass_points.get_data(i),j));
           }
       }
   }
   
   for(i=0;i<_chisquared->get_dim();i++){
       min0.set(i,_min_found.get_data(i));
       max0.set(i,_max_found.get_data(i));
   }
   
   for(i=0;i<_ricochet_particles.get_dim();i++){
       boundsChanged.set(i,0);
   }
   
   start_pts.set_name("node_ricochet_start_pts");
   end_pts.set_name("node_ricochet_end_pts");
   gradient.set_name("node_ricochet_gradient");
   trial.set_name("node_ricochet_trial");
   dir.set_name("node_ricochet_dir");
   
   distanceMin=1.0e-2;
   for(ix=0;ix<_ricochet_particles.get_dim();ix++){
       start_pts.add_row(_chisquared->get_pt(_ricochet_particles.get_data(ix))[0]);
       flow=2.0*exception_value;
       fhigh=-2.0*exception_value;

       updated=0;
       if(_ricochet_strikes.get_data(ix)>0){
           updated=kick_particle(ix,dir);
       }

       if(updated==0){
           node_gradient(_ricochet_particles.get_data(ix),gradient);
       
           gnorm=gradient.normalize();
           component=0.0;
           for(i=0;i<_chisquared->get_dim();i++){
               component+=_ricochet_velocities.get_data(ix,i)*gradient.get_data(i);
           }
       
           for(i=0;i<_chisquared->get_dim();i++){
               dir.set(i,_ricochet_velocities.get_data(ix,i)-2.0*component*gradient.get_data(i));
           }
       
           dirnorm=dir.normalize();
       }

       flow=_chisquared->get_fn(_ricochet_particles.get_data(ix));
       for(i=0;i<_chisquared->get_dim();i++){
           lowball.set(i,_chisquared->get_pt(_ricochet_particles.get_data(ix),i));
       }

       while(flow>=_chisquared->target()){
           for(i=0;i<_chisquared->get_dim();i++){
               elowball.set(i,_chisquared->get_pt(_centerdex,i));
               ehighball.set(i,lowball.get_data(i));
               edir.set(i,lowball.get_data(i)-elowball.get_data(i));
           }
           edir.normalize();
           eflow=_chimin;
           
           component=1.0;
           efhigh=flow;

           if(eflow>_chisquared->target() || eflow>efhigh){
               printf("WARNING eflow %e %e %e\n",
               eflow,efhigh,_chisquared->target());
               printf("%e\n",_chisquared->get_fn(_centerdex));
               exit(1);
           }
           
           iFound=node_bisection(elowball,eflow,ehighball,efhigh,1);
           for(i=0;i<_chisquared->get_dim();i++){
               lowball.set(i,_chisquared->get_pt(iFound,i));
           }
           flow=_chisquared->get_fn(iFound);
           
       }
       
       if(flow>=_chisquared->target()){
           printf("WARNING in node ricochet flow %e; target %e\n",flow,_chisquared->target());
           exit(1);
       }
       
       for(i=0;i<_chisquared->get_dim();i++){
           highball.set(i,lowball.get_data(i));
       }
       
       component=1.0;
       while(fhigh<_chisquared->target()){
           for(i=0;i<_chisquared->get_dim();i++){
               highball.add_val(i,component*dir.get_data(i));
           }
           evaluate(highball,&fhigh,&j);
           component*=2.0;
           
           if(fhigh<_chisquared->target()){
               for(i=0;i<_chisquared->get_dim();i++){
                   lowball.set(i,highball.get_data(i));
               }
               flow=fhigh;
           }
           
       }
       
       if(flow>_chisquared->target() || flow>fhigh){
           printf("WARNING in ricochet %e %e %e\n",
           flow,fhigh,_chisquared->target());
           exit(1);
       }
       
       iFound=node_bisection(lowball,flow,highball,fhigh,0);
       if(iFound>=0){
           distanceMoved.set(ix,node_distance(start_pts(start_pts.get_rows()-1)[0],_chisquared->get_pt(iFound)[0]));
           chiFound.set(ix,_chisquared->get_fn(iFound));
           add_to_boundary(iFound);
 
       }
       else{
           distanceMoved.set(ix,0.0);
           chiFound.set(ix,_chisquared->chimin());
       }

       _ricochet_grad_norm.add(_ricochet_discovery_dexes.get_data(ix),gnorm);
       _ricochet_dir_norm.add(_ricochet_discovery_dexes.get_data(ix),dirnorm);
       _ricochet_discoveries.add(_ricochet_discovery_dexes.get_data(ix),iFound);
       _ricochet_distances.add(_ricochet_discovery_dexes.get_data(ix),distanceMoved.get_data(ix));
       _ricochet_discovery_time.add(_ricochet_discovery_dexes.get_data(ix),_chisquared->get_called());
       end_pts.add(iFound);

       _ricochet_origins.set(ix,_ricochet_particles.get_data(ix));
       _ricochet_particles.set(ix,iFound);
       for(i=0;i<_chisquared->get_dim();i++){
           _ricochet_velocities.set(ix,i,dir.get_data(i));
           
           if(_chisquared->get_pt(iFound,i)>max0.get_data(i) || _chisquared->get_pt(iFound,i)<min0.get_data(i)){
               boundsChanged.set(ix,1);
           }
           
           if(_max_found.get_data(i)-_min_found.get_data(i)>0.0){
               _distance_traveled.add_val(i,
                      fabs((start_pts.get_data(ix,i)-_chisquared->get_pt(iFound,i))/(_max_found.get_data(i)-_min_found.get_data(i))));
           }
           
       }
       
   }
   
   printf("    done with actual ricochet: volume %e\n",volume());
   
   double ricochet_dd,ddmax;
   int iChosen;

   if(end_pts.get_dim()!=_ricochet_particles.get_dim() || end_pts.get_dim()!=start_pts.get_rows()){
       printf("WARNING end_pts.dim %d number of ricochet particles %d\n",
       end_pts.get_dim(),_ricochet_particles.get_dim());
       printf("start pts %d\n",start_pts.get_rows());
       
       exit(1);
   }

   iChosen=-1;
   double mu;
   int isAStrike;
   for(i=0;i<_ricochet_particles.get_dim();i++){
       mu=-2.0*exception_value;
       if(boundsChanged.get_data(i)==0){
           mu=ricochet_model(_chisquared->get_pt(_ricochet_particles.get_data(i))[0],kd_copy);
       }
       _ricochet_mu.add(_ricochet_discovery_dexes.get_data(i),mu);

       isAStrike=0;

       if(boundsChanged.get_data(i)==0 && mu>0.0 &&
           mu<1.1*_chisquared->target()-0.1*_chisquared->chimin() &&
           mu>_chisquared->chimin()){
           
           isAStrike=1;
       }
       
       if(distanceMoved.get_data(i)<distanceMin ||
          chiFound.get_data(i)<0.1*_chisquared->chimin()+0.9*_chisquared->target()){
          
          isAStrike=1;   
       }
       
       if(isAStrike==1){
           if(_ricochet_strikes.get_data(i)>0){
               _failed_kicks++;
           }
           _ricochet_strikes.add_val(i,1);
           _ricochet_strike_log.add(_ricochet_discovery_dexes.get_data(i),_ricochet_strikes.get_data(i));
       }
       else{
           if(_ricochet_strikes.get_data(i)>0){
               _successful_kicks++;
           }
           _ricochet_strikes.set(i,0);
           _ricochet_strike_log.add(_ricochet_discovery_dexes.get_data(i),0);
       }
       
       if(_ricochet_strikes.get_data(i)>=_allowed_ricochet_strikes){
           originate_particle(i,_ricochet_velocities(i)[0]);
           _ricochet_strikes.set(i,0);
       }
       
       if(iChosen<0 || mu>ddmax){
           iChosen=i;
           ddmax=mu;
           for(j=0;j<_chisquared->get_dim();j++){
               trial.set(j,0.5*(start_pts.get_data(i,j)+_chisquared->get_pt(end_pts.get_data(i),j)));
           }
       }
       
   }
   
   _ct_ricochet+=_chisquared->get_called()-ibefore;
   int r_called=_chisquared->get_called()-ibefore;
   
   int totalNeedKick=0;
   for(i=0;i<_ricochet_strikes.get_dim();i++){
       if(_ricochet_strikes.get_data(i)>0)totalNeedKick++;
   }
   
   int do_off_center,k;
   trial.reset();
   for(i=0;i<_ricochet_particles.get_dim();i++){
       do_off_center=1;
       for(j=0;j<_chisquared->get_dim();j++){
           trial.set(j,0.5*(_chisquared->get_pt(_ricochet_particles.get_data(i),j)+_chisquared->get_pt(_centerdex,j)));
       }
       evaluate(trial,&mu,&iFound);
       if(mu<=_chisquared->target()){
           do_off_center=0;
       }
       for(j=0;j<_off_center_origins.get_dim() && do_off_center==1;j++){
           for(k=0;k<_chisquared->get_dim();k++){
               trial.set(k,0.5*(_chisquared->get_pt(_ricochet_particles.get_data(i),k)+_chisquared->get_pt(_off_center_origins.get_data(j),k)));
           }
           evaluate(trial,&mu,&iFound);
           if(mu<=_chisquared->target()){
               do_off_center=0;
           }
       }
       
       if(do_off_center==1){
           off_center_compass(_ricochet_particles.get_data(i));
       }
   }
   
   printf("    ending ricochet with volume %e -- %d -- %d -- need kick %d\n\n",
   volume(),r_called,_ricochet_particles.get_dim(),totalNeedKick);
}

int node::get_n_particles(){
    return _ricochet_particles.get_dim();
}

int node::get_n_candidates(){
    return _ricochet_candidates.get_dim();
}

void node::print_ricochet_discoveries(char *nameRoot){
    char outname[letters];
    int i,j,ix;
    FILE *output;
    
    for(ix=0;ix<_ricochet_discoveries.get_rows();ix++){
        sprintf(outname,"%s_%d.txt",nameRoot,ix);
        output=fopen(outname,"w");
        for(i=0;i<_ricochet_discoveries.get_cols(ix);i++){
            for(j=0;j<_chisquared->get_dim();j++){
                fprintf(output,"%le ",_chisquared->get_pt(_ricochet_discoveries.get_data(ix,i),j));
            }
            fprintf(output,"%le %le %d -- %e %d ",
                _chisquared->get_fn(_ricochet_discoveries.get_data(ix,i)),
                _ricochet_mu.get_data(ix,i),
                _ricochet_strike_log.get_data(ix,i),
                _ricochet_distances.get_data(ix,i),_ricochet_discovery_time.get_data(ix,i));
            fprintf(output,"-- grad %e dir %e ",
                _ricochet_grad_norm.get_data(ix,i),
                _ricochet_dir_norm.get_data(ix,i));
            fprintf(output,"\n");
        }
        fclose(output);
    }
    
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

void arrayOfNodes::add(int cc, chisq_wrapper *gg){

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
    
    _data[_ct].find_bases();

    _ct++;

}

void arrayOfNodes::add(chisq_wrapper *g, int i){
    add(i,g);
}

int arrayOfNodes::get_dim(){
    return _ct;
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
