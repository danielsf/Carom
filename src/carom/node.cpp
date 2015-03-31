#include "node.h"

node::node(){
    initialize();
}

node::~node(){}

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
    _gradient_dex=-1;
    _bisection_tolerance=0.01;
    _min_changed=0;
    _active=1;
    _found_bases=0;
    _ct_ricochet=0;
    _ct_simplex=0;
    _calls_to_ricochet=0;
    _allowed_ricochet_strikes=4;
    _volume=0.0;
    _since_expansion=0;
    
    _compass_points.set_name("node_compass_points");
    _ricochet_candidates.set_name("node_ricochet_candidates");
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
    _distance_traveled.set_name("node_distance_traveled");
    _ricochet_particles.set_name("node_ricochet_particles");
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
}

void node::copy(const node &in){
    _centerdex=in._centerdex;
    _centerdex_basis=in._centerdex_basis;
    _gradient_dex=in._gradient_dex;
    _chimin=in._chimin;
    _chimin_bases=in._chimin_bases;
    _min_changed=in._min_changed;
    _active=in._active;
    _found_bases=in._found_bases;
    _ct_ricochet=in._ct_ricochet;
    _ct_simplex=in._ct_simplex;
    _calls_to_ricochet=in._calls_to_ricochet;
    _allowed_ricochet_strikes=in._allowed_ricochet_strikes;
    _ellipse_center=in._ellipse_center;
    _volume=in._volume;
    _since_expansion=in._since_expansion;
    
    int i,j;
    
    _chisquared=in._chisquared;
    
    _compass_points.reset();
    for(i=0;i<in._compass_points.get_dim();i++){
        _compass_points.set(i,in._compass_points.get_data(i));
    }
    
    _ricochet_candidates.reset();
    for(i=0;i<in._ricochet_candidates.get_dim();i++){
        _ricochet_candidates.set(i,in._ricochet_candidates.get_data(i));
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
    
    _distance_traveled.reset();
    for(i=0;i<in._distance_traveled.get_dim();i++){
        _distance_traveled.set(i,in._distance_traveled.get_data(i));
    }
    
    _ricochet_particles.reset();
    _ricochet_velocities.reset();
    _ricochet_particles.set_cols(in._ricochet_particles.get_cols());
    _ricochet_velocities.set_cols(in._ricochet_velocities.get_cols());
    for(i=0;i<in._ricochet_velocities.get_rows();i++){
        for(j=0;j<in._ricochet_velocities.get_cols();i++){
            _ricochet_velocities.set(i,j,in._ricochet_velocities.get_data(i,j));
        }
    }
    for(i=0;i<in._ricochet_particles.get_rows();i++){
        for(j=0;j<in._ricochet_particles.get_cols();j++){
            _ricochet_particles.set(i,j,in._ricochet_particles.get_data(i,j));
        }
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

int node::get_activity(){
    return _active;
}

int node::get_ct_ricochet(){
    return _ct_ricochet;
}

void node::set_center(int ix){
    _centerdex=ix;
    _gradient_dex=ix;
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
    
    if(_ricochet_particles.get_rows()!=_ricochet_velocities.get_rows()){
    
        printf("WARNING in node::%s\n",word);
        printf("ricochet particles %d\n",_ricochet_particles.get_rows());
        printf("ricochet velocities %d\n",_ricochet_velocities.get_rows());
        printf("ricochet strikes %d\n",_ricochet_strikes.get_dim());
        
        exit(1);
    
    }
    
    if(_ricochet_particles.get_rows()!=_ricochet_strikes.get_dim()){
        printf("WARNING in node::%s\n",word);
        printf("ricochet particles %d\n",_ricochet_particles.get_rows());
        printf("ricochet strikes %d\n",_ricochet_strikes.get_dim());
        exit(1);
    }
}

void node::evaluate(array_1d<double> &pt, double *value, int *dex){
    is_it_safe("evaluate");
    
    _chisquared->evaluate(pt,value,dex);
    
    int i;
    
    if(dex>=0){
        if(value[0]<_chimin){
            _chimin=value[0];
            _centerdex=dex[0];
            _gradient_dex=dex[0];
            _min_changed=1;
        }
        
        if(value[0]<=_chisquared->target()){
            for(i=0;i<pt.get_dim();i++){
                if(i>=_min_found.get_dim() || pt.get_data(i)<_min_found.get_data(i)){
                    _min_found.set(i,pt.get_data(i));
                }
                
                if(i>=_max_found.get_dim() || pt.get_data(i)>=_max_found.get_data(i)){
                    _max_found.set(i,pt.get_data(i));
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
    double dx,dxstart,xcenter,ycenter;
    
    int proceed,xpBound,xmBound,ypBound,ymBound;
    
    dxstart=1.0e-3;
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
    dx=1.0e-3;
    
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
    dxstart=1.0e-4;
    
    for(i=0;i<_chisquared->get_dim();i++){
        norm=_max_found.get_data(i)-_min_found.get_data(i);
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
            
            if(if1!=if2 || (if1<0 && if2<0)){
                grad.set(i,(y1-y2)/(x1-x2));
                if(if1<0 && if2<0){
                    if1=-1;
                    if2=-2;
                }
            }
            else{
                dx*=2.0;
            }
        }
        
        trial.set(i,_chisquared->get_pt(dex,i));
    }

}

int node::bisection(array_1d<double> &lowball, double flow, array_1d<double> &highball, double fhigh, int doSlope){

    is_it_safe("bisection");
    
    if(flow>fhigh){
        printf("WARNING in node bisection flow %e fhigh %e\n",flow,fhigh);
        exit(1);
    }
    
    double ftrial,threshold;
    array_1d<double> trial;
    trial.set_name("node_bisection_trial");
    
    if(_chisquared->get_deltachi()<0.0){
        threshold=_bisection_tolerance;
    }
    else{
        threshold=_bisection_tolerance*_chisquared->get_deltachi();
    }
    
    int took_a_step=0,ct,i,iout;
    double wgt;
    
    ct=0;
    iout=-1;
    while(ct<100 && (took_a_step==0 || _chisquared->target()-flow>threshold)){
        
        if(doSlope==1){
            wgt=(fhigh-_chisquared->target())/(fhigh-flow);
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
        if(ftrial<_chisquared->target()){
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

    return error;
    
}

void node::compass_search(){
    
    is_it_safe("compass_search");
    _compass_points.reset();
    _chisquared->set_iWhere(iCompass);
    
    int ibefore=_chisquared->get_called();
    
    int ix,i,j,iFound;
    double sgn,flow,fhigh,dx,ftrial,step;
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
                    lowball.set(i,_chisquared->get_pt(_centerdex,i));
                }
            }
            else{
                for(i=0;i<_chisquared->get_dim();i++){
                    trial.set(i,_chisquared->get_pt(_centerdex,i)+sgn*dx*_basis_vectors.get_data(ix,i));
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
                    lowball.set(i,_chisquared->get_pt(_centerdex,i));
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
                printf("WARNING in compass %e %e %e\n",
                flow,fhigh,_chisquared->target());
                exit(1);
            }
            iFound=bisection(lowball,flow,highball,fhigh,1);
            
            dx=0.0;
            for(i=0;i<_chisquared->get_dim();i++){
                dx+=_basis_vectors.get_data(ix,i)*(_chisquared->get_pt(_centerdex,i)-_chisquared->get_pt(iFound,i));
            }   
            
            dx=fabs(dx);
            
            if(iFound>=0){
                _compass_points.add(iFound);
                j=iFound;
                for(i=0;i<_chisquared->get_dim();i++){
                    trial.set(i,0.5*(_chisquared->get_pt(_centerdex,i)+_chisquared->get_pt(iFound,i)));
                }
                evaluate(trial,&ftrial,&iFound);
                if(iFound>=0){
                    _basis_associates.add(iFound);
                }
            
                for(i=0;i<_chisquared->get_dim();i++){
                    trial.set(i,0.75*_chisquared->get_pt(_centerdex,i)+0.25*_chisquared->get_pt(j,i));
                }
                evaluate(trial,&ftrial,&iFound);
                if(iFound>=0){
                    _basis_associates.add(iFound);
                }
            }
            
        }
        _basis_lengths.set(ix,blength);
    }
    
    printf("before off_diag %d\n",_chisquared->get_called()-ibefore);
    compass_off_diagonal();
    printf("leaving compass %d\n\n",_chisquared->get_called()-ibefore);
}

void node::compass_off_diagonal(){
    is_it_safe("compass_off_diagonal");
    _chisquared->set_iWhere(iCompass);

    int ix,iy;
    array_1d<double> trial,lowball,highball,dir;
    trial.set_name("node_off_diag_trial");
    lowball.set_name("node_off_diag_lowball");
    highball.set_name("node_off_diag_highball");
    dir.set_name("node_off_diag_dir");
    
    double dx,dy,dmin,step;
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
                            trial.set(i,_chisquared->get_pt(_centerdex,i)+dmin*dir.get_data(i));
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
                            lowball.set(i,_chisquared->get_pt(_centerdex,i));
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
                        printf("WARNING in off_diag %e %e %e\n",
                        flow,fhigh,_chisquared->target());
                        exit(1);
                    }
                    iFound=bisection(lowball,flow,highball,fhigh,1);
                    
                    if(iFound>=0){
                        dmin=0.0;
                        mu=0.0;
                        for(i=0;i<_chisquared->get_dim();i++){
                            mu+=(_chisquared->get_pt(_centerdex,i)-_chisquared->get_pt(iFound,i))*_basis_vectors.get_data(ix,i);
                        }
                        dmin+=mu*mu;
                        mu=0.0;
                        for(i=0;i<_chisquared->get_dim();i++){
                            mu+=(_chisquared->get_pt(_centerdex,i)-_chisquared->get_pt(iFound,i))*_basis_vectors.get_data(iy,i);
                        }
                        dmin+=mu*mu;
                        dmin=sqrt(dmin);
                        if(xweight<0.0 && yweight<0.0)dmin_nn=dmin;
                        if(xweight<0.0 && yweight>0.0)dmin_np=dmin;
                        
                        _compass_points.add(iFound);
                        for(i=0;i<_chisquared->get_dim();i++){
                            trial.set(i,0.5*(_chisquared->get_pt(_centerdex,i)+_chisquared->get_pt(iFound,i)));
                        }
                        evaluate(trial,&ftrial,&iFound);
                        if(iFound>=0){
                            _basis_associates.add(iFound);
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

void node::guess_bases(array_2d<double> &bases){
    is_it_safe("guess_bases");
    
    array_2d<double> covar;
    covar.set_name("node_guess_bases_covar");
    covar.set_cols(_chisquared->get_dim());
    bases.set_cols(_chisquared->get_dim());
    
    int ix,iy;
    double mu,covarmax=-1.0;
    for(ix=0;ix<_chisquared->get_dim();ix++){
        for(iy=ix;iy<_chisquared->get_dim();iy++){
            mu=node_second_derivative(_centerdex,ix,iy);
            if(fabs(mu)>covarmax){
                covarmax=fabs(mu);
            }
            covar.set(ix,iy,mu);
            if(ix!=iy){
                covar.set(iy,ix,mu);
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
    
    int i1=_chisquared->get_dim()/2;
    try{
        eval_symm(covar,evecs,evals,i1,_chisquared->get_dim(),1,0.001);
    }
    catch(int iex){
        printf("Guess failed on first batch of eigen vectors\n");
        throw -1;
    }
    
    for(ix=0;ix<i1;ix++){
        for(iy=0;iy<_chisquared->get_dim();iy++){
            bases.set(ix,iy,evecs.get_data(iy,ix));
        }
        bases(ix)->normalize();
    }
    
    printf("got first batch of guessed bases\n");
    
    try{
        eval_symm(covar,evecs,evals,_chisquared->get_dim()-i1,_chisquared->get_dim(),-1,0.001);
    }
    catch(int iex){
        printf("Guess failed on second batch of eigen vectors\n");
        throw -1;
    }
    
    for(ix=i1;ix<_chisquared->get_dim();ix++){
        for(iy=0;iy<_chisquared->get_dim();iy++){
            bases.set(ix,iy,evecs.get_data(iy,ix-i1));
        }
        bases(ix)->normalize();
    }
    
    printf("got second batch of guessed bases\n");
    
    validate_bases(bases,"node_guess_bases");
    
    printf("validated guessed bases\n");
}

void node::find_bases(){
    is_it_safe("find_bases");
    
    if(_basis_associates.get_dim()==0){
        compass_search();
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
   
    _found_bases++;
    _min_changed=0;
    _centerdex_basis=_centerdex;
    printf("done finding bases\n");
}

void node::off_center_compass(int iStart){
    
    if(iStart<0){
        return;
    }
    
    int i,goAhead;
    double dd,ddmin;
    
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
    _off_center_compass_points.reset();
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
            
            iFound=bisection(lowball,flow,highball,fhigh,1);
            
            if(iFound>=0){
                _off_center_compass_points.add(iFound);
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
    
    tree.nn_srch(pt,npts,neigh,dd);
    
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

void node::initialize_ricochet(){
    is_it_safe("initialize_ricochet");
    
    if(_compass_points.get_dim()==0){
        find_bases();
    }
    
    _ricochet_velocities.reset();
    _ricochet_particles.reset();
    _ricochet_strikes.reset();
    _ricochet_velocities.set_cols(_chisquared->get_dim());
    _ricochet_particles.set_cols(_chisquared->get_dim());
    _distance_traveled.reset();
    _ricochet_candidates.reset();
    
    array_1d<int> dexes,chosen_particles;
    array_1d<double> dd,ddsorted;
    array_2d<double> projected_candidates,projected_chosen;
    int i,j,ix,iChosen;
    double dist,dist_best,dist_local_best;
    dexes.set_name("node_initialize_ricochet_dexes");
    dd.set_name("node_initialize_ricochet_dd");
    ddsorted.set_name("node_initialize_ricochet_ddsorted");
    chosen_particles.set_name("node_initialize_ricochet_chosen_particles");
    projected_candidates.set_name("node_initialize_ricochet_projected_candidates");
    projected_chosen.set_name("node_initialize_ricochet_projected_chosen");

    for(i=0;i<_compass_points.get_dim();i++){
        if(_chisquared->get_fn(_compass_points.get_data(i))>0.5*(_chisquared->target()+_chisquared->chimin())){
            _ricochet_candidates.add(_compass_points.get_data(i));
        }
    }
    
    if(_ricochet_candidates.get_dim()<2*_chisquared->get_dim()){
        printf("\nBY THE WAY: JUST USING ALL OF THE CANDIDATES FOR RICOCHET\n\n");
        for(i=0;i<_ricochet_candidates.get_dim();i++){
            chosen_particles.add(_ricochet_candidates.get_data(i));
            _ricochet_candidates.remove(i);
            i--;
        }
    }
    else{
        iChosen=-1;
        dist_best=-2.0*exception_value;
        for(i=0;i<_ricochet_candidates.get_dim();i++){
            dist=fabs(apply_quadratic_model(_chisquared->get_pt(_ricochet_candidates.get_data(i))[0])-_chisquared->get_fn(_ricochet_candidates.get_data(i)));
            if(iChosen<0 || dist>dist_best){
                iChosen=i;
                dist_best=dist;
            }
        }
        if(iChosen<0){
            printf("WARNING could not choose a first ricochet particle\n");
            exit(1);
        }
        chosen_particles.add(_ricochet_candidates.get_data(iChosen));
        dd.reset();
        for(i=0;i<_chisquared->get_dim();i++){
            dd.set(i,_chisquared->get_pt(_ricochet_candidates.get_data(iChosen),i)-_chisquared->get_pt(_centerdex,i));
        }
        dist=dd.normalize();
        if(dist<1.0e-20){
            printf("WARNING first ricochet particle was a distance %e from center\n",dist);
            exit(1);
        }
        projected_chosen.add_row(dd);
        _ricochet_candidates.remove(iChosen);
        
        for(ix=0;ix<_ricochet_candidates.get_dim();ix++){
            for(i=0;i<_chisquared->get_dim();i++){
                dd.set(i,_chisquared->get_pt(_ricochet_candidates.get_data(ix),i)-_chisquared->get_pt(_centerdex,i));
            }
            dist=dd.normalize();
            if(dist<1.0e-20){
                _ricochet_candidates.remove(ix);
                ix--;
            }
            else{
                projected_candidates.add_row(dd);
            }
            
        }
        
        while(chosen_particles.get_dim()<2*_chisquared->get_dim() && _ricochet_candidates.get_dim()>0){
            if(_ricochet_candidates.get_dim()!=projected_candidates.get_rows()){
                printf("WARNING candidates %d projected %d\n",
                _ricochet_candidates.get_dim(), projected_candidates.get_rows());
                printf("chosen %d\n",chosen_particles.get_dim());
                exit(1);
            }
            
            iChosen=-1;
            for(ix=0;ix<_ricochet_candidates.get_dim();ix++){
                dist_local_best=-2.0*exception_value;
                for(i=0;i<chosen_particles.get_dim();i++){
                    dist=0.0;
                    for(j=0;j<_chisquared->get_dim();j++){
                        dist+=projected_candidates.get_data(ix,j)*projected_chosen.get_data(i,j);
                    }
                    if(dist>dist_local_best){
                        dist_local_best=dist;
                    }
                }
                if(iChosen<0 || dist_local_best<dist_best){
                    dist_best=dist_local_best;
                    iChosen=ix;
                }
            }
            
            chosen_particles.add(_ricochet_candidates.get_data(iChosen));
            projected_chosen.add_row(projected_candidates(iChosen)[0]);
            projected_candidates.remove_row(iChosen);
            _ricochet_candidates.remove(iChosen);
        
 
        }
    }

    _ricochet_discovery_dexes.reset();
    _ricochet_discoveries.reset();
    _ricochet_discovery_time.reset();
    _ricochet_distances.reset();
    _ricochet_grad_norm.reset();
    _ricochet_dir_norm.reset();
    _ricochet_mu.reset();
    _ricochet_strike_log.reset();
    for(i=0;i<chosen_particles.get_dim();i++){
        for(j=0;j<_chisquared->get_dim();j++){
            _ricochet_particles.set(i,j,_chisquared->get_pt(chosen_particles.get_data(i),j));
        }
        _ricochet_discovery_dexes.set(i,i);
        _ricochet_discoveries.set(i,0,chosen_particles.get_data(i));
        _ricochet_discovery_time.set(i,0,_chisquared->get_called());
        _ricochet_grad_norm.set(i,0,0.0);
        _ricochet_dir_norm.set(i,0,0.0);
        _ricochet_distances.set(i,0,0.0);
        _ricochet_mu.set(i,0,0.0);
        _ricochet_strike_log.set(i,0,0);
    }
    
    for(i=0;i<_chisquared->get_dim();i++){
        _distance_traveled.set(i,0.0);
    }

   int iOrigin;
   array_1d<double> grazing,radial;
   grazing.set_name("node_initialize_ricochet_grazing");
   radial.set_name("node_initialize_ricochet_radial");
   
   for(i=0;i<_ricochet_particles.get_rows();i++){
      iOrigin=-1;
      dist_best=2.0*exception_value;
      
      for(j=0;j<_chisquared->get_dim();j++){
          radial.set(j,_ricochet_particles.get_data(i,j)-_chisquared->get_pt(_centerdex,j));
      }
      radial.normalize();
      
      for(j=0;j<_ricochet_candidates.get_dim();j++){
           dist=_chisquared->distance(_ricochet_particles(i)[0],_ricochet_candidates.get_data(j));
           if(dist<dist_best){
              dist_best=dist;
              iOrigin=_ricochet_candidates.get_data(j);
           }
       }
       if(iOrigin<0){
           for(j=0;j<_chisquared->get_dim();j++){
               grazing.set(j,_ricochet_particles.get_data(i,j)-_chisquared->get_pt(_centerdex,j));
           }
       }
       else{
           for(j=0;j<_chisquared->get_dim();j++){
               grazing.set(j,radial.get_data(j));
           }
       }
       grazing.normalize();
       
       for(j=0;j<_chisquared->get_dim();j++){
           _ricochet_velocities.set(i,j,0.5*(grazing.get_data(j)+radial.get_data(j)));
       }
   }
   
 
    for(i=0;i<_ricochet_particles.get_rows();i++){
        _ricochet_strikes.set(i,0);
    }

    _volume=volume();

    FILE *output;
    output=fopen("ricochet_particles.sav","w");
    for(ix=0;ix<_ricochet_particles.get_rows();ix++){
        for(i=0;i<_chisquared->get_dim();i++){
            fprintf(output,"%e ",_ricochet_particles.get_data(ix,i));
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

void node::step_kick(int ix, double ratio, array_1d<double> &dir){

    int i,nearestParticle;
    double x1,x2,ddbest,ddmin,dd;
    nearestParticle=-1;
    ddmin=1.0e-10;
    
    for(i=0;i<_ricochet_candidates.get_dim();i++){
        dd=node_distance(_ricochet_candidates.get_data(i), _ricochet_particles(ix)[0]);
        if(dd>ddmin){
            if(nearestParticle<0 || dd<ddbest){
                ddbest=dd;
                nearestParticle=_ricochet_candidates.get_data(i);
            }
        }
    }
    
    int irow;
    for(irow=0;irow<_ricochet_discoveries.get_rows();irow++){
        if(irow!=_ricochet_discovery_dexes.get_data(ix)){
            for(i=0;i<_ricochet_discoveries.get_cols(irow);i++){
                dd=node_distance(_ricochet_discoveries.get_data(irow, i), _ricochet_particles(ix)[0]);
                if(dd>ddmin){
                    if(nearestParticle<0 || dd<ddbest){
                        ddbest=dd;
                        nearestParticle=_ricochet_discoveries.get_data(irow,i);
                    }
                }
            }
        }
    }

    for(i=0;i<_chisquared->get_dim();i++){
           x1=_ricochet_particles.get_data(ix,i);
           _ricochet_particles.set(ix,i,ratio*x1+(1.0-ratio)*_chisquared->get_pt(_centerdex,i));

     }
           
           
     if(nearestParticle>=0){
         for(i=0;i<_chisquared->get_dim();i++){
             dir.set(i,_ricochet_particles.get_data(ix,i)-_chisquared->get_pt(nearestParticle,i));
         }
     }
     else{
         for(i=0;i<_chisquared->get_dim();i++){
             dir.set(i,_chisquared->get_pt(_centerdex,i)-_ricochet_particles.get_data(ix,i));
         }
     }
     dir.normalize();
     
     //maybe should try reflecting about the gradient...
}


void node::origin_kick(int ix, array_1d<double> &dir){

    //choose new origin
    
    int iChosen=-1,iCandidate=-1;;
    int i,j;
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
    
    for(i=0;i<_chisquared->get_dim();i++){
        _ricochet_particles.set(ix,i,_chisquared->get_pt(iChosen,i));
    }
    
    _ricochet_candidates.remove(iCandidate);
    
    int irow,iOrigin;
    double dd,ddbest,ddmin;
    ddmin=1.0e-20;
    iOrigin=-1;
    for(i=0;i<_ricochet_candidates.get_dim();i++){
        dd=node_distance(_ricochet_candidates.get_data(i),_ricochet_particles(ix)[0]);
        if(dd>ddmin){
            if(iOrigin<0 || dd<ddbest){
                ddbest=dd;
                iOrigin=_ricochet_candidates.get_data(i);
            }
        }
    }
    
    for(irow=0;irow<_ricochet_discoveries.get_rows();irow++){
        for(i=0;i<_ricochet_discoveries.get_cols(irow);i++){
            dd=node_distance(_ricochet_discoveries.get_data(irow,i),_ricochet_particles(ix)[0]);
            if(dd>ddmin){
                if(iOrigin<0 || dd<ddbest){
                    ddbest=dd;
                    iOrigin=_ricochet_discoveries.get_data(irow,i);
                }
            }
        }
    }

    _ricochet_grad_norm.add(_ricochet_discovery_dexes.get_data(ix),-1.0);
    _ricochet_dir_norm.add(_ricochet_discovery_dexes.get_data(ix),-1.0);
    _ricochet_discoveries.add(_ricochet_discovery_dexes.get_data(ix),iChosen);
    _ricochet_distances.add(_ricochet_discovery_dexes.get_data(ix),-1.0);
    _ricochet_discovery_time.add(_ricochet_discovery_dexes.get_data(ix),_chisquared->get_called());
    _ricochet_mu.add(_ricochet_discovery_dexes.get_data(ix),-2.0*exception_value);
    _ricochet_strike_log.add(_ricochet_discovery_dexes.get_data(ix),-2);

    array_1d<double> gradient;
    gradient.set_name("node_origin_kick_gradient");
    _chisquared->find_gradient(_ricochet_particles(ix)[0],gradient);
    gradient.normalize();

    double component=0.0;
    for(i=0;i<_chisquared->get_dim();i++){
        dir.set(i,_ricochet_particles.get_data(ix,i)-_chisquared->get_pt(iOrigin,i));
        component+=dir.get_data(i)*gradient.get_data(i);
    }
    
    for(i=0;i<_chisquared->get_dim();i++){
        dir.subtract_val(i,2.0*component*gradient.get_data(i));
    }

}

void node::kick_particle(int ix, array_1d<double> &dir){
    /*if(_ricochet_strikes.get_data(ix)==1 || _ricochet_candidates.get_dim()==0){
        step_kick(ix,0.9,dir);
    }
    else{
        origin_kick(ix,dir);
    }*/
    step_kick(ix,(1.0-0.1*_ricochet_strikes.get_data(ix)),dir);
}

void node::search(){

    ricochet();

    if(_ct_simplex<_ct_ricochet){
        simplex_search();
    }
    
    int ibefore;
    if(_chimin_bases-_chimin>0.5*(_chisquared->target()-_chimin_bases)){
        ibefore=_chisquared->get_called();
        compass_search();
        find_bases();
        initialize_ricochet();
        _active=1;
        _ct_simplex+=_chisquared->get_called()-ibefore;
    }

    if(_active==0){
        find_bases();
    }

}


void node::simplex_search(){
    is_it_safe("simplex_search");
    
    if(_min_found.get_dim()!=_chisquared->get_dim() || _max_found.get_dim()!=_chisquared->get_dim()){
        printf("    leaving node simplex search because there are no bounds\n");
        return;
    }
    
    printf("    node simplex search\n");
    
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

    if(_ricochet_particles.get_rows()<=_chisquared->get_dim()+1){
        for(i=0;i<_ricochet_particles.get_rows();i++){
            seed.add_row(_ricochet_particles(i)[0]);
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
        for(i=0;i<_ricochet_particles.get_rows();i++){
            dexes.set(i,i);
            evaluate(_ricochet_particles(i)[0],&mu,&j);
            dmu.set(i,fabs(mu-apply_quadratic_model(_ricochet_particles(i)[0])));
        }
        sort_and_check(dmu,dmusorted,dexes);
        
        for(i=dexes.get_dim()-1;seed.get_rows()<_chisquared->get_dim()+1;i--){
            seed.add_row(_ricochet_particles(dexes.get_data(i))[0]);
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
    
}

void node::ricochet(){
    is_it_safe("ricochet");
    
    if(_ricochet_velocities.get_rows()==0){
        initialize_ricochet();
    }
    
    if(_ricochet_velocities.get_rows()!=_ricochet_particles.get_rows()){
        printf("WARNING in node n_velocities %d n_particles %d\n",
        _ricochet_velocities.get_rows(),_ricochet_particles.get_rows());
        
        exit(1);
    }
    
    if(_ricochet_velocities.get_cols()!=_chisquared->get_dim() ||
       _ricochet_particles.get_cols()!=_chisquared->get_dim()){
   
       printf("WARNING in node ricochet dim %d %d shld be %d\n",
       _ricochet_velocities.get_cols(),_ricochet_particles.get_cols(),
       _chisquared->get_dim());
       
       exit(1);
       
   }
   
   _chisquared->set_iWhere(iRicochet);
   _calls_to_ricochet++;
   
   int ibefore=_chisquared->get_called();
   
   printf("    starting ricochet with volume %e and pts %d\n",_volume,
   _ricochet_particles.get_rows());
   
   int ix,i,j,iFound;
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
   
   for(i=0;i<_ricochet_particles.get_rows();i++){
       boundsChanged.set(i,0);
   }
   
   start_pts.set_name("node_ricochet_start_pts");
   end_pts.set_name("node_ricochet_end_pts");
   gradient.set_name("node_ricochet_gradient");
   trial.set_name("node_ricochet_trial");
   dir.set_name("node_ricochet_dir");
   
   distanceMin=1.0e-2;
   for(ix=0;ix<_ricochet_particles.get_rows();ix++){
       start_pts.add_row(_ricochet_particles(ix)[0]);
       flow=2.0*exception_value;
       fhigh=-2.0*exception_value;
       if(_ricochet_strikes.get_data(ix)==0){
           try{
               _chisquared->find_gradient(_ricochet_particles(ix)[0],gradient);
           }
           catch(int iex){
               printf("ricochet failed to get gradient\n");
               exit(1);
               //code to do a brute force gradient if necessary
           }
       
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
       else{
           kick_particle(ix,dir);
       }

       evaluate(_ricochet_particles(ix)[0],&flow,&i);
       for(i=0;i<_chisquared->get_dim();i++){
           lowball.set(i,_ricochet_particles.get_data(ix,i));
       }

       if(flow>=_chisquared->target()){
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
               exit(1);
           }
           
           iFound=bisection(elowball,eflow,ehighball,efhigh,1);
           for(i=0;i<_chisquared->get_dim();i++){
               lowball.set(i,_chisquared->get_pt(iFound,i));
           }
           flow=_chisquared->get_fn(iFound);
           
       }
       
       if(flow>=_chisquared->target()){
           printf("WARNING in node ricochet flow %e\n",flow);
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
       
       iFound=bisection(lowball,flow,highball,fhigh,0);
       if(iFound>=0){
           distanceMoved.set(ix,node_distance(start_pts(start_pts.get_rows()-1)[0],_chisquared->get_pt(iFound)[0]));
           chiFound.set(ix,_chisquared->get_fn(iFound));
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
       for(i=0;i<_chisquared->get_dim();i++){
           _ricochet_particles.set(ix,i,_chisquared->get_pt(iFound,i));
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

   if(end_pts.get_dim()!=_ricochet_particles.get_rows() || end_pts.get_dim()!=start_pts.get_rows()){
       printf("WARNING end_pts.dim %d number of ricochet particles %d\n",
       end_pts.get_dim(),_ricochet_particles.get_rows());
       printf("start pts %d\n",start_pts.get_rows());
       
       exit(1);
   }

   array_1d<int> rejectThis;
   rejectThis.set_name("node_ricochet_rejectThis");
   for(i=0;i<_ricochet_particles.get_rows();i++){
       rejectThis.set(i,0);
   }

   iChosen=-1;
   double mu;
   int isAStrike;
   for(i=0;i<_ricochet_particles.get_rows();i++){
       mu=-2.0*exception_value;
       if(boundsChanged.get_data(i)==0){
           mu=ricochet_model(_ricochet_particles(i)[0],kd_copy);
       }
       _ricochet_mu.add(_ricochet_discovery_dexes.get_data(i),mu);

       isAStrike=0;

       if(boundsChanged.get_data(i)==0 &&
           mu<1.1*_chisquared->target()-0.1*_chisquared->chimin() &&
           mu>0.9*_chisquared->target()+0.1*_chisquared->chimin() &&
           mu>0.0 ){
           
           isAStrike=1;
       }
       
       if(distanceMoved.get_data(i)<distanceMin ||
          chiFound.get_data(i)<0.1*_chisquared->chimin()+0.9*_chisquared->target()){
          
          isAStrike=1;   
       }
       
       if(isAStrike==1){
           _ricochet_strikes.add_val(i,1);
           _ricochet_strike_log.add(_ricochet_discovery_dexes.get_data(i),_ricochet_strikes.get_data(i));
       }
       else{
           _ricochet_strikes.set(i,0);
           _ricochet_strike_log.add(_ricochet_discovery_dexes.get_data(i),0);
       }
       
       if(_ricochet_strikes.get_data(i)>=_allowed_ricochet_strikes){
           rejectThis.set(i,1); 
       }
       
       if(iChosen<0 || mu>ddmax){
           iChosen=i;
           ddmax=mu;
           for(j=0;j<_chisquared->get_dim();j++){
               trial.set(j,0.5*(start_pts.get_data(i,j)+_chisquared->get_pt(end_pts.get_data(i),j)));
           }
       }
       
   }
   
   
   double ftrial;
   array_1d<double> old_start;
   old_start.set_name("node_ricochet_old_start");
   int iMove=-1;
   
   /*
   if(iChosen>=0){
       ddmax=_chisquared->distance(_centerdex,end_pts.get_data(iChosen));
       evaluate(trial,&ftrial,&iFound);
       if(iFound>=0){
           off_center_compass(iFound);
           
           for(i=0;i<_off_center_compass_points.get_dim();i++){
               ricochet_dd=_chisquared->distance(_off_center_compass_points.get_data(i),_centerdex);
               if(ricochet_dd>ddmax){
                   ddmax=ricochet_dd;
                   iMove=_off_center_compass_points.get_data(i);
               }
           }
           
           if(iMove>=0){
               for(i=0;i<_chisquared->get_dim();i++){
                   old_start.set(i,_ricochet_particles.get_data(iChosen,i));
                   _ricochet_particles.set(iChosen,i,_chisquared->get_pt(iMove,i));
                   _ricochet_velocities.set(iChosen,i,_chisquared->get_pt(iMove,i)-_chisquared->get_pt(end_pts.get_data(iChosen),i));
                   
               }
               x1=node_distance(old_start, _chisquared->get_pt(iMove)[0]);
               
               _ricochet_grad_norm.add(_ricochet_discovery_dexes.get_data(iChosen),-1.0);
               _ricochet_dir_norm.add(_ricochet_discovery_dexes.get_data(iChosen),-1.0);
               _ricochet_discoveries.add(_ricochet_discovery_dexes.get_data(iChosen),iMove);
               _ricochet_distances.add(_ricochet_discovery_dexes.get_data(iChosen),x1);
               _ricochet_discovery_time.add(_ricochet_discovery_dexes.get_data(iChosen),_chisquared->get_called());
               _ricochet_mu.add(_ricochet_discovery_dexes.get_data(iChosen),-2.0*exception_value);
               _ricochet_strike_log.add(_ricochet_discovery_dexes.get_data(iChosen),-1);
               if(rejectThis.get_data(iChosen)==1){
                   for(i=0;i<_chisquared->get_dim();i++){
                       if(_ricochet_particles.get_data(iChosen,i)>max0.get_data(i) || _ricochet_particles.get_data(iChosen,i)<min0.get_data(i)){
                           rejectThis.set(iChosen,0);
                       }
                   }
               }
 
               if(rejectThis.get_data(iChosen)==1){
                   
                   mu=ricochet_model(_ricochet_particles(iChosen)[0],kd_copy);
                   if(mu<0.0 || mu>1.1*_chisquared->target()-0.1*_chisquared->chimin()){
                       rejectThis.set(iChosen,0);
                   }
               }
           }
       }
   }
   */
   
   for(i=0;i<_ricochet_particles.get_rows();i++){
       if(rejectThis.get_data(i)==1){
           _ricochet_particles.remove_row(i);
           _ricochet_velocities.remove_row(i);
           _ricochet_strikes.remove(i);
           _ricochet_discovery_dexes.remove(i);
           rejectThis.remove(i);
           i--;
       }
   }
   
   if(_ricochet_particles.get_rows()==0){
       printf("deactivating because no particles are worth it\n");
       _active=0;
   }
   
   double volume1=volume();

   if(volume1>_volume){
       _volume=volume1;
       _since_expansion=0;
   }
   else{
       _since_expansion+=_chisquared->get_called()-ibefore;
   }
 
   if(_since_expansion>1000){
       printf("deactivating because we have not expanded\n");
       _active=0;
   }
 
   _ct_ricochet+=_chisquared->get_called()-ibefore;
   int r_called=_chisquared->get_called()-ibefore;
   
   int totalNeedKick=0;
   for(i=0;i<_ricochet_strikes.get_dim();i++){
       if(_ricochet_strikes.get_data(i)>0)totalNeedKick++;
   }
   
   printf("    ending ricochet with volume %e -- %d -- %d -- need kick %d\n\n",
   volume1,r_called,_ricochet_particles.get_rows(),totalNeedKick);
}

int node::get_n_particles(){
    return _ricochet_particles.get_rows();
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
