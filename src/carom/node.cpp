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
    _bisection_tolerance=0.01;
    _ricochet_since_expansion=0;
    _min_changed=0;
    _active=1;
    _found_bases=0;
    _ct_ricochet=0;
    _calls_to_ricochet=0;
    
    _compass_points.set_name("node_compass_points");
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
    _ricochet_particles.set_name("node_ricochet_particles");
    _ricochet_velocities.set_name("node_ricochet_velocities");
    _ricochet_rr.set_name("node_ricochet_rr");
}

void node::copy(const node &in){
    _centerdex=in._centerdex;
    _chimin=in._chimin;
    _ricochet_since_expansion=in._ricochet_since_expansion;
    _min_changed=in._min_changed;
    _active=in._active;
    _found_bases=in._found_bases;
    _ct_ricochet=in._ct_ricochet;
    _calls_to_ricochet=in._calls_to_ricochet;
    
    
    int i,j;
    
    _chisquared=in._chisquared;
    
    _compass_points.reset();
    for(i=0;i<in._compass_points.get_dim();i++){
        _compass_points.set(i,in._compass_points.get_data(i));
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
    
    _ricochet_particles.reset();
    _ricochet_velocities.reset();
    _ricochet_rr.reset();
    _ricochet_particles.set_cols(in._ricochet_velocities.get_cols());
    _ricochet_velocities.set_cols(in._ricochet_velocities.get_cols());
    for(i=0;i<in._ricochet_velocities.get_rows();i++){
        _ricochet_rr.set(i,in._ricochet_rr.get_data(i));
        for(j=0;j<in._ricochet_velocities.get_cols();i++){
            _ricochet_particles.set(i,j,in._ricochet_particles.get_data(i,j));
            _ricochet_velocities.set(i,j,in._ricochet_velocities.get_data(i,j));
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
    if(_chisquared!=NULL){
        _chimin=_chisquared->get_fn(ix);
    }
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
    
    if(_ricochet_particles.get_rows()!=_ricochet_velocities.get_rows() ||
       _ricochet_rr.get_dim()!=_ricochet_particles.get_rows()){
    
        printf("WARNING in node::%s\n",word);
        printf("ricochet particles %d\n",_ricochet_particles.get_rows());
        printf("ricochet velocities %d\n",_ricochet_velocities.get_rows());
        printf("ricochet rr %d\n",_ricochet_rr.get_dim());
        
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
    
    /////////////////testing
    for(ix=0;ix<_chisquared->get_dim();ix++){
        mu=0.0;
        for(i=0;i<_chisquared->get_dim();i++){
            mu+=bases_out.get_data(ix,i)*bases_out.get_data(ix,i);
        }
        if(fabs(mu-1.0)>1.0e-6){
            printf("WARNING in node perturb bases, square norm %e\n",mu);
            exit(1);
        }
        
        for(jx=ix+1;jx<_chisquared->get_dim();jx++){
            mu=0.0;
            for(i=0;i<_chisquared->get_dim();i++){
                mu+=bases_out.get_data(ix,i)*bases_out.get_data(jx,i);
            }
            
            if(fabs(mu)>1.0e-6){
                printf("WARNING in node perturb bases, dot product %e\n",mu);
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
                if(sgn<0.0){
                    if(dx<0.0){
                        printf("WARNING dx is %e when should be positive for blength\n",dx);
                        exit(1);
                    }
                    blength=dx;
                }
                else{
                    if(dx>0.0){
                        printf("WARNING dx is %e when should be negative for blength\n",dx);
                        exit(1);
                    }
                    if(-1.0*dx<blength)blength=-1.0*dx;
                }
            }
            
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

void node::find_bases(){
    is_it_safe("find_bases");
    
    if(_basis_associates.get_dim()==0){
        compass_search();
    }
    
    _ellipse_center=_centerdex;
    _min_changed=0;
    
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
    while(ct<2000 && stdev>stdevlim && aborted<max_abort){
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
        _ellipse_center=_centerdex;
        compass_search();
    }
    
    _found_bases++;
    printf("done finding bases\n");
}

void node::off_center_compass(int iStart){
    
    int ibefore=_chisquared->get_called();
    
    array_1d<double> lowball,highball,trial;
    double flow,fhigh,ftrial,dx;
    int iFound;
    
    lowball.set_name("node_off_center_compass_lowball");
    highball.set_name("node_off_center_compass_highball");
    trial.set_name("node_off_center_compass_trial");
    
    int ix,i;
    double sgn;
    for(ix=0;ix<_chisquared->get_dim();ix++){
        dx=1.0;
        for(sgn=-1.0;sgn<1.1;sgn+=2.0){
            flow=2.0*exception_value;
            fhigh=-2.0*exception_value;
            
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
                dx*=2.0;
            }
            
            iFound=bisection(lowball,flow,highball,fhigh,1);
            
            if(sgn<0.0){
                dx=0.0;
                for(i=0;i<_chisquared->get_dim();i++){
                    dx+=(_chisquared->get_pt(iStart,i)-_chisquared->get_pt(iFound,i))*_basis_vectors.get_data(ix,i);
                }
                if(dx<0.0){
                    dx*=-1.0;
                }
            }
        }
    }
    
    printf("done with off-center compass %d -- %d %d\n",_chisquared->get_called()-ibefore,iStart,_centerdex);

}

double node::ricochet_model(array_1d<double> &pt, kd_tree &tree){
    is_it_safe("ricochet_model");

    int npts=5;
    double ell;
    array_2d<double> covar,covarin;
    array_1d<int> neigh;
    array_1d<double> dd;
    
    covar.set_name("node_ricochet_model_covar");
    covarin.set_name("node_ricochet_model_covarin");
    neigh.set_name("node_ricochet_model_neigh");
    dd.set_name("node_ricochet_model_dd");
    
    tree.nn_srch(pt,npts,neigh,dd);
    ell=dd.get_data(npts/2);
    
    double mu,nugget;
    int i,j,k;
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
    
    double fbar=0.0;
    for(i=0;i<npts;i++){
        fbar+=_chisquared->get_fn(neigh.get_data(i));
    }
    fbar=fbar/double(npts);
    
    mu=fbar;
    for(i=0;i<npts;i++){
        for(j=0;j<npts;j++){
            mu+=qq.get_data(i)*covarin.get_data(i,j)*(_chisquared->get_fn(neigh.get_data(j))-fbar);
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
    
    _ricochet_since_expansion=0;
    _ricochet_velocities.reset();
    _ricochet_particles.reset();
    _ricochet_velocities.set_cols(_chisquared->get_dim());
    _ricochet_particles.set_cols(_chisquared->get_dim());
    
    array_1d<int> dexes;
    array_1d<double> dd,ddsorted;
    int i,j;
    dexes.set_name("node_initialize_ricochet_dexes");
    dd.set_name("node_initialize_ricochet_dd");
    ddsorted.set_name("node_initialize_ricochet_ddsorted");
    
    for(i=0;i<_compass_points.get_dim();i++){
        dd.set(i,_chisquared->distance(_centerdex,_compass_points.get_data(i)));
        dexes.set(i,_compass_points.get_data(i));
    }
    
    sort_and_check(dd,ddsorted,dexes);
    
    int iUse;
    for(i=0;i<dexes.get_dim() && i<2*_chisquared->get_dim();i++){
        iUse=dexes.get_data(dexes.get_dim()-1-i);
        _ricochet_rr.set(i,ricochet_distance(_centerdex,iUse));
        for(j=0;j<_chisquared->get_dim();j++){
            _ricochet_particles.set(i,j,_chisquared->get_pt(iUse,j));
            _ricochet_velocities.set(i,j,_chisquared->get_pt(iUse,j)-_chisquared->get_pt(_centerdex,j));
        }
    }

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
   double volume0=volume();
   
   printf("    starting ricochet with volume %e\n",volume0);
   double flow,fhigh,eflow,efhigh;
   array_1d<double> lowball,highball,elowball,ehighball,edir;
   
   lowball.set_name("node_ricochet_lowball");
   highball.set_name("node_ricochet_highball");
   elowball.set_name("node_ricochet_elowball");
   ehighball.set_name("node_ricochet_ehighball"); 
   edir.set_name("node_ricochet_edir");
   
   kd_tree kd_copy(_chisquared->get_tree()[0]);
    
   int ix,i,j,iFound;
   double dx,x1,x2,y1,y2,component;
   array_1d<double> gradient,trial,dir;
   array_1d<int> end_pts;
   
   array_2d<double> start_pts;
   array_1d<double> ricochet_max,ricochet_min;
   
   ricochet_max.set_name("node_ricochet_max");
   ricochet_min.set_name("node_ricochet_min");
   
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
   
   start_pts.set_name("node_ricochet_start_pts");
   end_pts.set_name("node_ricochet_end_pts");
   gradient.set_name("node_ricochet_gradient");
   trial.set_name("node_ricochet_trial");
   dir.set_name("node_ricochet_dir");
   
   for(ix=0;ix<_ricochet_particles.get_rows();ix++){
       flow=2.0*exception_value;
       fhigh=-2.0*exception_value;
       try{
           _chisquared->find_gradient(_ricochet_particles(ix)[0],gradient);
       }
       catch(int iex){
           printf("ricochet failed to get gradient\n");
           exit(1);
           //code to do a brute force gradient if necessary
       }
       
       gradient.normalize();
       component=0.0;
       for(i=0;i<_chisquared->get_dim();i++){
           component+=_ricochet_velocities.get_data(ix,i)*gradient.get_data(i);
       }
       
       for(i=0;i<_chisquared->get_dim();i++){
           dir.set(i,_ricochet_velocities.get_data(ix,i)-2.0*component*gradient.get_data(i));
       }
       
       _chisquared->evaluate(_ricochet_particles(ix)[0],&flow,&i);
       for(i=0;i<_chisquared->get_dim();i++){
           lowball.set(i,_ricochet_particles.get_data(ix,i));
       }
       if(flow>=_chisquared->target()){
           for(i=0;i<_chisquared->get_dim();i++){
               elowball.set(i,_chisquared->get_pt(_centerdex,i));
               ehighball.set(i,elowball.get_data(i));
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
       
       start_pts.add_row(lowball);
       iFound=bisection(lowball,flow,highball,fhigh,0);
       end_pts.add(iFound);
       for(i=0;i<_chisquared->get_dim();i++){
           _ricochet_particles.set(ix,i,_chisquared->get_pt(iFound,i));
           _ricochet_velocities.set(ix,i,dir.get_data(i));
       }
       
   }
     
   double ricochet_dd,ddmax;
   int iChosen;

   if(end_pts.get_dim()!=_ricochet_particles.get_rows()){
       printf("WARNING end_pts.dim %d number of ricochet particles %d\n",
       end_pts.get_dim(),_ricochet_particles.get_rows());
       
       exit(1);
   }

   iChosen=-1;
   for(i=0;i<end_pts.get_dim();i++){
       ricochet_dd=ricochet_distance(_ellipse_center,end_pts.get_data(i));
       
       if(iChosen<0 || ricochet_dd>ddmax){
           ddmax=ricochet_dd;
           iChosen=end_pts.get_data(i);
           for(j=0;j<_chisquared->get_dim();j++){
               trial.set(j,0.5*(start_pts.get_data(i,j)+_chisquared->get_pt(end_pts.get_data(i),j)));
           }
       }
       _ricochet_rr.set(i,ricochet_dd);
       
       if(ricochet_dd<1.0+0.1*double(_calls_to_ricochet)){
           _ricochet_particles.remove_row(i);
           _ricochet_velocities.remove_row(i);
           _ricochet_rr.remove(i);
           end_pts.remove(i);
           start_pts.remove_row(i);
           i--;
       }
   }
   
   double ftrial;
   if(iChosen>=0){
       evaluate(trial,&ftrial,&iFound);
       if(iFound>=0){
           off_center_compass(iFound);
       }
   }
   
   if(_ricochet_particles.get_rows()==0){
       _active=0;
   }
   
   double volume1=volume();
   
   _ct_ricochet+=_chisquared->get_called()-ibefore;
   int r_called=_chisquared->get_called()-ibefore;
   
   if(volume1>1.001*volume0){
       _ricochet_since_expansion=0;
   }
   else{
       _ricochet_since_expansion++;
   }
   
   if(_ricochet_since_expansion>2){
       _active=0;
   }
   
   if(_active==0 && _found_bases<2){
       find_bases();
   }
   
   printf("    ending ricochet with volume %e -- %d -- %d\n\n",
   volume1,r_called,_ricochet_particles.get_rows());
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
