#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "chisq.h"

void chisquared::allot_arrays(){
    int ix,iy;
    
    _bases.set_dim(_dim,_dim);
    _widths.set_dim(_ncenters,_dim);
    _centers.set_dim(_ncenters,_dim);     
    _nboundary.set_dim(_dim*_dim);
    _boundary_room.set_dim(_dim*_dim);

    
    for(ix=0;ix<_dim*_dim;ix++)_nboundary.set(ix,0);
    for(ix=0;ix<_dim*_dim;ix++)_boundary_room.set(ix,3);
    
    _boundary=new double**[_dim*_dim];
    for(ix=0;ix<_dim*_dim;ix++){
        _boundary[ix]=new double*[_boundary_room.get_data(ix)];
        for(iy=0;iy<_boundary_room.get_data(ix);iy++){
	    _boundary[ix][iy]=new double[3];
	}
        
    }
    
    _mins.set_dim(_dim);
    _maxs.set_dim(_dim);

    for(ix=0;ix<_dim;ix++){
        _mins.set(ix,2.0*exception_value);
	_maxs.set(ix,-2.0*exception_value);
    } 
    
    _bases.set_name("chisq_bases");
    _widths.set_name("chisq_widths");
    _centers.set_name("chisq_centers");
    _nboundary.set_name("chisq_nboundary");
    _boundary_room.set_name("chisq_boundary_room");
    _mins.set_name("chisq_mins");
    _maxs.set_name("chisq_maxs");
    
}

void chisquared::set_max(int dex, double nn){
    _maxs.set(dex,nn);
}

void chisquared::set_min(int dex, double nn){
    _mins.set(dex,nn);
}

double chisquared::get_min(int dex){
    return _mins.get_data(dex);
}

double chisquared::get_max(int dex){
    return _maxs.get_data(dex);
}

double chisquared::get_time_spent(){
    return _time_spent;
}

void chisquared::reset_boundary(){
    int i;
    for(i=0;i<_dim*_dim;i++)_nboundary.set(i,0);
}

void chisquared::make_bases(int seed){
    make_bases(seed, 1);
}

void chisquared::make_bases(int seed, int doCenters){

    int do_random_bases=1;
    
    if(seed==0){
        seed=int(time(NULL));
    }
    else if(seed<0){
        seed*=-1;
        do_random_bases=0;
    }
    
    if(_dice==NULL){
        _dice=new Ran(seed);
    }
    
    double nn;
    int i,j,ii,jj,goon;
    
    _centers.set_where("chisq_make_bases");
    _bases.set_where("chisq_make_bases");
    _widths.set_where("chisq_make_bases");
    
    if(do_random_bases==1){
        for(ii=0;ii<_dim;ii++){
            goon=1;
	    while(goon==1){
                goon=0;
	        for(i=0;i<_dim;i++)_bases.set(ii,i,_dice->doub()-0.5);
	        for(jj=0;jj<ii;jj++){
	            nn=0.0;
		    for(i=0;i<_dim;i++)nn+=_bases.get_data(ii,i)*_bases.get_data(jj,i);
		    for(i=0;i<_dim;i++)_bases.subtract_val(ii,i,nn*_bases.get_data(jj,i));
	        }
	    
	        nn=0.0;
	        for(i=0;i<_dim;i++){
		    nn+=power(_bases.get_data(ii,i),2);
	        }
	        if(nn<1.0e-20)goon=1;
	        nn=sqrt(nn);
	        for(i=0;i<_dim;i++){
	            _bases.divide_val(ii,i,nn);
	        }
	    }
        }
    }//if do_random_bases==1
    else{
        for(i=0;i<_dim;i++){
            for(jj=0;jj<_dim;jj++){
                if(i==jj)_bases.set(i,jj,1.0);
                else _bases.set(i,jj,0.0);
            }
        }
    }
    
    double normerr,ortherr;
    for(ii=0;ii<_dim;ii++){
        nn=0.0;
	for(i=0;i<_dim;i++)nn+=power(_bases.get_data(ii,i),2);
	nn=fabs(1.0-nn);
	if(ii==0 || nn>normerr)normerr=nn;
	
	for(jj=ii+1;jj<_dim;jj++){
	   nn=0.0;
	   for(i=0;i<_dim;i++)nn+=_bases.get_data(ii,i)*_bases.get_data(jj,i);
	   nn=fabs(nn);
	   if((ii==0 && jj==1) || nn>ortherr)ortherr=nn;
	}
    }
    
    printf("normerr %e ortherr %e\n",normerr,ortherr);
    if(normerr>1.0e-3 || ortherr>1.0e-3){
        death_knell("normerr or ortherr too large");
    }
    
    
    if(doCenters==1){
        make_centersRandom();
    }
    
    _centers.set_where("nowhere");
    _bases.set_where("nowhere");
    _widths.set_where("nowhere");
    
    _time_spent=0.0;
    _called=0;
    
    //printf("set centers and widths %d %d\n",_dim,_ncenters);
    
}  

void chisquared::make_centersRandom(){ 
    
    /*make centers for non_random bases*/
    
    double nn;
    int i,j,ii,jj,goon;
    
    double rr,theta,dx,dy;
    array_1d<double> trial_center,trial_pt;
    int acceptable,iterations=0;;
    
    trial_center.set_dim(_dim);
    trial_pt.set_dim(_dim);

    goon=0;
    while(goon==0){
        iterations++;
        
        if(iterations>5000){
            printf("WARNING; chisq was unable to construct an acceptable function\n");
            printf("ncenters %d dim %d\n",_ncenters,_dim);
            
            exit(1);
        }
        
	for(ii=0;ii<_ncenters;ii++){
            for(i=0;i<_dim;i++)_centers.set(ii,i,1.0e30);
	    for(i=0;i<_dim;i++)_widths.set(ii,i,1.0e-5);
        }
	
	goon=1;
        for(ii=0;ii<_ncenters;ii++){
	    
	    acceptable=1;
	    
	    for(i=0;i<_dim;i++)trial_center.set(i,0.0);
	    
            for(i=0;i<_dim;i++){
                trial_center.add_val(i,normal_deviate(_dice,0.0,30.0));
	    }
            
	    if(ii>0){
	        rr=normal_deviate(_dice,40.0,20.0);
	        theta=_dice->doub()*2.0*pi;
	        
                if(cos(theta)<0.0)dx=-2.0;
                else dx=2.0;
                
                if(sin(theta)<0.0)dy=-2.0;
                else dy=2.0;
                
		trial_center.set(0,_centers.get_data(0,0)+(rr*cos(theta)+dx)*_widths.get_data(0,0));
		trial_center.set(1,_centers.get_data(0,1)+(rr*sin(theta)+dy)*_widths.get_data(0,1));
	 
	    } 

	    for(i=0;i<_dim;i++){
		trial_pt.set(i,0.0);
	        for(j=0;j<_dim;j++)trial_pt.add_val(i,trial_center.get_data(j)*_bases.get_data(j,i));
	    }
	    
	    rr=(*this)(trial_pt);
	    
	    if(rr<100.0)acceptable=0;
	    
	    if(acceptable==1){
	        for(i=0;i<_dim;i++){
		    _centers.set(ii,i,trial_center.get_data(i));
		    _widths.set(ii,i,fabs(normal_deviate(_dice,3.0-(1.8/21.0)*double(_dim),0.05))+0.05);
	        }
	    }
	    else ii--;
	
        }
	
	
	
	
	
	for(i=0;i<_dim && acceptable==1;i++){
	    for(ii=0;ii<_ncenters && acceptable==1;ii++){
                if(_centers.get_data(ii,i)-4.0*_widths.get_data(ii,i)<-100.0)acceptable=0;
                if(_centers.get_data(ii,i)+4.0*_widths.get_data(ii,i)>100.0)acceptable=0;
                
	        for(jj=ii+1;jj<_ncenters && acceptable==1;jj++){
		    nn=fabs(_centers.get_data(ii,i)-_centers.get_data(jj,i));
		    if(nn<2.0*_widths.get_data(ii,i) || nn<2.0*_widths.get_data(jj,i)){
                        acceptable=0;
                        //printf("centers %d %d -- %d -- %e %e -- %e %e\n",ii,jj,i,centers.get_data(ii,i),_widths.get_data(ii,i),
                        //centers.get_data(jj,i),_widths.get_data(jj,i));
                    }
		}
	    }
	}
        if(acceptable==0)goon=0;
	
    }
    

}

void chisquared::add_to_boundary(array_1d<double> &alpha, int ix, int iy,double chitest){
    
    array_2d<double> buffer;
    array_1d<double> pt;
    double rr;
    int ipt,room;
    int i,j;
    
    buffer.set_name("chisq_add_to_boundary_buffer");
    pt.set_name("chisq_add_to_boundary_pt");
    
    //printf("adding %e\n",chitest);
    if(ix>iy){
        i=ix;
	ix=iy;
	iy=i;
    }
    
    ipt=ix*_dim+iy;
    room=_boundary_room.get_data(ipt);
    
    if(_nboundary.get_data(ipt)>=room){
  
	buffer.set_dim(room,3);
	for(i=0;i<room;i++){
	    for(j=0;j<3;j++)buffer.set(i,j,_boundary[ipt][i][j]);
	    delete [] _boundary[ipt][i];
	}
        delete [] _boundary[ipt];
	_boundary_room.add_val(ipt,100);
	_boundary[ipt]=new double*[_boundary_room.get_data(ipt)];
	for(i=0;i<_boundary_room.get_data(ipt);i++)_boundary[ipt][i]=new double[3];
	
	for(i=0;i<room;i++){
	    for(j=0;j<3;j++)_boundary[ipt][i][j]=buffer.get_data(i,j);
	}
	buffer.reset();
    }
    
    _boundary[ipt][_nboundary.get_data(ipt)][0]=alpha.get_data(ix);
    _boundary[ipt][_nboundary.get_data(ipt)][1]=alpha.get_data(iy);
    _boundary[ipt][_nboundary.get_data(ipt)][2]=chitest;
    
    _nboundary.add_val(ipt,1);
    
    rr=0.0;
    for(i=0;i<_dim;i++)rr+=alpha.get_data(i)*alpha.get_data(i);
    rr=sqrt(rr);
    if(rr>_rr_max)_rr_max=rr;
    
    pt.set_dim(_dim);
    for(i=0;i<_dim;i++)pt.set(i,0.0);
    for(i=0;i<_dim;i++){
        for(j=0;j<_dim;j++){
	    pt.add_val(i,alpha.get_data(j)*_bases.get_data(j,i));
	}
    }
    
    double nn=0.0;
    for(i=0;i<_dim;i++)nn+=pt.get_data(i)*pt.get_data(i);
    nn=sqrt(nn);
    if(fabs(nn-rr)>1.0e-4){
        printf("WARNING nn %e rr %e\n",nn,rr);
    }
    
    //not sure why this is here
    //I think I was using it for analysis of cartoons
    /*for(i=0;i<_dim;i++){
        if(pt.get_data(i)<_mins.get_data(i))_mins.set(i,pt.get_data(i));
	if(pt.get_data(i)>_maxs.get_data(i))_maxs.set(i,pt.get_data(i));
    }*/
    
}

void chisquared::print_mins_maxs(){
    int i;
    double nn;
    nn=0.0;
    printf("mins and maxs\n");
    for(i=0;i<_dim;i++){
        printf("p%d -- %e %e -- %e\n",i,_mins.get_data(i),_maxs.get_data(i),_maxs.get_data(i)-_mins.get_data(i));
	nn+=power(_maxs.get_data(i)-_mins.get_data(i),2);
    }
    printf("\nfiducial distance %e\n",sqrt(nn));
}

double chisquared::get_rr_max(){
    return _rr_max;
}

int chisquared::get_n_boundary(int ix, int iy){
    
    int i;
    if(ix>iy){
        i=ix;
	ix=iy;
	iy=i;
    }
    
    return _nboundary.get_data(ix*_dim+iy);
}

double chisquared::get_width(int ic, int ix){
  
    return _widths.get_data(ic,ix);
   
}

double chisquared::get_center(int ic, int ix){
    return _centers.get_data(ic,ix);
}

double chisquared::get_real_center(int ic, int ix){
    if(ic>=_ncenters || ix>=_dim){
        return exception_value;
    }
    
    int i;
    double ans=0.0;
    for(i=0;i<_dim;i++){
        ans+=_centers.get_data(ic,i)*_bases.get_data(i,ix);
    }
    return ans;
    
}

double chisquared::distance_to_center(int ic, array_1d<double> &pt){
    if(ic<0 || ic>=_centers.get_rows()){
        printf("WARHNG asked for center %d but max %d\n",
	ic,_centers.get_rows());
	
	exit(1);
    }
    
    array_1d<double> projected;
    
    int i;
    
    for(i=0;i<_centers.get_cols();i++){
        projected.set(i,project_to_basis(i,pt));
    }
    
    double dd;
    dd=0.0;
    for(i=0;i<_centers.get_cols();i++){
        dd+=power(projected.get_data(i)-_centers.get_data(ic,i),2);
    }
    dd=sqrt(dd);
    return dd;
}

double chisquared::get_boundary(int ix, int iy, int ipt, int idim){

    int i;
    if(ix>iy){
        i=ix;
	ix=iy;
	iy=i;
    }
    
    if(ix<0 || iy<0 || ix>=_dim || iy>=_dim){
        printf("WARNING asked for boundary slog %d %d but dim %d\n",ix,iy,_dim);
    }
    
    if(idim>=3 || ipt>=_nboundary.get_data(ix*_dim+iy))return exception_value;
    return _boundary[ix*_dim+iy][ipt][idim];
}

void chisquared::death_knell(char *word)const{
    printf("%s\n",word);
    exit(1);
}

int chisquared::get_dim(){
    return _dim;
}

int chisquared::get_ncenters(){
    return _ncenters;
}

chisquared::chisquared(){
    death_knell("meaningless constructor");
};

chisquared::chisquared(int id){
    _ncenters=1;
    _dim=id;
    
    _time_spent=0.0;
    
    _boundary=NULL;
    _dice=NULL;
    
    _rr_max=-1.0;
    _called=0;
    
    allot_arrays();
};

chisquared::chisquared(int id, int ic){
    _dim=id;
    _ncenters=ic;
    
    _time_spent=0.0;
 
    _boundary=NULL;   
    _dice=NULL;
     
    _rr_max=-1.0;
    _called=0;
    
    allot_arrays();
};



chisquared::~chisquared(){
    int i,ix,iy;
    
    if(_boundary!=NULL){
        for(ix=0;ix<_dim;ix++){
	    for(iy=ix+1;iy<_dim;iy++){
	            for(i=0;i<_boundary_room.get_data(ix*_dim+iy);i++){
		        delete [] _boundary[ix*_dim+iy][i];
		    }
		    delete [] _boundary[ix*_dim+iy];
		
	    }
	}
	delete [] _boundary;
    }
    
    if(_dice!=NULL){
        delete _dice;
    }

}

int chisquared::get_called(){
    return _called;
}

void chisquared::reset_timer(){
    _called=0;
    _time_spent=0.0;
}

void chisquared::decrement_called(){
    _called--;
}

double chisquared::operator()(array_1d<double> &v){
    death_knell("meaningless operator");
    return -1.0;
}

void chisquared::build_boundary(double rr){
    death_knell("meaningless build_boundary");
}

void chisquared::get_basis(int ix, array_1d<double> &v){
    int i;
    if(ix<_dim){
        for(i=0;i<_dim;i++)v.set(i,_bases.get_data(ix,i));
    }
    else{
        printf("WARNING called get_basis with %d %d\n",ix,_dim);
	exit(1);
    }
}

double chisquared::project_to_basis(int ix, array_1d<double> &vv) const{
    int i;
    double nn=1.0e30;
    if(ix<_dim){
        nn=0.0;
	for(i=0;i<_dim;i++)nn+=vv.get_data(i)*_bases.get_data(ix,i);
    }
    else{
        printf("WARNING called project_to_basis with %d %d\n",ix,_dim);
    }
    
    return nn;
}

s_curve::~s_curve(){}

s_curve::s_curve() : chisquared(6), _trig_factor(10.0){make_bases(22);
}

s_curve::s_curve(int id) : chisquared(id), _trig_factor(10.0){make_bases(22);}

s_curve::s_curve(int id, int ic) : chisquared(id,ic), _trig_factor(10.0){
        make_bases(22);
        _widths.set(0,0,0.5);
        _widths.set(0,1,0.5);
        
        if(ic>1)_widths.multiply_val(1,1,0.25);
        if(ic>2)_widths.multiply_val(2,0,0.25);
}


ellipses::~ellipses(){}

ellipses::ellipses() : chisquared(22){make_bases(13);}

ellipses::ellipses(int id) : chisquared(id){make_bases(13);}

ellipses::ellipses(int id, int ic) : chisquared(id,ic){
    //printf("constructed _dim %d _ncenters %d\n",_dim,_ncenters);
    make_bases(13);
}

ellipses_integrable::~ellipses_integrable(){}

ellipses_integrable::ellipses_integrable() : ellipses(){
    make_bases(-13);
    
    /*int i,j;
    for(i=1;i<_ncenters;i++){
        for(j=0;j<_dim;j++){
            _widths.set(i,j,_widths.get_data(0,j));
        }
    }*/

}

ellipses_integrable::ellipses_integrable(int id) : ellipses(id){
    make_bases(-13);

    /*int i,j;
    for(i=1;i<_ncenters;i++){
        for(j=0;j<_dim;j++){
            _widths.set(i,j,_widths.get_data(0,j));
        }
    }*/
    
}

ellipses_integrable::ellipses_integrable(int id, int ic) : ellipses(id,ic){
    make_bases(-13);

    /*int i,j;
    for(i=1;i<_ncenters;i++){
        for(j=0;j<_dim;j++){
            _widths.set(i,j,_widths.get_data(0,j));
        }
    }*/


}

linear_ellipses::~linear_ellipses(){}

linear_ellipses::linear_ellipses() : chisquared(22){make_bases(17);}

linear_ellipses::linear_ellipses(int id) : chisquared(id){make_bases(17);}

linear_ellipses::linear_ellipses(int id, int ic) : chisquared(id,ic){
    make_bases(17);
}

void s_curve::get_trough_points(array_2d<double> &outpoints){
    /*
    Get the points against which we will measure performance of APS that are near
    the centers of low chisquared
    */
    
    //first do the points along the S curve
    array_1d<double> a0;
    int i,j;
    for(i=0;i<_dim;i++){
        a0.set(i,_centers.get_data(0,i));
    }
    
    array_1d<double> vv;
    double theta,chival,xth,yth;
    
    for(theta=-1.0*pi;theta<=1.0*pi+0.001;theta+=0.4*pi){
        xth=_trig_factor*sin(theta)+_centers.get_data(0,0);
        if(theta!=0.0){
	    yth=_trig_factor*theta*(cos(theta)-1.0)/fabs(theta)+_centers.get_data(0,1);
        }
        else{
            yth=_centers.get_data(0,1);
        }
	a0.set(0,xth);
        a0.set(1,yth);
        
        for(i=0;i<_dim;i++){
            vv.set(i,0.0);
            for(j=0;j<_dim;j++){
                vv.add_val(i,a0.get_data(j)*_bases.get_data(j,i));
            }
        }
        chival=(*this)(vv);
        if(chival>0.01){
            printf("WARNING in S trough chival %e\n",chival);
            exit(1);
        }
        
        outpoints.add_row(vv);
    }
    
    int k;
    for(i=1;i<_ncenters;i++){
        for(j=0;j<_dim;j++){
            vv.set(j,0.0);
            for(k=0;k<_dim;k++){
                vv.add_val(j,_centers.get_data(i,k)*_bases.get_data(k,j));
            }
        }
        
        chival=(*this)(vv);
        if(chival>0.01){
            printf("WARNING in center %d trough chival %e\n",i,chival);
            exit(1);
        }
        outpoints.add_row(vv);
    }
    
    _called=0;

}

void s_curve::get_border_points(array_2d<double> &outpoints){

    if(_dice==NULL){
        death_knell("you called build_boundary before making bases");
    }

    int ix,iy,ic,ir,i,j;
    
    double tol=0.5,br=33.93;
    double theta,dxdth,dydth,x0,y0;
    double grad[2],norm,dfabsdth,ds,chitest;
    
    array_1d<double> pt,alpha;
    
    alpha.set_dim(_dim);
    pt.set_dim(_dim);
    
    alpha.set_name("s_curve_border_alpha");
    pt.set_name("s_curve_border_pt");
    
    _centers.set_where("s_curve_border");
    _bases.set_where("s_curve_border");
    _widths.set_where("s_curve_border");
    
    ///////////////below we will do ix=0, iy=1 for the 0th center (the S curve itself)
    ix=0;
    iy=1;
    
    array_1d<double> alphaUp,alphaDown,alphaNearest,vvNearest;
    double fUp,fDown,fNearest,ftrial;
    double fWorst,err;
    
    fWorst=-1.0;
    
    if(_widths.get_data(0,0)<_widths.get_data(0,1))ds=0.1*_widths.get_data(0,0);
    else ds=0.1*_widths.get_data(0,1);
    
    for(theta=-1.0*pi;theta<=1.01*pi;theta+=0.5*pi){
        if(theta<0.0)dfabsdth=-1.0;
        else dfabsdth=1.0;
    
        x0=_centers.get_data(0,0)+_trig_factor*sin(theta);
        if(theta!=0.0){
            y0=_centers.get_data(0,1)+_trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
        }
        else{
            y0=_centers.get_data(0,1);
        }
        
        for(ir=0;ir<2;ir++){

            dxdth=_trig_factor*cos(theta);
            
            if(theta!=0.0){
	        dydth=_trig_factor*((cos(theta)-1.0)/fabs(theta)-sin(theta)*theta/fabs(theta)
	                   -(cos(theta)-1.0)*dfabsdth/theta);
            }
            else{
                dydth=0.0;
            }
	
            if(ir==0){
	        grad[0]=dydth;
	        grad[1]=-1.0*dxdth;
	    }
	    else{
	        grad[0]=-1.0*dydth;
	        grad[1]=dxdth;
	    }
            norm=grad[0]*grad[0]+grad[1]*grad[1];
	    norm=sqrt(norm);
	
	    for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(0,i));
	    
            for(i=2;i<_dim;i++){
                alphaUp.set(i,_centers.get_data(0,i));
                alphaDown.set(i,_centers.get_data(0,i));
            }
            
            alphaDown.set(0,x0);
            alphaDown.set(1,y0);
            alphaUp.set(0,x0);
            alphaUp.set(1,y0);
            for(i=0;i<_dim;i++){
                pt.set(i,0.0);
                for(j=0;j<_dim;j++)pt.add_val(i,alphaDown.get_data(j)*_bases.get_data(j,i));
            }
            fDown=(*this)(pt);
            if(fDown>br){
                printf("WARNING fDown %e br %e\n",fDown,br);
                exit(1);
            }
            
            fUp=fDown;
            while(fUp<=br){
                alphaUp.add_val(0,grad[0]*ds/norm);
                alphaUp.add_val(1,grad[1]*ds/norm);
                
                for(i=0;i<_dim;i++){
                    pt.set(i,0.0);
                    for(j=0;j<_dim;j++){
                        pt.add_val(i,alphaUp.get_data(j)*_bases.get_data(j,i));
                    }
                }
                fUp=(*this)(pt);
            }
            
            if(fUp-br<br-fDown){
                fNearest=fUp;
                for(i=0;i<_dim;i++){
                    alphaNearest.set(i,alphaUp.get_data(i));
                    vvNearest.set(i,pt.get_data(i));
                }
            }
            else{
                fNearest=fDown;
                for(i=0;i<_dim;i++){
                    alphaNearest.set(i,alphaDown.get_data(i));
                    vvNearest.set(i,pt.get_data(i));
                }
            }
            
            while(fabs(br-fNearest)>tol){
                for(i=0;i<_dim;i++)alpha.set(i,0.5*(alphaUp.get_data(i)+alphaDown.get_data(i)));
                for(i=0;i<_dim;i++){
                    pt.set(i,0.0);
                    for(j=0;j<_dim;j++){
                        pt.add_val(i,alpha.get_data(j)*_bases.get_data(j,i));
                    }
                }
                
                ftrial=(*this)(pt);
                
                if(ftrial<br){
                    for(i=0;i<_dim;i++)alphaDown.set(i,alpha.get_data(i));
                }
                else{
                    for(i=0;i<_dim;i++)alphaUp.set(i,alpha.get_data(i));
                }
                
                if(fabs(br-ftrial)<fabs(br-fNearest)){
                    fNearest=ftrial;
                    for(i=0;i<_dim;i++){
                        alphaNearest.set(i,alpha.get_data(i));
                        vvNearest.set(i,pt.get_data(i));
                    }
                }
            }
            
            err=fabs(fNearest-br);
            if(err>fWorst)fWorst=err;
            
	    outpoints.add_row(vvNearest);
		
        }
    }
    
    printf("after 0,1 fworst %e -- %d\n",fWorst,outpoints.get_rows());
    
    for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(0,i));
    
    double th,s2x,s2y,aa,rr;
    int iabort,abort_max=20;
    for(theta=-1.0*pi;theta<=1.5*pi;theta+=2.0*pi){
        for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(0,i));
    
        x0=_centers.get_data(0,0)+_trig_factor*sin(theta);
        if(theta!=0.0){
            y0=_centers.get_data(0,1)+_trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
        }
        else{
            y0=_centers.get_data(0,1);
        }
        
        s2x=_widths.get_data(0,0)*_widths.get_data(0,0);
        s2y=_widths.get_data(0,1)*_widths.get_data(0,1);
        
	if(theta<0.0){
	    th=0.0;
	}
	else{
	   th=pi;
	}

        aa=cos(th)*cos(th)/s2x+sin(th)*sin(th)/s2y;
	rr=br/aa;
	rr=sqrt(rr);
            
        iabort=0;
        chitest=br+10.0*tol;
        while(iabort<abort_max && fabs(chitest-br)>tol){
            
	    alpha.set(0,x0+rr*cos(th));
	    alpha.set(1,y0+rr*sin(th));
	
	    for(i=0;i<_dim;i++){
		pt.set(i,0.0);
	        for(j=0;j<_dim;j++)pt.add_val(i,alpha.get_data(j)*_bases.get_data(j,i));
	    }
    
            chitest=(*this)(pt);
	    
	    if(fabs(chitest-br)<tol){
                err=fabs(chitest-br);
                if(err>fWorst)fWorst=err;
	        outpoints.add_row(pt);
	    }
                
            if(chitest>br)rr*=0.95;
            else rr*=1.02;
            iabort++;
        }
        err=fabs(chitest-br);
        if(err>fWorst)fWorst=err;

    }
    
    printf("after caps fWorst %e -- %d\n",fWorst,outpoints.get_rows());
    
    /////////////////////now do ix=0,1; iy>=2
    double xx,yy,yth,xth;

    for(iy=2;iy<_dim;iy++){for(ir=0;ir<2;ir++){
        for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(0,i));
        for(theta=-1.0*pi;theta<=1.0*pi;theta+=0.5*pi){
           alpha.set(0,_centers.get_data(0,0)+_trig_factor*sin(theta));
           
           if(theta!=0.0){
               alpha.set(1,_centers.get_data(0,1)
                   +_trig_factor*theta*(cos(theta)-1.0)/fabs(theta));
       
           }
           else{
               alpha.set(1,_centers.get_data(0,1));
           }
           
           if(ir==0)alpha.set(iy,_centers.get_data(0,iy)+sqrt(br)*_widths.get_data(0,iy));
           else alpha.set(iy,_centers.get_data(0,iy)-sqrt(br)*_widths.get_data(0,iy));
       
           for(i=0;i<_dim;i++){
	       pt.set(i,0.0);
	       for(j=0;j<_dim;j++)pt.add_val(i,alpha.get_data(j)*_bases.get_data(j,i));
           } 
           chitest=(*this)(pt);
           
	   if(fabs(chitest-br)<5.0){
               err=fabs(chitest-br);
               if(err>fWorst)fWorst=err;
               outpoints.add_row(pt);
	   }
        }
    }}
   
    printf("after (0,1), >2 fWorst %e -- %d\n",fWorst,outpoints.get_rows());
    
    double thetamin,thetamax,dtheta,newtheta;

    for(ix=0;ix<2;ix++){
        if(ix==1){
            thetamin=-1.0*pi;
	    thetamax=1.5*pi;
	    dtheta=2.0*pi;
        }
        else{
            thetamin=-0.5*pi;
	    thetamax=0.6*pi;
	    dtheta=1.0*pi;
    
        }


        s2x=_widths.get_data(0,ix)*_widths.get_data(0,ix);
        for(iy=2;iy<_dim;iy++){
            s2y=_widths.get_data(0,iy)*_widths.get_data(0,iy);
            for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(0,i));
            for(theta=thetamin;theta<thetamax;theta+=dtheta){
	     
	         if(ix==0){
	             xth=_centers.get_data(0,0)+_trig_factor*sin(theta);
                     if(theta!=0.0){
	                 alpha.set(1,_centers.get_data(0,1)
		            +_trig_factor*theta*(cos(theta)-1.0)/fabs(theta));
                     }
                     else{
                         alpha.set(1,_centers.get_data(0,1));
                     }
	         }
	         else{
                     if(theta!=0.0){
	                 xth=_centers.get_data(0,1)+_trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
                     }
                     else{
                         xth=_centers.get_data(0,1);
                     }
	             alpha.set(0,_centers.get_data(0,0)+_trig_factor*sin(theta));
	         }
	     
	     
	         for(th=0.0;th<1.99*pi;th+=0.5*pi){
	             aa=(cos(th)*cos(th)/s2x+sin(th)*sin(th)/s2y);
		     rr=br/aa;
		     rr=sqrt(rr);
		 
		     alpha.set(iy,_centers.get_data(0,iy)+rr*sin(th));
		     alpha.set(ix,xth+rr*cos(th));
	         
                     if(ix==0){
                         newtheta=find_theta_from_x(alpha.get_data(ix));
                         if(newtheta!=0.0){
                             alpha.set(1,_centers.get_data(0,1)+_trig_factor*newtheta*(cos(newtheta)-1.0)/fabs(newtheta));
                         }
                         else{
                             alpha.set(1,_centers.get_data(0,1));
                         }
                     }
                     else{
                         newtheta=find_theta_from_y(alpha.get_data(ix));
                         alpha.set(0,_centers.get_data(0,0)+_trig_factor*sin(newtheta));
                     }
		 
		     for(i=0;i<_dim;i++){
                         pt.set(i,0.0);
		         for(j=0;j<_dim;j++){
		             pt.add_val(i,alpha.get_data(j)*_bases.get_data(j,i));
		         }
		     }
		     chitest=(*this)(pt);
		     if(fabs(chitest-br)<0.1){
                         err=fabs(chitest-br);
                         if(err>fWorst)fWorst=err;
		         outpoints.add_row(pt);
		     }
		 }
	         
	     
	    }
        }
    }
    
    printf("after caps fWorst %e -- %d\n",fWorst,outpoints.get_rows());
    
    for(ic=0;ic<_ncenters;ic++){
        for(ix=0;ix<_dim;ix++){
	    for(iy=ix+1;iy<_dim;iy++){
	        if(ic>0 || ix>1){
		    
		    for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(ic,i));
		    s2x=power(_widths.get_data(ic,ix),2);
		    s2y=power(_widths.get_data(ic,iy),2);
		    for(theta=0.0;theta<=1.9*pi;theta+=0.5*pi){
		        aa=power(cos(theta),2)/s2x+power(sin(theta),2)/s2y;
			rr=br/aa;
			rr=sqrt(rr);
			alpha.set(ix,_centers.get_data(ic,ix)+rr*cos(theta));
			alpha.set(iy,_centers.get_data(ic,iy)+rr*sin(theta));
		        
			for(i=0;i<_dim;i++){
			    pt.set(i,0.0);
			    for(j=0;j<_dim;j++)pt.add_val(i,alpha.get_data(j)*_bases.get_data(j,i));
			}
			
			chitest=(*this)(pt);
			if(fabs(chitest-br)<5.0){
			    //if(ic==1 && iy==1)printf("adding %d %d %d\n",ix,iy,ic);
			    
                            err=fabs(chitest-br);
                            if(err>fWorst)fWorst=err;
			    outpoints.add_row(pt);
			}
		    }
		
		}
	    }
	}
    }
    printf("and finally fWorst %e -- %d\n",fWorst,outpoints.get_rows());
    _centers.set_where("nowhere");
    _widths.set_where("nowhere");
    _bases.set_where("nowhere");
  
    
}

double s_curve::operator()(array_1d<double> &in_pt){
    
    if(_dice==NULL){
         death_knell("you called operator before making bases");
    } 
    
    
    double before=double(time(NULL));
    
    _called++;

    
    array_1d<double> dd,pt;
    dd.set_dim(_ncenters);
    pt.set_dim(_dim);
    
    pt.set_name("s_curve_operator_pt");
    dd.set_name("s_curve_operator_dd");
    
    _centers.set_where("s_curve_operator");
    _widths.set_where("s_curve_operator");
    _bases.set_where("s_curve_operator");
    
    
    double theta,xth,yth,dth,dthmin;
    int i;
    
    for(i=0;i<_dim;i++)pt.set(i,project_to_basis(i,in_pt));
    
    //for(i=0;i<_dim;i++)printf("%e %e\n",pt[i],_centers[0][i]);
    dd.set(0,0.0);
    for(i=2;i<_dim;i++)dd.add_val(0,power((pt.get_data(i)-_centers.get_data(0,i))/_widths.get_data(0,i),2));
    
    dthmin=-1.0;
    for(theta=-1.0*pi;theta<=1.0*pi;theta+=0.01){
        xth=_trig_factor*sin(theta)+_centers.get_data(0,0);
        if(theta!=0.0){
	    yth=_trig_factor*theta*(cos(theta)-1.0)/fabs(theta)+_centers.get_data(0,1);
        }
        else{
            yth=_centers.get_data(0,1);
        }
	dth=power((xth-pt.get_data(0))/_widths.get_data(0,0),2)+power((yth-pt.get_data(1))/_widths.get_data(0,1),2);
	
	if(dthmin<0.0 || dth<dthmin)dthmin=dth;
    }
    dd.add_val(0,dthmin);
    
    int ii;
    for(ii=1;ii<_ncenters;ii++){
        dd.set(ii,0.0);
	for(i=0;i<_dim;i++){
	    dd.add_val(ii,power((pt.get_data(i)-_centers.get_data(ii,i))/_widths.get_data(ii,i),2));
	}
    }
    
    double ddmin;
    for(ii=0;ii<_ncenters;ii++){
        if(ii==0 || dd.get_data(ii)<ddmin)ddmin=dd.get_data(ii);
    }
    
    _centers.set_where("nowhere");
    _bases.set_where("nowhere");
    _widths.set_where("nowhere");
    
    _time_spent+=double(time(NULL))-before;
    
    return ddmin;
}

double s_curve::distance_to_center(int ic, array_1d<double> &in_pt) {
    
    if(_dice==NULL){
         death_knell("you called operator before making bases");
    } 
    
    if(ic<0 || ic>=_centers.get_rows()){
        printf("WARNING asked for center %d but max %d\n",
	ic,_centers.get_rows());
	
	exit(1);
    }
    

    
    array_1d<double> pt;
    
    pt.set_dim(_dim);
    
    pt.set_name("s_curve_distance_to_center_pt");
   
    _centers.set_where("s_curve_operator");
    _widths.set_where("s_curve_operator");
    _bases.set_where("s_curve_operator");
    
    
    double theta,xth,yth,dth,dthmin,dd;
    int i;
    
    for(i=0;i<_dim;i++)pt.set(i,project_to_basis(i,in_pt));
    
    
    //for(i=0;i<_dim;i++)printf("%e %e\n",pt[i],_centers[0][i]);
    dd=0.0;
    
    if(ic==0){
        for(i=2;i<_dim;i++)dd+=power((pt.get_data(i)-_centers.get_data(0,i))/_widths.get_data(0,i),2);
    
        dthmin=-1.0;
        for(theta=-1.0*pi;theta<=1.0*pi;theta+=0.01){
            xth=_trig_factor*sin(theta)+_centers.get_data(0,0);
            if(theta!=0.0){
	        yth=_trig_factor*theta*(cos(theta)-1.0)/fabs(theta)+_centers.get_data(0,1);
            }
            else{
                yth=_centers.get_data(0,1);
            }
	    dth=power((xth-pt.get_data(0))/_widths.get_data(0,0),2)+power((yth-pt.get_data(1))/_widths.get_data(0,1),2);
	
	    if(dthmin<0.0 || dth<dthmin)dthmin=dth;
        }
        dd+=dth;
    }
    else{

	for(i=0;i<_dim;i++){
	    dd+=power((pt.get_data(i)-_centers.get_data(ic,i))/_widths.get_data(ic,i),2);
	}
    
    }
    
    dd=sqrt(dd);
    
    _centers.set_where("nowhere");
    _bases.set_where("nowhere");
    _widths.set_where("nowhere");
    
    return dd;
}

double s_curve::find_theta_from_x(double x){

    double dx=x-_centers.get_data(0,0);
    
    double ratio=dx/_trig_factor;
    if(ratio>1.0)ratio=1.0;
    else if(ratio<-1.0)ratio=-1.0;
    
    double naive=asin(ratio);
    
    return naive;
}

double s_curve::find_theta_from_y(double y){
    double dy=y-_centers.get_data(0,1);
    double sgn;
    if(dy<0.0)sgn=1.0;
    else sgn=-1.0;
    
    double ratio=dy/_trig_factor;
    if(ratio>2.0)ratio=2.0;
    else if(ratio<-2.0)ratio=-2.0;
    
    double naive=acos(1.0-fabs(ratio));
    return sgn*naive;
}

void s_curve::build_boundary(double br){
    if(_dice==NULL){
        death_knell("you called build_boundary before making bases");
    }

    reset_boundary();

    int ix,iy,ic,ir,i,j;
    
    double tol=0.5;
    double theta,dxdth,dydth,x0,y0;
    double grad[2],norm,dfabsdth,ds,chitest;
    
    array_1d<double> pt,alpha;
    
    alpha.set_dim(_dim);
    pt.set_dim(_dim);
    
    alpha.set_name("s_curve_build_boundary_alpha");
    pt.set_name("s_curve_build_boundary_pt");
    
    _centers.set_where("s_curve_build_boundary");
    _bases.set_where("s_curve_build_boundary");
    _widths.set_where("s_curve_build_boundary");
    
    ///////////////below we will do ix=0, iy=1 for the 0th center (the S curve itself)
    ix=0;
    iy=1;
    
    array_1d<double> alphaUp,alphaDown,alphaNearest;
    double fUp,fDown,fNearest,ftrial;
    double fWorst,err;
    
    fWorst=-1.0;
    
    if(_widths.get_data(0,0)<_widths.get_data(0,1))ds=0.1*_widths.get_data(0,0);
    else ds=0.1*_widths.get_data(0,1);
    
    for(theta=-1.0*pi;theta<=1.0*pi;theta+=0.01){
        if(theta<0.0)dfabsdth=-1.0;
        else dfabsdth=1.0;
    
        x0=_centers.get_data(0,0)+_trig_factor*sin(theta);
        if(theta!=0.0){
            y0=_centers.get_data(0,1)+_trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
        }
        else{
            y0=_centers.get_data(0,1);
        }
    
        for(ir=0;ir<2;ir++){

            dxdth=_trig_factor*cos(theta);
            if(theta!=0.0){
	        dydth=_trig_factor*((cos(theta)-1.0)/fabs(theta)-sin(theta)*theta/fabs(theta)
	               -(cos(theta)-1.0)*dfabsdth/theta);
            }
            else{
                dydth=0.0;
            }
	
            if(ir==0){
	        grad[0]=dydth;
	        grad[1]=-1.0*dxdth;
	    }
	    else{
	        grad[0]=-1.0*dydth;
	        grad[1]=dxdth;
	    }
            norm=grad[0]*grad[0]+grad[1]*grad[1];
	    norm=sqrt(norm);
	
	    for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(0,i));
	    
            for(i=2;i<_dim;i++){
                alphaUp.set(i,_centers.get_data(0,i));
                alphaDown.set(i,_centers.get_data(0,i));
            }
            
            alphaDown.set(0,x0);
            alphaDown.set(1,y0);
            alphaUp.set(0,x0);
            alphaUp.set(1,y0);
            for(i=0;i<_dim;i++){
                pt.set(i,0.0);
                for(j=0;j<_dim;j++)pt.add_val(i,alphaDown.get_data(j)*_bases.get_data(j,i));
            }
            fDown=(*this)(pt);
            if(fDown>br){
                printf("WARNING fDown %e br %e\n",fDown,br);
                exit(1);
            }
            
            fUp=fDown;
            while(fUp<=br){
                alphaUp.add_val(0,grad[0]*ds/norm);
                alphaUp.add_val(1,grad[1]*ds/norm);
                
                for(i=0;i<_dim;i++){
                    pt.set(i,0.0);
                    for(j=0;j<_dim;j++){
                        pt.add_val(i,alphaUp.get_data(j)*_bases.get_data(j,i));
                    }
                }
                fUp=(*this)(pt);
            }
            
            if(fUp-br<br-fDown){
                fNearest=fUp;
                for(i=0;i<_dim;i++)alphaNearest.set(i,alphaUp.get_data(i));
            }
            else{
                fNearest=fDown;
                for(i=0;i<_dim;i++)alphaNearest.set(i,alphaDown.get_data(i));
            }
            
            while(fabs(br-fNearest)>tol){
                for(i=0;i<_dim;i++)alpha.set(i,0.5*(alphaUp.get_data(i)+alphaDown.get_data(i)));
                for(i=0;i<_dim;i++){
                    pt.set(i,0.0);
                    for(j=0;j<_dim;j++){
                        pt.add_val(i,alpha.get_data(j)*_bases.get_data(j,i));
                    }
                }
                
                ftrial=(*this)(pt);
                
                if(ftrial<br){
                    for(i=0;i<_dim;i++)alphaDown.set(i,alpha.get_data(i));
                }
                else{
                    for(i=0;i<_dim;i++)alphaUp.set(i,alpha.get_data(i));
                }
                
                if(fabs(br-ftrial)<fabs(br-fNearest)){
                    fNearest=ftrial;
                    for(i=0;i<_dim;i++)alphaNearest.set(i,alpha.get_data(i));
                }
            }
            
            err=fabs(fNearest-br);
            if(err>fWorst)fWorst=err;
            add_to_boundary(alphaNearest,ix,iy,fNearest);
	    
            
	    /*for(i=0;i<2;i++){
	        printf("%e ",alpha[i]);
	    }
	    printf("%e\n",chitest);*/
	    
		
        }
    }
    
    printf("after 0,1 fworst %e\n",fWorst);
    
    for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(0,i));
    
    double th,s2x,s2y,aa,rr,thmin,thmax;
    int iabort,abort_max=20;
    for(theta=-1.0*pi;theta<=1.5*pi;theta+=2.0*pi){
        for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(0,i));
    
        x0=_centers.get_data(0,0)+_trig_factor*sin(theta);
        if(theta!=0.0){
            y0=_centers.get_data(0,1)+_trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
        }
        else{
            y0=_centers.get_data(0,1);
        }
        
        s2x=_widths.get_data(0,0)*_widths.get_data(0,0);
        s2y=_widths.get_data(0,1)*_widths.get_data(0,1);
        
	if(theta<0.0){
	    thmin=-0.5*pi;
	    thmax=0.5*pi;
	}
	else{
	   thmin=0.5*pi;
	   thmax=1.5*pi;
	}
	
        for(th=thmin;th<thmax;th+=0.01){
        
            aa=cos(th)*cos(th)/s2x+sin(th)*sin(th)/s2y;
	    rr=br/aa;
	    rr=sqrt(rr);
            
            iabort=0;
            chitest=br+10.0*tol;
            while(iabort<abort_max && fabs(chitest-br)>tol){
            
	        alpha.set(0,x0+rr*cos(th));
	        alpha.set(1,y0+rr*sin(th));
	
	        for(i=0;i<_dim;i++){
		    pt.set(i,0.0);
	            for(j=0;j<_dim;j++)pt.add_val(i,alpha.get_data(j)*_bases.get_data(j,i));
	        }
    
                chitest=(*this)(pt);
	    
	        if(fabs(chitest-br)<tol){
                    err=fabs(chitest-br);
                    if(err>fWorst)fWorst=err;
	            add_to_boundary(alpha,ix,iy,chitest);
	        }
                
                if(chitest>br)rr*=0.95;
                else rr*=1.02;
                iabort++;
            }
            err=fabs(chitest-br);
            if(err>fWorst)fWorst=err;
	    /*for(i=0;i<2;i++)printf("%e ",alpha[i]);
	    printf("%e\n",chitest);*/
        }

    }
    
    printf("after caps fWorst %e\n",fWorst);
    
    /////////////////////now do ix=0,1; iy>=2
    double xx,yy,yth,xth;

    for(iy=2;iy<_dim;iy++){for(ir=0;ir<2;ir++){
        for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(0,i));
        for(theta=-1.0*pi;theta<=1.0*pi;theta+=0.01){
           alpha.set(0,_centers.get_data(0,0)+_trig_factor*sin(theta));
           if(theta!=0.0){
               alpha.set(1,_centers.get_data(0,1)
                   +_trig_factor*theta*(cos(theta)-1.0)/fabs(theta));
           }
           else{
               alpha.set(1,_centers.get_data(0,1));
           }
       
           if(ir==0)alpha.set(iy,_centers.get_data(0,iy)+sqrt(br)*_widths.get_data(0,iy));
           else alpha.set(iy,_centers.get_data(0,iy)-sqrt(br)*_widths.get_data(0,iy));
       
           for(i=0;i<_dim;i++){
	       pt.set(i,0.0);
	       for(j=0;j<_dim;j++)pt.add_val(i,alpha.get_data(j)*_bases.get_data(j,i));
           } 
           chitest=(*this)(pt);
           
	   if(fabs(chitest-br)<5.0){
               err=fabs(chitest-br);
               if(err>fWorst)fWorst=err;
	       add_to_boundary(alpha,0,iy,chitest);
	       add_to_boundary(alpha,1,iy,chitest);
	   }
        }
    }}
   
    printf("after (0,1), >2 fWorst %e\n",fWorst);
    
    double thetamin,thetamax,dtheta,newtheta;

    for(ix=0;ix<2;ix++){
        if(ix==1){
            thetamin=-1.0*pi;
	    thetamax=1.5*pi;
	    dtheta=2.0*pi;
        }
        else{
            thetamin=-0.5*pi;
	    thetamax=0.6*pi;
	    dtheta=1.0*pi;
    
        }


        s2x=_widths.get_data(0,ix)*_widths.get_data(0,ix);
        for(iy=2;iy<_dim;iy++){
            s2y=_widths.get_data(0,iy)*_widths.get_data(0,iy);
            for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(0,i));
            for(theta=thetamin;theta<thetamax;theta+=dtheta){
	     
	         if(ix==0){
	             xth=_centers.get_data(0,0)+_trig_factor*sin(theta);
                     if(theta!=0.0){
	                 alpha.set(1,_centers.get_data(0,1)
		            +_trig_factor*theta*(cos(theta)-1.0)/fabs(theta));
                     }
                     else{
                         alpha.set(1,_centers.get_data(0,1));
                     }
	         }
	         else{
                     if(theta!=0.0){
	                 xth=_centers.get_data(0,1)+_trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
                     }
                     else{
                         xth=_centers.get_data(0,1);
                     }
                     
	             alpha.set(0,_centers.get_data(0,0)+_trig_factor*sin(theta));
	         }
	     
	     
	         for(th=0.0;th<2.0*pi;th+=0.01){
	             aa=(cos(th)*cos(th)/s2x+sin(th)*sin(th)/s2y);
		     rr=br/aa;
		     rr=sqrt(rr);
		 
		     alpha.set(iy,_centers.get_data(0,iy)+rr*sin(th));
		     alpha.set(ix,xth+rr*cos(th));
	             
                     if(ix==0){
                         newtheta=find_theta_from_x(alpha.get_data(ix));
                         if(newtheta!=0.0){
                             alpha.set(1,_centers.get_data(0,1)+_trig_factor*newtheta*(cos(newtheta)-1.0)/fabs(newtheta));
                         }
                         else{
                             alpha.set(1,_centers.get_data(0,1));
                         }
                     }
                     else{
                         newtheta=find_theta_from_y(alpha.get_data(ix));
                         alpha.set(0,_centers.get_data(0,0)+_trig_factor*sin(newtheta));
                     }
                     
		 
		     for(i=0;i<_dim;i++){
			 pt.set(i,0.0);
		         for(j=0;j<_dim;j++){
		             pt.add_val(i,alpha.get_data(j)*_bases.get_data(j,i));
		         }
		     }
		     chitest=(*this)(pt);
		     if(fabs(chitest-br)<1.0){
                         err=fabs(chitest-br);
                         if(err>fWorst)fWorst=err;
		         add_to_boundary(alpha,ix,iy,chitest);
		     }
		 
	         }
	     
	    }
        }
    }
    
    printf("after caps fWorst %e\n",fWorst);
    
    for(ic=0;ic<_ncenters;ic++){
        for(ix=0;ix<_dim;ix++){
	    for(iy=ix+1;iy<_dim;iy++){
	        if(ic>0 || ix>1){
		    
		    for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(ic,i));
		    s2x=power(_widths.get_data(ic,ix),2);
		    s2y=power(_widths.get_data(ic,iy),2);
		    for(theta=0.0;theta<=2.0*pi;theta+=0.01){
		        aa=power(cos(theta),2)/s2x+power(sin(theta),2)/s2y;
			rr=br/aa;
			rr=sqrt(rr);
			alpha.set(ix,_centers.get_data(ic,ix)+rr*cos(theta));
			alpha.set(iy,_centers.get_data(ic,iy)+rr*sin(theta));
		        
			for(i=0;i<_dim;i++){
			    pt.set(i,0.0);
			    for(j=0;j<_dim;j++)pt.add_val(i,alpha.get_data(j)*_bases.get_data(j,i));
			}
			
			chitest=(*this)(pt);
			if(fabs(chitest-br)<5.0){
			    //if(ic==1 && iy==1)printf("adding %d %d %d\n",ix,iy,ic);
			    
                            err=fabs(chitest-br);
                            if(err>fWorst)fWorst=err;
			    add_to_boundary(alpha,ix,iy,chitest);
			}
		    }
		
		}
	    }
	}
    }
    printf("and finally fWorst %e\n",fWorst);
    _centers.set_where("nowhere");
    _widths.set_where("nowhere");
    _bases.set_where("nowhere");
  
    
}

double ellipses::operator()(array_1d<double> &in_pt){
    
    if(_dice==NULL){
         death_knell("you called operator before making bases");
    } 
    
    double before=double(time(NULL));
    
    int ii,ix;
    double dd,ddmin,nn;
    
    _centers.set_where("ellipse_operator");
    _widths.set_where("ellipse_operator");
    _bases.set_where("ellipse_operator");
    
    _called++;
    
    for(ii=0;ii<_ncenters;ii++){
        dd=0.0;
        for(ix=0;ix<_dim;ix++){
	    nn=project_to_basis(ix,in_pt);
	    dd+=power((_centers.get_data(ii,ix)-nn)/_widths.get_data(ii,ix),2);
	    //printf("%e\n",_centers[ii][ix]-nn);
	}
	if(ii==0 || dd<ddmin)ddmin=dd;
    }
    
    _centers.set_where("nowhere");
    _widths.set_where("nowhere");
    _bases.set_where("nowhere");
    
    _time_spent+=double(time(NULL))-before;
    
    return ddmin;
    
}

void ellipses::build_boundary(double br){

    if(_dice==NULL){
         death_knell("you called build_boundary before making bases");
    } 
    
    reset_boundary();
    
    double tol=0.5;
    int ix,iy,ic,i,j;
    double theta,rr,aa,chitest;
    
    /*
    printf("in build_boundary _centers are\n");
    for(i=0;i<_dim;i++){
        for(j=0;j<_ncenters;j++){
            printf("%e ",_centers.get_data(j,i));
        }
        printf("\n");
    }
    */
    
    array_1d<double> alpha,pt;
    alpha.set_dim(_dim);
    pt.set_dim(_dim);
    
    pt.set_name("ellipse_build_boundary_pt");
    alpha.set_name("ellipse_build_boundary_alpha");
    
    _centers.set_where("ellipse_build_boundary");
    _bases.set_where("ellipse_build_boundary");
    _widths.set_where("ellipse_build_boundary");
    
    
    for(ic=0;ic<_ncenters;ic++){
        for(ix=0;ix<_dim;ix++){
            for(iy=ix+1;iy<_dim;iy++){
	        for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(ic,i));
	        for(theta=0.0;theta<=2.0*pi;theta+=0.01){
	            aa=power(cos(theta)/_widths.get_data(ic,ix),2)+power(sin(theta)/_widths.get_data(ic,iy),2);
		    rr=br/aa;
		    rr=sqrt(rr);
		    alpha.set(ix,_centers.get_data(ic,ix)+rr*cos(theta));
		    alpha.set(iy,_centers.get_data(ic,iy)+rr*sin(theta));
		    
		    
		    for(i=0;i<_dim;i++){
		        pt.set(i,0.0);
			for(j=0;j<_dim;j++)pt.add_val(i,alpha.get_data(j)*_bases.get_data(j,i));
		    }
		    
		    chitest=(*this)(pt);
		    if(fabs(chitest-br)<tol){
		        add_to_boundary(alpha,ix,iy,chitest);
		    }
		    else{
		        printf("failed to add %d %e rr %e\n",ic,chitest,rr);
			exit(1);
		    }
	        }
	    }
        }
    }
    
    _centers.set_where("nowhere");
    _bases.set_where("nowhere");
    _widths.set_where("nowhere");
    
}

double linear_ellipses::operator()(array_1d<double> &in_pt){

    if(_dice==NULL){
         death_knell("you called operator before making bases");
    } 
    
    double before=double(time(NULL));
    
    int ii,ix;
    double dd,ddmin,nn;
    
    _centers.set_where("linear_ellipse_operator\n");
    _bases.set_where("linear_ellipse_operator\n");
    _widths.set_where("linear_ellipse_operator\n");
    
    _called++;
    
    for(ii=0;ii<_ncenters;ii++){
        dd=0.0;
	for(ix=0;ix<_dim;ix++){
	    nn=project_to_basis(ix,in_pt);
	    dd+=power((_centers.get_data(ii,ix)-nn)/_widths.get_data(ii,ix),2);
	}
	
	if(ii==0 || dd<ddmin)ddmin=dd;
    }
    
    _centers.set_where("nowhere");
    _bases.set_where("nowhere");
    _widths.set_where("nowhere");
    
    _time_spent+=double(time(NULL))-before;
    
    return sqrt(ddmin);
}

void linear_ellipses::build_boundary(double br){
    if(_dice==NULL){
         death_knell("you called build_boundary before making bases");
    } 
    
    reset_boundary();
    
    double tol=0.5;
    int ic,ix,iy,i,j;
 
    double theta,rr,aa,chitest;
    
    array_1d<double> alpha,pt;
    
    alpha.set_dim(_dim);
    pt.set_dim(_dim);
    
    alpha.set_name("linear_ellipse_build_boundary_alpha");
    pt.set_name("linear_ellipse_build_boundary_pt");
    
    _centers.set_where("linear_ellipse_build_boundary");
    _bases.set_where("linear_ellipse_build_boundary");
    _widths.set_where("linear_ellipse_build_boundary");
    
 
    
    for(ic=0;ic<_ncenters;ic++){
        for(ix=0;ix<_dim;ix++){
	    for(iy=ix+1;iy<_dim;iy++){
	        for(i=0;i<_dim;i++)alpha.set(i,_centers.get_data(ic,i));
		for(theta=0.0;theta<=2.0*pi;theta+=0.01){
		    aa=power(cos(theta)/_widths.get_data(ic,ix),2)+power(sin(theta)/_widths.get_data(ic,iy),2);
		    rr=br*br/aa;
		    rr=sqrt(rr);
		    
		    alpha.set(ix,_centers.get_data(ic,ix)+rr*cos(theta));
		    alpha.set(iy,_centers.get_data(ic,iy)+rr*sin(theta));
		    for(i=0;i<_dim;i++){
		        pt.set(i,0.0);
			for(j=0;j<_dim;j++)pt.add_val(i,alpha.get_data(j)*_bases.get_data(j,i));
		    }
		    chitest=(*this)(pt);
		    if(fabs(chitest-br)<tol){
		        add_to_boundary(alpha,ix,iy,chitest);
		    }
		    else{
		       printf("failed to add %d %e\n",ic,chitest);
		       exit(1);
		    }
		}
	    }
	}
    }
    
    _centers.set_where("nowhere");
    _bases.set_where("nowhere");
    _widths.set_where("nowhere");
    
}

void ellipses_integrable::integrate_boundary(int ix1, int ix2, double lim, char *filename){
    integrate_boundary(ix1,ix2,lim,filename,5.0);
}

void ellipses_integrable::integrate_boundary(int ix1, int ix2, double lim, char *filename, double spread){
    
    int i,j;
    /*
    printf("in integrate_boundary, centers are\n");
    for(i=0;i<_dim;i++){
        for(j=0;j<_ncenters;j++){
            printf("%e ",_centers.get_data(j,i));
        }
        printf("\n");
    }
    */
    
    for(i=0;i<_dim;i++){
        for(j=0;j<_dim;j++){
            if(i==j){
                if(fabs(_bases.get_data(i,j)-1.0)>1.0e-6){
                    printf("WARNING bases %d %d %e\n",i,j,_bases.get_data(i,j));
                    exit(1);
                }
            }
            else{
                if(fabs(_bases.get_data(i,j))>1.0e-6){
                    printf("WARNING bases %d %d %e\n",i,j,_bases.get_data(i,j));
                    exit(1);
                }
            }
        }
    }
    
    /*for(i=0;i<_dim;i++){
        for(j=1;j<_ncenters;j++){
            if(fabs(_widths.get_data(0,i)-_widths.get_data(j,i))/_widths.get_data(0,i)>1.0e-6){
                printf("WARNING _widths do not agree %e %e %d %d\n",
                _widths.get_data(0,i),_widths.get_data(j,i),j,i);
                exit(1);
            }
        }
    }*/
    
    int icenter,k,imin,ct=0,duplicate;
    array_1d<double> chiarr,xarr,yarr,trial;
    array_1d<int> dexes;
    double xx,yy,xwidth,ywidth,nn,mm,ddmin,ddshld,total=0.0;
    double xbound,ybound,chival,xstep,ystep;
    
    array_1d<double> l_marginalized_volume;
    
    for(icenter=0;icenter<_ncenters;icenter++){
        l_marginalized_volume.set(icenter,0.0);
        for(i=0;i<_dim;i++){
            if(i!=ix1 && i!=ix2){
                l_marginalized_volume.add_val(icenter,log(_widths.get_data(icenter,i)));
            }
        }
    }
    
    /*for(icenter=0;icenter<_ncenters;icenter++){
        printf("l_vol %e \n",l_marginalized_volume.get_data(icenter));
    }*/
    
    for(icenter=0;icenter<_ncenters;icenter++){
        if(icenter==0 || _widths.get_data(icenter,ix1)>xbound)xbound=_widths.get_data(icenter,ix1);
        if(icenter==0 || _widths.get_data(icenter,ix1)<xwidth)xwidth=_widths.get_data(icenter,ix1);
        
        if(icenter==0 || _widths.get_data(icenter,ix2)>ybound)ybound=_widths.get_data(icenter,ix2);
        if(icenter==0 || _widths.get_data(icenter,ix2)<ywidth)ywidth=_widths.get_data(icenter,ix2);
    }
    
    dexes.set_dim(100);
    xarr.set_dim(100);
    yarr.set_dim(100);
    chiarr.set_dim(100);
    
    array_1d<int> doubled_up;
    doubled_up.set_name("doubled_up");
    
    xstep=0.1*xwidth;
    ystep=0.1*ywidth;
     
    for(icenter=0;icenter<_ncenters;icenter++){
        //xwidth=_widths.get_data(icenter,ix1);
        //ywidth=_widths.get_data(icenter,ix2);
        
        //xbound=_widths.get_data(icenter,ix1);
        //ybound=_widths.get_data(icenter,ix2);
        
        for(xx=_centers.get_data(icenter,ix1)-spread*xbound;xx<_centers.get_data(icenter,ix1)+(spread+0.1)*xbound;xx+=xstep){
            for(yy=_centers.get_data(icenter,ix2)-spread*ybound;yy<_centers.get_data(icenter,ix2)+(spread+0.1)*ybound;yy+=ystep){
                
                for(k=0;k<_ncenters;k++){
                    for(i=0;i<_dim;i++)trial.set(i,_centers.get_data(k,i));
                    trial.set(ix1,xx);
                    trial.set(ix2,yy);
                    nn=0.0;
                    for(i=0;i<_dim;i++){
                        nn+=power((trial.get_data(i)-_centers.get_data(k,i))/_widths.get_data(k,i),2);
                    }
                    
                    if(k==0 || nn<ddmin){
                        ddmin=nn;
                        imin=k;
                    }
                    if(k==icenter)ddshld=nn;
                }
                
                /*if(imin!=icenter){
                    
                    printf("dim %d %d\n",ix1,ix2);
                    printf("CURIOUS icenter %d imin %d\n",imin,icenter);
                    
                    printf("%e %e -- %e %e %e %e -- %e %e %e %e\n",
                    xx,yy,
                    _centers.get_data(icenter,ix1),_widths.get_data(icenter,ix1),
                    _centers.get_data(icenter,ix2),_widths.get_data(icenter,ix2),
                    _centers.get_data(imin,ix1),_widths.get_data(imin,ix1),
                    _centers.get_data(imin,ix2),_widths.get_data(imin,ix2));
                    
                    printf("%e %e\n",ddmin,ddshld);
                    
                    
                    exit(1);
                }*/
                
                for(i=0;i<_dim;i++)trial.set(i,_centers.get_data(icenter,i));
                trial.set(ix1,xx);
                trial.set(ix2,yy);
                
                chival=(*this)(trial)-2.0*l_marginalized_volume.get_data(icenter);
                
                
                if(icenter==imin){
                    duplicate=-1;
                }
                else{
                    duplicate=-1;
                    for(i=0;i<xarr.get_dim() && duplicate<0;i++){
                        if(fabs(xx-xarr.get_data(i))<xstep && fabs(yy-yarr.get_data(i))<ystep){
                            duplicate=i;
                        }
                    }
               
                }
                
                if(duplicate<0){
                    xarr.set(ct,xx);
                    yarr.set(ct,yy);
                    chiarr.set(ct,chival);
                    dexes.set(ct,ct);
                    ct++;
                }
                else{
                    nn=exp(-0.5*chiarr.get_data(duplicate))+exp(-0.5*chival);
                    chiarr.set(duplicate,-2.0*log(nn));
                
                }
                total+=exp(-0.5*chival);
                
                if(dexes.get_room()==ct){
                   //printf("adding room %d %d\n",dexes.get_dim(),ct);
                   dexes.add_room(1000);
                   xarr.add_room(1000);
                   yarr.add_room(1000);
                   chiarr.add_room(1000);
                   //printf("done %d\n",dexes.get_dim());
                   //exit(1);
                }
                
                //if(ct%1000==0)printf("    %d\n",ct);
                
            }//loop over yy
        }//loop overxx
    }
    
    if(dexes.get_dim()!=ct){
        printf("WARNING dexes %d but ct %d\n",dexes.get_dim(),ct);
        exit(1);
    }
    printf("time to sort %d -- %d\n",dexes.get_dim(),ct);
    
    array_1d<double> chisorted;
    sort_and_check(chiarr,chisorted,dexes);
    
    array_2d<double> scatter_data;
    scatter_data.set_cols(2);
    
    double sum=0.0;
    FILE *output;
    
    
    for(i=0;i<dexes.get_dim() && sum<lim*total;i++){
        j=dexes.get_data(i);
        sum+=exp(-0.5*chiarr.get_data(j));
        scatter_data.set(i,0,xarr.get_data(j));
        scatter_data.set(i,1,yarr.get_data(j));
        
    }

    
    //printf("moving on to boundary\n");
    
    array_1d<double> min,max;
    min.set(0,0.0);
    min.set(1,0.0);
    max.set(0,0.1*xwidth);
    max.set(1,0.1*ywidth);
    
    kd_tree scatter_tree(scatter_data,min,max);
    
    array_2d<double> boundary_pts;
    array_1d<int> neigh;
    array_1d<double> ddneigh;
  
    for(i=0;i<scatter_data.get_rows();i++){
        scatter_tree.nn_srch(*scatter_data(i),5,neigh,ddneigh);
        
        if(ddneigh.get_data(4)>1.001){
            boundary_pts.add_row(*scatter_data(i));
        }
    }
    
    kd_tree boundary_tree(boundary_pts,min,max);
    
    array_1d<double> origin,last_pt;
    int last_dex=0;
    
    origin.set(0,boundary_tree.get_pt(0,0));
    origin.set(1,boundary_tree.get_pt(0,1));
    last_pt.set(0,origin.get_data(0));
    last_pt.set(1,origin.get_data(1));
    
    output=fopen(filename,"w");

    while(boundary_tree.get_pts()>1){
        fprintf(output,"%e %e\n",last_pt.get_data(0),last_pt.get_data(1));
        
        boundary_tree.remove(last_dex);
        
        boundary_tree.nn_srch(last_pt,1,neigh,ddneigh);
        last_dex=neigh.get_data(0);
        last_pt.set(0,boundary_tree.get_pt(last_dex,0));
        last_pt.set(1,boundary_tree.get_pt(last_dex,1));
    }
    
    fprintf(output,"%e %e\n",last_pt.get_data(0),last_pt.get_data(1));
    fprintf(output,"%e %e\n",origin.get_data(0),origin.get_data(1));
    
    fclose(output);
    
    _called=0;

}
