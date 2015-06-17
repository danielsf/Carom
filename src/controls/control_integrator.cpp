#include "controls/control_integrator.h"

control_integrator::control_integrator(function_wrapper &chisq,
                     array_1d<double> &min, array_1d<double> &max,
                     array_1d<double> &dx, char *name_root){

    _chisq=&chisq;
    int i;
    for(i=0;i<min.get_dim();i++){
        _min.set(i,min.get_data(i));
        _max.set(i,max.get_data(i));
        _dx.set(i,dx.get_data(i));
    }
    
    for(i=0;i<letters-1 && name_root[i]!=0;i++){
        _name_root[i]=name_root[i];
    }
    _name_root[i]=0;
    
    _chiLim3d=124.0;
    _chi_min=2.0*exception_value;

}

void control_integrator::set_chiLim3d(double dd){
    _chiLim3d=dd;
}

void control_integrator::run_analysis(){
    
    printf("nameroot %s\n",_name_root);
    
    find_chi_min();
    
    array_1d<double> pt;
    
    array_1d<double> chi_vals;
    array_2d<double> coordinates;
    
    chi_vals.set_name("chi_vals");
    coordinates.set_name("coordinates");
    coordinates.set_dim(1000000,3);
    chi_vals.set_dim(1000000);
    coordinates.reset_preserving_room();
    chi_vals.reset_preserving_room();

    double xx,yy,zz;
    int ix,iy,iz;
    
    double d_threshold=20.0;
    double mu;
    double start=double(time(NULL));

    

    int ct=0;
    
    double foundMin=exception_value;
    int total_pts=1;
    array_1d<int> grid_ct;
    grid_ct.set_name("grid_ct");
    
    for(ix=0;ix<_min.get_dim();ix++){
        iy=int((_max.get_data(ix)-_min.get_data(ix))/_dx.get_data(ix));
        grid_ct.set(ix,iy);
        total_pts*=iy;
    }
    printf("total_pts %d\n",total_pts);
    
    int ipt;
    array_1d<int> idx;
    idx.set_name("idx");
    for(ipt=0;ipt<total_pts;ipt++){
        expand_grid(ipt,grid_ct,idx);
        for(ix=0;ix<_min.get_dim();ix++){
            pt.set(ix,_min.get_data(ix)+idx.get_data(ix)*_dx.get_data(ix));
        }
    
        mu=_chisq[0](pt);
                
        if(mu<foundMin)foundMin=mu;
                
        if(mu<_chiLim3d+d_threshold){
            coordinates.add_row(pt);
            chi_vals.add(mu);
        }
                
        if(ipt%1000000==0 && ipt>0){
            printf("ct %d good %d in %e -- %e\n",ipt,chi_vals.get_dim(),double(time(NULL))-start,foundMin);
        }
            
    }
    
    printf("foundMin %e\n",foundMin);
    
    start=double(time(NULL));
    array_1d<double> chi_sorted,likelihood_sorted;
    array_1d<int> dexes;
    int i;
    double total=0.0;
    dexes.set_dim(chi_vals.get_dim());
    likelihood_sorted.set_dim(chi_vals.get_dim());
    for(i=0;i<chi_vals.get_dim();i++)dexes.set(i,i);
    sort_and_check(chi_vals, chi_sorted, dexes);
    printf("sorting took %e\n",double(time(NULL))-start);
    
    double trial_chi_min;
    trial_chi_min=chi_sorted.get_data(0);
    if(trial_chi_min<_chi_min){
        _chi_min=trial_chi_min;
    }
    
    
    for(i=chi_sorted.get_dim()-1;i>=0;i--){
        chi_sorted.subtract_val(i,_chi_min);
        likelihood_sorted.set(i,exp(-0.5*chi_sorted.get_data(i)));
        total+=exp(-0.5*chi_sorted.get_data(i));
    }
    
    printf("total %e\n",total);
    
    char name[letters];
    
    for(ix=0;ix<3;ix++){
        for(iy=ix+1;iy<3;iy++){
            sprintf(name,"%s_%d_%d",_name_root,ix,iy);    
            write_output(ix,iy,chi_sorted,likelihood_sorted,dexes,coordinates,total,
                         _min,_max,_dx,_chiLim3d,name);

        }
    }

}

void control_integrator::convert_to_boundary(array_2d<double> &pts, double dx, double dy, int dim, array_2d<double> &output){
    array_1d<double> min,max;
    min.set_name("boundary_min");
    max.set_name("boundary_max");
    int i;
    min.set(0,0.0);
    max.set(0,dx);
    min.set(1,0.0);
    max.set(1,dy);
    array_1d<int> neigh;
    array_1d<double> dd;
    
    neigh.set_name("boundary_neigh");
    dd.set_name("boundary_dd");
    
    array_2d<double> boundary;
    boundary.set_name("boundary_boundary");
    
    //first we need to select all of the points that are on the boundary;
    kd_tree *point_tree;
    
    
    
    printf("making point_tree with %d %d\n",pts.get_rows(),pts.get_cols());
    
    double maxdx=-1.0;
    
    point_tree=new kd_tree(pts,min,max);
    for(i=0;i<pts.get_rows();i++){
        point_tree->nn_srch(i,5,neigh,dd);
        
        if(dd.get_data(4)>maxdx){
            maxdx=dd.get_data(4);
        }
        
        if(dd.get_data(4)>1.001){
            boundary.add_row(pts(i)[0]);
        }
    }
    delete point_tree;
    
    printf("making boundary with %d %d -- %e\n",boundary.get_rows(), boundary.get_cols(),maxdx);
    
    point_tree = new kd_tree(boundary,min,max);
    
    int last_dex;
    array_1d<double> last_point, original_point;
    last_point.set_name("boundary_last_point");
    original_point.set_name("boundary_original_point");
    
    for(i=0;i<dim;i++){
        original_point.set(i,boundary.get_data(0,i));
        last_point.set(i,boundary.get_data(0,i));
    }
    last_dex=0;
    
    output.reset();
    output.set_cols(dim);
    output.add_row(original_point);
    
    while(point_tree->get_pts()>1){
        //printf("removing %d %e %e\n",last_dex,point_tree->get_pt(last_dex,0),point_tree->get_pt(last_dex,1));
        point_tree->remove(last_dex);
        point_tree->nn_srch(last_point,1,neigh,dd);
        output.add_row(point_tree->get_pt(neigh.get_data(0))[0]);
        last_dex=neigh.get_data(0);
        point_tree->get_pt(last_dex, last_point);
    }
    output.add_row(point_tree->get_pt(0)[0]);
    output.add_row(original_point);
    
    delete point_tree;
    
}

int control_integrator::get_dex(double value, double min, double dx, int max){
    int base=int((value-min)/dx);
    
    double error=min+base*dx-value;
    if(error<-0.5*dx)base--;
    else if(error>0.5*dx)base++;
    else if(fabs(error)>dx){
        printf("WARNING in get_dex val %e min %e dx %e base %d\n",
        value,min,dx,base);
        printf("error %e;   %e\n",error,min+base*dx);
        exit(1);
    }
    
    if(base<0){
        return 0;
    }
    else if(base>max){
        return max;
    }
    
    return base;
}


void control_integrator::write_output(int xdex, int ydex,
                    array_1d<double> &chi_vals_sorted,
                    array_1d<double> &likelihood_sorted,
                    array_1d<int> &dexes, array_2d<double> &coordinates,
                    double totalLikelihood,
                    array_1d<double> &min, array_1d<double> &max, array_1d<double> &dx,
                    double chiLim3d,
                    char *outname_root){


    double delta_chi_2d=6.0;
    double delta_chi_3d=7.815;

    printf("minchi %e\n",chi_vals_sorted.get_data(0));
    array_2d<double> frequentist3D,frequentist2D,frequentist2Deff,bayesian3D,bayesian2D;
    array_2d<double> frequentist3Drelative;
    frequentist3D.set_name("freq_3d");
    frequentist3Drelative.set_name("freq_3drel");
    frequentist2D.set_name("freq_2d");
    frequentist2Deff.set_name("freq_2deff");
    bayesian3D.set_name("bayes_3d");
    bayesian2D.set_name("bayes_2d");
    
    frequentist3D.set_cols(2);
    frequentist3Drelative.set_cols(2);
    frequentist2D.set_cols(2);
    frequentist2Deff.set_cols(2);
    bayesian3D.set_cols(2);
    bayesian2D.set_cols(2);

    int xct,yct;    
    xct=int((max.get_data(xdex)-min.get_data(xdex))/dx.get_data(xdex));
    yct=int((max.get_data(ydex)-min.get_data(ydex))/dx.get_data(ydex));
    
    array_2d<double> minChi2Grid,marginalized_likelihood;
    minChi2Grid.set_name("prepare_output_minChi2Grid");
    marginalized_likelihood.set_name("prepare_output_marginalized_likelihood");
    minChi2Grid.set_dim(xct,yct);
    marginalized_likelihood.set_dim(xct,yct);
    int i,j;
    for(i=0;i<xct;i++){
        for(j=0;j<yct;j++){
            minChi2Grid.set(i,j,2.0*exception_value);
            marginalized_likelihood.set(i,j,0.0);
        }
    }

    int row,true_dex,ix,iy;
    for(row=chi_vals_sorted.get_dim()-1;row>=0;row--){
        true_dex=dexes.get_data(row);
        ix=get_dex(coordinates.get_data(true_dex,xdex),min.get_data(xdex),dx.get_data(xdex),xct-1);
        iy=get_dex(coordinates.get_data(true_dex,ydex),min.get_data(ydex),dx.get_data(ydex),yct-1);
        
        
        if(chi_vals_sorted.get_data(row)<minChi2Grid.get_data(ix,iy)){
            minChi2Grid.set(ix,iy,chi_vals_sorted.get_data(row));
        }
        marginalized_likelihood.add_val(ix,iy,likelihood_sorted.get_data(row));
    }
    
    double max_likelihood=-2.0*exception_value;
    for(ix=0;ix<xct;ix++){
        for(iy=0;iy<yct;iy++){
            if(marginalized_likelihood.get_data(ix,iy)>max_likelihood){
                max_likelihood=marginalized_likelihood.get_data(ix,iy);
            }
        }
    }
    
    printf("max_likelihood %e\n",max_likelihood);
    
    double min_chi2_eff=-2.0*log(max_likelihood);
    
    
    
    //////assemble the raw frequentist grids
    array_1d<double> trial,found_min,found_max;
    trial.set_name("trial");
    found_min.set_name("found_min");
    found_max.set_name("found_max");
    
    for(ix=0;ix<2;ix++){
        found_min.set(ix,2.0*exception_value);
        found_max.set(ix,-2.0*exception_value);
    }
    
    for(ix=0;ix<xct;ix++){
        trial.set(0,min.get_data(xdex)+ix*dx.get_data(xdex));
        for(iy=0;iy<yct;iy++){
            trial.set(1,min.get_data(ydex)+iy*dx.get_data(ydex));
            
            if(minChi2Grid.get_data(ix,iy)+_chi_min<=chiLim3d){
                frequentist3D.add_row(trial);
                
                for(i=0;i<2;i++){
                    if(trial.get_data(i)<found_min.get_data(i))found_min.set(i,trial.get_data(i));
                    if(trial.get_data(i)>found_max.get_data(i))found_max.set(i,trial.get_data(i));
                }
                
            }
            
            if(minChi2Grid.get_data(ix,iy)<=chi_vals_sorted.get_data(0)+delta_chi_3d){
                frequentist3Drelative.add_row(trial);
            }
            
            if(minChi2Grid.get_data(ix,iy)<=chi_vals_sorted.get_data(0)+delta_chi_2d){
                frequentist2D.add_row(trial);
            }
        }
    }
    
    ////assemble the raw bayesian grids
    array_1d<double> marginalized_likelihood_line;
    array_2d<double> marginalized_likelihood_dexes;
    marginalized_likelihood_line.set_name("marg_like_line");
    marginalized_likelihood_dexes.set_name("marg_like_dex");
    marginalized_likelihood_dexes.set_cols(2);
    
    row=0;
    for(ix=0;ix<xct;ix++){
        for(iy=0;iy<yct;iy++){
           marginalized_likelihood_line.set(row,marginalized_likelihood.get_data(ix,iy));
           marginalized_likelihood_dexes.set(row,0,ix);
           marginalized_likelihood_dexes.set(row,1,iy);
           row++;
        }
    }
    
    array_1d<double> sorted;
    array_1d<int> sorted_dexes;
    sorted.set_name("output_sorted");
    sorted_dexes.set_name("output_sorted_dexes");
    for(i=0;i<marginalized_likelihood_line.get_dim();i++){
        sorted_dexes.set(i,i);
    }
    
    sort_and_check(marginalized_likelihood_line,sorted,sorted_dexes);
    double sum=0.0;
    row=0;
    for(i=sorted.get_dim()-1;i>=0 && sum<0.95*totalLikelihood;i--){
        sum+=sorted.get_data(i);
        
        ix=marginalized_likelihood_dexes.get_data(sorted_dexes.get_data(i),0);
        iy=marginalized_likelihood_dexes.get_data(sorted_dexes.get_data(i),1);
        
        bayesian2D.set(row,0,min.get_data(xdex)+ix*dx.get_data(xdex));
        bayesian2D.set(row,1,min.get_data(ydex)+iy*dx.get_data(ydex));
        row++;
    }

    printf("Bayesian2d %d %e %e\n",bayesian2D.get_rows(),totalLikelihood,sum);
    
    array_2d<int> used;
    used.set_name("output_used");
    used.set_dim(xct,yct);
    for(ix=0;ix<xct;ix++){
        for(iy=0;iy<yct;iy++){
            used.set(ix,iy,0);
        }
    }
    
    sum=0.0;
    row=0;
    for(i=0;i<chi_vals_sorted.get_dim() && sum<0.95*totalLikelihood;i++){
        sum+=likelihood_sorted.get_data(i);
        true_dex=dexes.get_data(i);
        
        ix=get_dex(coordinates.get_data(true_dex,xdex),min.get_data(xdex),dx.get_data(xdex),xct-1);
        iy=get_dex(coordinates.get_data(true_dex,ydex),min.get_data(ydex),dx.get_data(ydex),yct-1);
        
        if(used.get_data(ix,iy)==0){
            bayesian3D.set(row,0,coordinates.get_data(true_dex,xdex));
            bayesian3D.set(row,1,coordinates.get_data(true_dex,ydex));
            used.set(ix,iy,1);
            row++;
        }
        
    }

    printf("Bayesian3D %d %e %e\n",bayesian3D.get_rows(),totalLikelihood,sum);

    double effective_chi;
    for(ix=0;ix<xct;ix++){
        trial.set(0,min.get_data(xdex)+ix*dx.get_data(xdex));
        for(iy=0;iy<yct;iy++){
            trial.set(1,min.get_data(xdex)+iy*dx.get_data(ydex));
            
            effective_chi=-2.0*log(marginalized_likelihood.get_data(ix,iy))-min_chi2_eff;
            if(effective_chi<0.0){
                printf("WARNING eff %e min %e\n",-2.0*log(marginalized_likelihood.get_data(ix,iy)),min_chi2_eff);
                exit(1);
            }
            if(effective_chi<=delta_chi_2d){
                frequentist2Deff.add_row(trial);
            }
            
            
        }
    }

    printf("freq2deff %d %e\n",frequentist2Deff.get_rows(),min_chi2_eff);

    ///write outputs
    FILE *output;
    
    char outname[2*letters];
    array_2d<double> boundary;
    boundary.set_name("output_boundary");
    
    convert_to_boundary(frequentist3D,dx.get_data(xdex),dx.get_data(ydex),2,boundary);
    
    sprintf(outname,"%s_frequentist3D.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING could not open frequentist3d\n");
        exit(1);
    }
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);

    convert_to_boundary(frequentist3Drelative,dx.get_data(xdex),dx.get_data(ydex),2,boundary);
    
    sprintf(outname,"%s_frequentist3Drelative.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING could not open freq3drel\n");
        exit(1);
    }
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);

    convert_to_boundary(frequentist2D,dx.get_data(xdex),dx.get_data(ydex),2,boundary);
    
    sprintf(outname,"%s_frequentist2D.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING could not open freq2d\n");
        exit(1);
    }
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);

    convert_to_boundary(frequentist2Deff,dx.get_data(xdex),dx.get_data(ydex),2,boundary);
    
    sprintf(outname,"%s_frequentist2Deff.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING coudl not open freq2deff\n");
        exit(1);
    }
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);
    
    
    convert_to_boundary(bayesian2D,dx.get_data(xdex),dx.get_data(ydex),2,boundary);
    
    sprintf(outname,"%s_bayesian2D.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING could not open bayes2d\n");
        exit(1);
    }
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);
    
    
    convert_to_boundary(bayesian3D,dx.get_data(xdex),dx.get_data(ydex),2,boundary);
    
    sprintf(outname,"%s_bayesian3D.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING could not open bayes3d\n");
        exit(1);
    }
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);
    
    
    printf("2d ct %d\n",xct*yct);
    printf("min found max found\n");
    printf("%e %e %e %e\n",min.get_data(xdex),found_min.get_data(0),max.get_data(xdex),found_max.get_data(0));
    printf("%e %e %e %e\n",min.get_data(ydex),found_min.get_data(1),max.get_data(ydex),found_max.get_data(1));

}

void control_integrator::find_chi_min(){
    array_1d<double> center;
    array_2d<double> seed;
    Ran dice(42);
    
    int ix,iy;
    for(ix=0;ix<_min.get_dim()+1;ix++){
        for(iy=0;iy<_min.get_dim();iy++){
            seed.set(ix,iy,_min.get_data(iy)+dice.doub()*(_max.get_data(iy)-_min.get_data(iy)));
        }
    }
    
    
    simplex_minimizer f_min;
    f_min.set_chisquared(_chisq);
    f_min.set_minmax(_min,_max);
    f_min.set_dice(&dice);
    f_min.use_gradient();
    f_min.find_minimum(seed, center);
    
    double trial_min=f_min.get_minimum();
    printf("simplex found %e\n",trial_min);
    if(trial_min<_chi_min){
        _chi_min=trial_min;
    }
}
