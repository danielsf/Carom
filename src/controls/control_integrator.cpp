#include "controls/control_integrator.h"

control_integrator::control_integrator(function_wrapper &chisq,
                     array_1d<double> &min, array_1d<double> &max,
                     array_1d<double> &dx, char *name_root){

    _min.set_name("control_min");
    _max.set_name("control_max");
    _dx.set_name("control_dx");

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

    _chi_lim_freq=124.4;
    _chi_min=2.0*exception_value;

    _confidence_limit=0.95;

}

double control_integrator::get_min(int ii){
    return _min.get_data(ii);
}

double control_integrator::get_max(int ii){
    return _max.get_data(ii);
}

void control_integrator::set_chi_lim_freq(double dd){
    _chi_lim_freq=dd;
}

void control_integrator::run_analysis(){
    run_analysis(0.95);
}

void control_integrator::run_analysis(double cc){

    chisquared_distribution distro;

    _confidence_limit=cc;

    _delta_chi_bayes=distro.confidence_limit(double(_min.get_dim()),_confidence_limit);
    _chi_lim_freq=distro.confidence_limit(100.0,_confidence_limit);

    printf("nameroot %s\n",_name_root);

    array_1d<double> pt;

    array_1d<double> chi_vals,good_min,good_max;
    array_2d<int> coordinates;

    good_min.set_name("good_min");
    good_max.set_name("good_max");
    chi_vals.set_name("chi_vals");
    coordinates.set_name("coordinates");
    coordinates.set_dim(1000000,_min.get_dim());
    chi_vals.set_dim(1000000);
    coordinates.reset_preserving_room();
    chi_vals.reset_preserving_room();

    double xx,yy,zz;
    int ix,iy,iz;

    for(ix=0;ix<_min.get_dim();ix++){
        good_min.set(ix,2.0*exception_value);
        good_max.set(ix,-2.0*exception_value);
    }

    double d_threshold=20.0;
    double mu;
    double start=double(time(NULL));



    int ct=0;

    double foundMin=exception_value;
    long int total_pts=1;
    array_1d<int> grid_ct;
    grid_ct.set_name("grid_ct");

    for(ix=0;ix<_min.get_dim();ix++){
        iy=int((_max.get_data(ix)-_min.get_data(ix))/_dx.get_data(ix));
        grid_ct.set(ix,iy);
        total_pts*=iy;
    }
    printf("total_pts %ld\n",total_pts);

    FILE *output;

    output=fopen("output/scratch/control_scatter.sav","w");

    long int ipt;
    array_1d<int> idx;
    idx.set_name("idx");
    for(ipt=0;ipt<total_pts;ipt++){
        expand_grid(ipt,grid_ct,idx);
        for(ix=0;ix<_min.get_dim();ix++){
            pt.set(ix,_min.get_data(ix)+idx.get_data(ix)*_dx.get_data(ix));
        }

        mu=_chisq[0](pt);

        if(mu<foundMin)foundMin=mu;

        if(mu<_chi_lim_freq+d_threshold){
            for(ix=0;ix<_min.get_dim();ix++){
                if(pt.get_data(ix)<good_min.get_data(ix)){
                    good_min.set(ix,pt.get_data(ix));
                }
                if(pt.get_data(ix)>good_max.get_data(ix)){
                    good_max.set(ix,pt.get_data(ix));
                }
            }
            coordinates.add_row(idx);
            chi_vals.add(mu);
        }

        if(mu<119.8){
            for(ix=0;ix<_min.get_dim();ix++){
                fprintf(output,"%e ",pt.get_data(ix));
            }
            fprintf(output,"\n");
        }


        if(ipt%1000000==0 && ipt>0){
            printf("ct %ld good %d in %e %e -- %e\n",
            ipt,chi_vals.get_dim(),double(time(NULL))-start,
            (double(time(NULL))-start)/double(ipt),foundMin);
        }

    }

    fclose(output);
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
    printf("before simplex chi min is %e\n",_chi_min);
    find_chi_min(good_min,good_max);


    for(i=chi_sorted.get_dim()-1;i>=0;i--){
        chi_sorted.subtract_val(i,_chi_min);
        likelihood_sorted.set(i,exp(-0.5*chi_sorted.get_data(i)));
        total+=exp(-0.5*chi_sorted.get_data(i));
    }

    printf("total %e\n",total);

    char name[letters];

    for(ix=0;ix<_min.get_dim();ix++){
        for(iy=ix+1;iy<_min.get_dim();iy++){
            sprintf(name,"%s_%d_%d",_name_root,ix,iy);
            write_output(ix,iy,chi_sorted,likelihood_sorted,dexes,coordinates,total,
                         _chi_lim_freq,name);

        }
    }

    printf("_delta_chi_bayes: %e\n",_delta_chi_bayes);
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
                    array_1d<int> &dexes, array_2d<int> &coordinates,
                    double totalLikelihood,
                    double chiLimFullD,
                    char *outname_root){



    double delta_chi_2d=6.0;

    printf("minchi %e\n",chi_vals_sorted.get_data(0));
    array_2d<double> frequentistFullD,frequentist2D,frequentist2Deff,bayesianFullD,bayesian2D;
    array_2d<double> frequentistFullDrelative;
    frequentistFullD.set_name("freq_FullD");
    frequentistFullDrelative.set_name("freq_FullDrel");
    frequentist2D.set_name("freq_2d");
    frequentist2Deff.set_name("freq_2deff");
    bayesianFullD.set_name("bayes_FullD");
    bayesian2D.set_name("bayes_2d");

    frequentistFullD.set_cols(2);
    frequentistFullDrelative.set_cols(2);
    frequentist2D.set_cols(2);
    frequentist2Deff.set_cols(2);
    bayesianFullD.set_cols(2);
    bayesian2D.set_cols(2);

    int xct,yct;
    int i,j;
    for(i=0;i<coordinates.get_rows();i++){
        if(i==0 || coordinates.get_data(i,xdex)+1>xct){
            xct=coordinates.get_data(i,xdex)+1;
        }

        if(i==0 || coordinates.get_data(i,ydex)+1>yct){
            yct=coordinates.get_data(i,ydex)+1;
        }
    }

    array_2d<double> minChi2Grid,marginalized_likelihood;
    minChi2Grid.set_name("prepare_output_minChi2Grid");
    marginalized_likelihood.set_name("prepare_output_marginalized_likelihood");
    minChi2Grid.set_dim(xct,yct);
    marginalized_likelihood.set_dim(xct,yct);
    for(i=0;i<xct;i++){
        for(j=0;j<yct;j++){
            minChi2Grid.set(i,j,2.0*exception_value);
            marginalized_likelihood.set(i,j,0.0);
        }
    }

    int row,true_dex,ix,iy;
    for(row=chi_vals_sorted.get_dim()-1;row>=0;row--){
        true_dex=dexes.get_data(row);
        ix=coordinates.get_data(true_dex,xdex);
        iy=coordinates.get_data(true_dex,ydex);


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
        trial.set(0,_min.get_data(xdex)+ix*_dx.get_data(xdex));
        for(iy=0;iy<yct;iy++){
            trial.set(1,_min.get_data(ydex)+iy*_dx.get_data(ydex));

            if(minChi2Grid.get_data(ix,iy)+_chi_min<=chiLimFullD){
                frequentistFullD.add_row(trial);

                for(i=0;i<2;i++){
                    if(trial.get_data(i)<found_min.get_data(i))found_min.set(i,trial.get_data(i));
                    if(trial.get_data(i)>found_max.get_data(i))found_max.set(i,trial.get_data(i));
                }

            }

            if(minChi2Grid.get_data(ix,iy)<=_delta_chi_bayes){
                frequentistFullDrelative.add_row(trial);
            }

            if(minChi2Grid.get_data(ix,iy)<=delta_chi_2d){
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
    for(i=sorted.get_dim()-1;i>=0 && sum<_confidence_limit*totalLikelihood;i--){
        sum+=sorted.get_data(i);

        ix=marginalized_likelihood_dexes.get_data(sorted_dexes.get_data(i),0);
        iy=marginalized_likelihood_dexes.get_data(sorted_dexes.get_data(i),1);

        bayesian2D.set(row,0,_min.get_data(xdex)+ix*_dx.get_data(xdex));
        bayesian2D.set(row,1,_min.get_data(ydex)+iy*_dx.get_data(ydex));
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
    for(i=0;i<chi_vals_sorted.get_dim() && sum<_confidence_limit*totalLikelihood;i++){
        sum+=likelihood_sorted.get_data(i);
        true_dex=dexes.get_data(i);

        ix=coordinates.get_data(true_dex,xdex);
        iy=coordinates.get_data(true_dex,ydex);

        if(used.get_data(ix,iy)==0){
            bayesianFullD.set(row,0,_min.get_data(xdex)+ix*_dx.get_data(xdex));
            bayesianFullD.set(row,1,_min.get_data(ydex)+iy*_dx.get_data(ydex));
            used.set(ix,iy,1);
            row++;
        }

    }

    printf("BayesianFullD %d %e %e\n",bayesianFullD.get_rows(),totalLikelihood,sum);

    double effective_chi;
    for(ix=0;ix<xct;ix++){
        trial.set(0,_min.get_data(xdex)+ix*_dx.get_data(xdex));
        for(iy=0;iy<yct;iy++){
            trial.set(1,_min.get_data(xdex)+iy*_dx.get_data(ydex));

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

    //points wwith chisq < a priori chi^2_lim
    convert_to_boundary(frequentistFullD,_dx.get_data(xdex),_dx.get_data(ydex),boundary);

    sprintf(outname,"%s_frequentistFullD.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING could not open frequentistFullD\n");
        exit(1);
    }
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);

    //points with chisq < chisq_min + delta_chisq(dof)
    sprintf(outname,"%s_frequentistFullDrelative_scatter.txt",outname_root);
    output=fopen(outname,"w");
    for(ix=0;ix<frequentistFullDrelative.get_rows();ix++){
        fprintf(output,"%e %e\n",frequentistFullDrelative.get_data(ix,0),frequentistFullDrelative.get_data(ix,1));
    }
    fclose(output);

    convert_to_boundary(frequentistFullDrelative,_dx.get_data(xdex),_dx.get_data(ydex),boundary);

    sprintf(outname,"%s_frequentistFullDrelative.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING could not open freqFullDrel\n");
        exit(1);
    }
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);

    //point with chisq < chisq_min + 6
    convert_to_boundary(frequentist2D,_dx.get_data(xdex),_dx.get_data(ydex),boundary);

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

    //points wiht -2 ln(marginalized_likelihood) < -2.0*ln(max_marginalized_likelihood) + 6
    convert_to_boundary(frequentist2Deff,_dx.get_data(xdex),_dx.get_data(ydex),boundary);

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


    //sum of marginalized likelihood < _confidence_limit * total
    convert_to_boundary(bayesian2D,_dx.get_data(xdex),_dx.get_data(ydex),boundary);

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


    //sum of individual likelihood in full D space < _confidence_limit * total
    convert_to_boundary(bayesianFullD,_dx.get_data(xdex),_dx.get_data(ydex),boundary);

    sprintf(outname,"%s_bayesianFullD.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING could not open bayesFullD\n");
        exit(1);
    }
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);


    printf("2d ct %d\n",xct*yct);
    printf("min found max found\n");
    printf("%e %e %e %e\n",_min.get_data(xdex),found_min.get_data(0),_max.get_data(xdex),found_max.get_data(0));
    printf("%e %e %e %e\n",_min.get_data(ydex),found_min.get_data(1),_max.get_data(ydex),found_max.get_data(1));
    printf("confidence limit %e\n",_confidence_limit);
    printf("delta_chi_bayes %e\n",_delta_chi_bayes);
    printf("chi_lim_freq %e\n",_chi_lim_freq);

}

void control_integrator::find_chi_min(array_1d<double> &min,
                                      array_1d<double> &max){
    array_1d<double> center;
    array_2d<double> seed;
    Ran dice(42);

    seed.set_cols(min.get_dim());

    int ix,iy;
    for(ix=0;ix<_min.get_dim()+1;ix++){
        for(iy=0;iy<_min.get_dim();iy++){
            seed.set(ix,iy,min.get_data(iy)+dice.doub()*(max.get_data(iy)-min.get_data(iy)));
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
