#include "ellipse.h"
#include "exampleLikelihoods.h"

int main(int iarg, char *argv[]){

char in_name[letters];
char out_name[letters];

sprintf(in_name, "output/workspace/jellyBean_d12_s90_output.sav");
sprintf(in_name, "output/workspace/lump_d12_s13_output.sav");
sprintf(out_name, "output/scratch/ellipse_samples.txt");

int dim=12;
int n_cols=16;
double dchisq=21.03;

array_1d<double> pt;
array_2d<double> good_pts;

FILE *input;
input=fopen(in_name, "r");

char word[letters];
int i,j;
double xx;
for(i=0;i<n_cols+1;i++){
    fscanf(input,"%s",word);
}
double chi_min=2.0*exception_value;
while(fscanf(input,"%le",&xx)>0){
    for(i=1;i<dim;i++){
        fscanf(input,"%le",&xx);
    }
    fscanf(input,"%le",&xx);
    if(xx<chi_min){
        chi_min=xx;
    }
    for(i=dim+1;i<n_cols;i++){
        fscanf(input,"%le",&xx);
    }
}
fclose(input);

printf("chi_min %e\n",chi_min);

input=fopen(in_name,"r");
for(i=0;i<n_cols+1;i++){
    fscanf(input,"%s",word);
}
printf("last word %s\n",word);
while(fscanf(input,"%le",&xx)>0){
    pt.set(0,xx);
    for(i=1;i<dim;i++){
        fscanf(input,"%le",&xx);
        pt.set(i,xx);
    }
    fscanf(input,"%le",&xx);
    if(xx<=chi_min+dchisq){
        good_pts.add_row(pt);
    }
    for(i=dim+1;i<n_cols;i++){
        fscanf(input,"%le",&xx);
    }
}

ellipse_list good_ellipse_list;

ellipse *local_ellipse;

local_ellipse = new ellipse;
local_ellipse[0].use_geo_center();
local_ellipse[0].build(good_pts);
local_ellipse[0].careful_set_radii(good_pts);
printf("first ellipse %d;%d %d\n",
local_ellipse[0].dim(),good_pts.get_rows(),good_pts.get_cols());

//local_ellipse[0].multiply_radii(1.05);
good_ellipse_list.add(local_ellipse[0]);

delete local_ellipse;

int ix,iy;
array_2d<double> good_pts_2d;
for(ix=0;ix<dim;ix++){
    for(iy=ix+1;iy<dim;iy++){
        pt.reset();
        good_pts_2d.reset_preserving_room();
        for(i=0;i<good_pts.get_rows();i++){
            pt.set(0,good_pts.get_data(i,ix));
            pt.set(1,good_pts.get_data(i,iy));
            good_pts_2d.add_row(pt);
        }
        local_ellipse = new ellipse;
        local_ellipse[0].use_geo_center();
        local_ellipse[0].build(good_pts_2d);
        local_ellipse[0].careful_set_radii(good_pts_2d);
        //local_ellipse[0].multiply_radii(1.05);
        good_ellipse_list.add(local_ellipse[0]);
        delete local_ellipse;
    }
}

if(good_ellipse_list(0)->dim()!=dim){
    printf("WARNING macro ellipse dim %d\n",good_ellipse_list(0)->dim());
    exit(1);
}

i=0;
for(ix=0;ix<dim;ix++){
    for(iy=ix+1;iy<dim;iy++){
        i++;
        if(good_ellipse_list(1)->dim()!=2){
            printf("WARNING sub ellipse %d %d; %d dim %d\n",
            ix,iy,i,good_ellipse_list(1)->dim());
        }
    }
}

FILE *output;
int i_ell=0;
pt.reset();
double theta,yy,x_true,y_true;
for(ix=0;ix<dim;ix++){
    for(iy=ix+1;iy<dim;iy++){
        i_ell++;
        local_ellipse = good_ellipse_list(i_ell);
        /*xx=0.0;
        for(i=0;i<2;i++){
            xx+=local_ellipse->bases(0,i)*local_ellipse->bases(1,i);
        }
        printf("dot product %d %d %e\n",ix,iy,xx);*/
        sprintf(out_name,"trial_ellipses/ell_%d_%d.txt",ix,iy);
        output=fopen(out_name,"w");
        pt.set(0,local_ellipse->center(0));
        pt.set(1,local_ellipse->center(1));

        pt.add_val(0,local_ellipse->radii(0)*local_ellipse->bases(0,0));
        pt.add_val(1,local_ellipse->radii(0)*local_ellipse->bases(0,1));
        fprintf(output,"%e %e\n",pt.get_data(0),pt.get_data(1));

        pt.set(0,local_ellipse->center(0));
        pt.set(1,local_ellipse->center(1));

        pt.subtract_val(0,local_ellipse->radii(0)*local_ellipse->bases(0,0));
        pt.subtract_val(1,local_ellipse->radii(0)*local_ellipse->bases(0,1));
        fprintf(output,"%e %e\n",pt.get_data(0),pt.get_data(1));

        pt.set(0,local_ellipse->center(0));
        pt.set(1,local_ellipse->center(1));

        pt.add_val(0,local_ellipse->radii(1)*local_ellipse->bases(1,0));
        pt.add_val(1,local_ellipse->radii(1)*local_ellipse->bases(1,1));
        fprintf(output,"%e %e\n",pt.get_data(0),pt.get_data(1));

        pt.set(0,local_ellipse->center(0));
        pt.set(1,local_ellipse->center(1));

        pt.subtract_val(0,local_ellipse->radii(1)*local_ellipse->bases(1,0));
        pt.subtract_val(1,local_ellipse->radii(1)*local_ellipse->bases(1,1));
        fprintf(output,"%e %e\n",pt.get_data(0),pt.get_data(1));


        for(theta=0.0;theta<2.0*pi;theta+=0.05*pi){
            xx = local_ellipse->radii(0)*cos(theta);
            yy = local_ellipse->radii(1)*sin(theta);

            x_true = xx*local_ellipse->bases(0,0)+yy*local_ellipse->bases(1,0);
            y_true = xx*local_ellipse->bases(0,1)+yy*local_ellipse->bases(1,1);
            x_true += local_ellipse->center(0);
            y_true += local_ellipse->center(1);
            fprintf(output,"%e %e\n",x_true,y_true);

        }

        fclose(output);

    }
}

ellipse_sampler sampler;
sampler.initialize(dim,56);

gaussianJellyBean12 chifn;

pt.reset();
array_1d<double> sphere_pt,current_pt;
array_1d<double> pt_2d;

int degen=0;
int accept_it;
int n_samples=0;
double chi_old;
double chi_new;

Ran dice(99);

sprintf(out_name, "output/scratch/ellipse_samples.txt");

ellipse *macro_ellipse;
macro_ellipse = good_ellipse_list(0);

int calls=0;
int duds=0;
printf("looking for first point\n");
accept_it=0;
while(accept_it==0){
    sampler.get_pt(sphere_pt);
    calls++;
    for(i=0;i<dim;i++){
         pt.set(i,macro_ellipse->center(i));
    }
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            pt.add_val(j,sphere_pt.get_data(i)*macro_ellipse->radii(i)*macro_ellipse->bases(i,j));
        }
    }
    accept_it=1;
    i_ell=0;
    for(ix=0;ix<dim && accept_it==1;ix++){
        for(iy=ix+1;iy<dim && accept_it==1;iy++){
            i_ell++;
            local_ellipse=good_ellipse_list(i_ell);
            pt_2d.set(0,pt.get_data(ix));
            pt_2d.set(1,pt.get_data(iy));
            if(local_ellipse->contains(pt_2d)==0){
                accept_it=0;
                duds++;
            }
        }
    }
}

printf("got first point %d %d\n",calls,duds);

chi_old=chifn(pt);
for(i=0;i<dim;i++){
    current_pt.set(i,pt.get_data(i));
}

output=fopen(out_name,"w");

int new_pts=0;
double chi_min_sampled=chi_old;
while(n_samples<5000000){

    sampler.get_pt(sphere_pt);
    calls++;
    for(i=0;i<dim;i++){
         pt.set(i,macro_ellipse->center(i));
    }
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            pt.add_val(j,sphere_pt.get_data(i)*macro_ellipse->radii(i)*macro_ellipse->bases(i,j));
        }
    }
    accept_it=1;
    i_ell=0;
    for(ix=0;ix<dim && accept_it==1;ix++){
        for(iy=ix+1;iy<dim && accept_it==1;iy++){
            i_ell++;
            local_ellipse=good_ellipse_list(i_ell);
            pt_2d.set(0,pt.get_data(ix));
            pt_2d.set(1,pt.get_data(iy));
            if(local_ellipse->contains(pt_2d)==0){
                accept_it=0;
                duds++;
            }
        }
    }
    if(accept_it==1){
        n_samples++;
        chi_new=chifn(pt);
        if(chi_new<chi_min_sampled){
            chi_min_sampled=chi_new;
        }
        accept_it=0;
        if(chi_new<chi_old){
            accept_it=1;
        }
        else{
            if(dice.doub()<exp(-0.5*(chi_new-chi_old))){
                accept_it=1;
            }
        }
        if(accept_it==1){
            fprintf(output,"%e %d",chi_old,degen);
            for(i=0;i<dim;i++){
                fprintf(output," %e",current_pt.get_data(i));
            }
            fprintf(output,"\n");
            degen=1;
            for(i=0;i<dim;i++){
                current_pt.set(i,pt.get_data(i));
            }
            chi_old=chi_new;
            new_pts++;
        }
        else{
            degen++;
        }
    }

    if(calls>0 && calls%10000==0){
        printf("calls %d samples %d duds %d new_pts %d; %.3e\n",
        calls,n_samples,duds,new_pts,chi_min_sampled);
    }

}

fprintf(output,"%e %d",chi_old,degen);
for(i=0;i<dim;i++){
    fprintf(output," %e",current_pt.get_data(i));
}
fprintf(output,"\n");

fclose(output);

}
