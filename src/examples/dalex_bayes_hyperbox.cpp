#include "hyperbox.h"

void pixellate(const array_1d<double> &pt,
               const array_1d<double> &dx,
               const array_1d<double> &min,
               array_1d<int> &px){

    int i,j;
    for(i=0;i<pt.get_dim();i++){
        for(j=0;min.get_data(i)+j*dx.get_data(i)<pt.get_data(i);j++);
        if(min.get_data(i)+j*dx.get_data(i)-pt.get_data(i)>0.5*dx.get_data(i)){
            j--;
        }
        if(j<0){
            printf("WARNING pixellate set j %d\n",j);
            exit(1);
        }
        if(fabs(min.get_data(i)+j*dx.get_data(i)-pt.get_data(i))>0.5*dx.get_data(i)){
            printf("pixellate failed\n");
            exit(1);
        }
        px.set(i,j);
    }

}

int main(int iargc, char *argv[]){

    double pixel_factor=0.1;
    int i,j,k,dim;
    char in_name[letters];
    char out_name[letters];
    double delta_chisq=-1.0;
    in_name[0]=0;
    out_name[0]=0;
    dim=-1;
    for(i=1;i<iargc;i++){
        if(argv[i][0]=='-'){
            switch(argv[i][1]){
                case 'i':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        in_name[j]=argv[i][j];
                    }
                    in_name[j]=0;
                    break;
                case 'o':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        out_name[j]=argv[i][j];
                    }
                    out_name[j]=0;
                    break;
                case 'd':
                    i++;
                    dim=atoi(argv[i]);
                    break;
                case 'c':
                    i++;
                    delta_chisq=atof(argv[i]);
                    break;
            }
        }
    }

    if(dim<0){
        printf("need to specify dim\n");
        exit(1);
    }
    if(in_name[0]==0){
        printf("need to specify in_name\n");
        exit(1);
    }
    if(out_name[0]==0){
        printf("need to specify out_name\n");
        exit(1);
    }
    if(delta_chisq<0.0){
        printf("need to specify delta_chisq\n");
        exit(1);
    }

    int n_cols=0;
    char word[letters];
    word[0]=0;
    FILE *in_file;
    in_file = fopen(in_name, "r");
    while(compare_char("log", word)==0){
        fscanf(in_file,"%s",word);
        if(compare_char("#",word)==0){
            n_cols++;
        }
    }


    array_1d<double> pt;
    pt.set_name("pt");
    array_1d<double> dalex_chisq;
    array_2d<double> dalex_pts;
    double xx;
    dalex_chisq.set_name("dalex_chisq");
    dalex_pts.set_name("dalex_pts");

    printf("n_cols %d\n",n_cols);
    while(fscanf(in_file,"%le",&xx)>0){
        pt.set(0,xx);
        for(i=1;i<dim;i++){
            fscanf(in_file,"%le",&xx);
            pt.set(i,xx);
        }
        dalex_pts.add_row(pt);
        fscanf(in_file,"%le",&xx);
        dalex_chisq.add(xx);
        for(i=dim+1;i<n_cols;i++){
            fscanf(in_file,"%le",&xx);
        }
    }

    fclose(in_file);
    printf("read in data\n");

    double chisq_min=exception_value;
    for(i=0;i<dalex_chisq.get_dim();i++){
        if(dalex_chisq.get_data(i)<chisq_min){
            chisq_min=dalex_chisq.get_data(i);
        }
    }

    array_1d<double> dx;
    dx.set_name("dx");
    array_1d<double> xmin;
    array_1d<double> good_xmin,good_xmax;
    xmin.set_name("xmin");
    good_xmin.set_name("good_xmin");
    good_xmax.set_name("good_xmax");
    for(i=0;i<dalex_pts.get_rows();i++){
        if(dalex_chisq.get_data(i)<chisq_min+delta_chisq){
            for(j=0;j<dim;j++){
                if(j>=good_xmin.get_dim() || dalex_pts.get_data(i,j)<good_xmin.get_data(j)){
                    good_xmin.set(j,dalex_pts.get_data(i,j));
                }
                if(j>=good_xmax.get_dim() || dalex_pts.get_data(i,j)>good_xmax.get_data(j)){
                    good_xmax.set(j,dalex_pts.get_data(i,j));
                }
            }
        }

        for(j=0;j<dim;j++){
            if(j>=xmin.get_dim() || dalex_pts.get_data(i,j)<xmin.get_data(j)){
                xmin.set(j,dalex_pts.get_data(i,j));
            }
        }
    }

    printf("set min max %d %d %e\n",dim,good_xmax.get_dim(),chisq_min);

    for(i=0;i<dim;i++){
        dx.set(i,pixel_factor*(good_xmax.get_data(i)-good_xmin.get_data(i)));
        if(dx.get_data(i)<0.0){
            printf("WARNING dx is %e\n",dx.get_data(i));
            exit(1);
        }
    }

    int ix, iy;

    array_2d<int> pixel_list;
    pixel_list.set_name("pixel_list");
    array_1d<int> pixel;
    pixel.set_name("pixel");
    int is_valid,is_same;
    asymm_array_2d<int> pixel_mapping;
    pixel_mapping.set_name("pixel_mapping");
    for(i=0;i<dalex_pts.get_rows();i++){
        pixellate(dalex_pts(i),dx,xmin,pixel);
        is_valid=1;
        for(j=0;j<pixel_list.get_rows();j++){
            is_same=1;
            for(k=0;k<dim;k++){
                if(pixel_list.get_data(j,k)!=pixel.get_data(k)){
                    is_same=0;
                    break;
                }
            }
            if(is_same==1){
                is_valid=0;
                pixel_mapping.add(j,i);
                break;
            }
        }
        if(is_valid==1){
            pixel_list.add_row(pixel);
            pixel_mapping.set(pixel_list.get_rows()-1,0,i);
        }
    }

    printf("n pixels %d\n",pixel_list.get_rows());
    FILE *out_file;
    out_file=fopen("pixel_test.txt","w");
    for(i=0;i<pixel_list.get_rows();i++){
        for(j=0;j<dim;j++){
            fprintf(out_file,"%le ",xmin.get_data(j)+pixel_list.get_data(i,j)*dx.get_data(j));
        }
        fprintf(out_file,"\n");
    }
    fclose(out_file);

    hyperbox_list hb_list;
    hyperbox hb;
    hb_list.set_room(dalex_pts.get_rows());

    array_2d<double> box_pts;
    box_pts.set_name("box_pts");
    array_1d<double> box_min,box_max;
    box_min.set_name("box_min");
    box_max.set_name("box_max");

    double t_start=double(time(NULL));
    int n_box_pts=0;
    int pt_dex;
    pt.reset_preserving_room();
    for(i=0;i<pixel_list.get_rows();i++){
        box_pts.reset_preserving_room();
        for(j=0;j<dim;j++){
            box_min.set(j,xmin.get_data(j)+pixel_list.get_data(i,j)*dx.get_data(j)-0.5*dx.get_data(j));
            box_max.set(j,xmin.get_data(j)+pixel_list.get_data(i,j)*dx.get_data(j)+0.5*dx.get_data(j));
            if(box_min.get_data(j)>box_max.get_data(j)){
               printf("setting min/max backwards\n");
               printf("%e  %e\n",box_min.get_data(j),box_max.get_data(j));
               printf("%e %d %e\n",xmin.get_data(j),pixel_list.get_data(i,j),dx.get_data(j));
               exit(1);
            }
        }
        for(j=0;j<pixel_mapping.get_cols(i);j++){
            pt_dex = pixel_mapping.get_data(i,j);
            for(k=0;k<dim;k++){
                pt.set(k,dalex_pts.get_data(pt_dex,k));
            }
            pt.set(dim, dalex_chisq.get_data(pt_dex)-chisq_min);
            box_pts.add_row(pt);
        }
        if(box_pts.get_rows()==0){
            printf("trying to initialize empty box\n");
            exit(1);
        }
        for(j=0;j<box_pts.get_rows();j++){
            for(k=0;k<dim;k++){
                if(box_pts.get_data(j,k)<box_min.get_data(k) ||
                   box_pts.get_data(j,k)>box_max.get_data(k)){

                    printf("hyper box assignment failed\n");
                    printf("%e %e %e\n",box_min.get_data(k),
                    box_pts.get_data(j,k),box_max.get_data(k));
                    exit(1);
                }
            }
        }
        n_box_pts+=box_pts.get_rows();
        hb.build(box_pts,box_min,box_max);
        hb_list.add(hb);
    }

    printf("hyperboxes %d; points %d\n",hb_list.ct(),n_box_pts);

    for(i=0;i<hb_list.ct();i++){
        if(hb_list(i)->dim()!=dim){
            printf("WARNING initial box has dim %d\n",
            hb_list(i)->dim());
            exit(1);
        }
    }

    int all_clear=0;
    array_1d<double> min1,min2,max1,max2;
    array_2d<double> pts1,pts2;
    min1.set_name("min1");
    min2.set_name("min2");
    max1.set_name("max1");
    max2.set_name("max2");
    pts1.set_name("pts1");
    pts2.set_name("pts2");
    while(all_clear==0){
        all_clear=1;
        for(i=0;i<hb_list.ct();i++){
            if(hb_list(i)->n_pts()>1){
                all_clear=0;
                hb_list(i)->split(pts1,min1,max1,pts2,min2,max2);
                hb_list(i)->build(pts1,min1,max1);
                hb.build(pts2,min2,max2);
                hb_list.add(hb);
            }
        }
    }

    for(i=0;i<hb_list.ct();i++){
        if(hb_list(i)->n_pts()>1){
            printf("have box with %d pts\n",hb_list(i)->n_pts());
            exit(1);
        }
    }

    printf("%d boxes\n",hb_list.ct());

    dalex_chisq.reset_preserving_room();
    for(i=0;i<hb_list.ct();i++){
        dalex_chisq.set(i,hb_list(i)->pts(0,dim));
    }
    array_1d<double> ln_posterior;
    ln_posterior.set_name("ln_posterior");
    array_1d<double> ln_vol_arr;
    ln_vol_arr.set_name("ln_vol_arr");
    double ln_vol;
    for(i=0;i<hb_list.ct();i++){
        ln_vol=0.0;
        for(j=0;j<dim;j++){
            if(hb_list(i)->max(j)-hb_list(i)->min(j)>1.0e-120){
                ln_vol+=log(hb_list(i)->max(j)-hb_list(i)->min(j));
            }
            else{
                ln_vol=-69.0*dim;
                break;
            }
        }
        ln_posterior.set(i,-0.5*hb_list(i)->pts(0,dim)+ln_vol);
        ln_vol_arr.set(i,ln_vol);
    }

    array_1d<double> sorted_ln_posterior;
    array_1d<int> ln_posterior_dex;
    sorted_ln_posterior.set_name("sorted_ln_posterior");
    ln_posterior_dex.set_name("ln_posterior_dex");
    for(i=0;i<ln_posterior.get_dim();i++){
        ln_posterior_dex.set(i,i);
    }
    sort(ln_posterior,sorted_ln_posterior,ln_posterior_dex);
    double total_prob=0.0;
    for(i=0;i<ln_posterior.get_dim();i++){
        total_prob+=exp(sorted_ln_posterior.get_data(i));
    }

    array_1d<double> sorted_chisq;
    array_1d<int> chisq_dex;
    sorted_chisq.set_name("sorted_chisq");
    chisq_dex.set_name("chisq_dex");
    for(i=0;i<dalex_chisq.get_dim();i++){
        chisq_dex.set(i,i);
    }

    sort(dalex_chisq,sorted_chisq,chisq_dex);

    double local_prob=0.0;
    int dex;

    out_file=fopen(out_name,"w");
    for(i=0;i<dalex_chisq.get_dim();i++){
        dex=chisq_dex.get_data(i);
        local_prob+=exp(ln_posterior.get_data(dex));
        fprintf(out_file,"%e %e %e ",
                local_prob/total_prob,
                dalex_chisq.get_data(dex)+chisq_min,
                ln_vol_arr.get_data(dex));
        for(j=0;j<dim;j++){
            fprintf(out_file,"%e %e ",hb_list(dex)->min(j),hb_list(dex)->max(j));
        }
        fprintf(out_file,"\n");
    }
    fclose(out_file);

}
