#include "hyperbox.h"
#include "kd.h"
#include "exampleLikelihoods.h"

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

class hyperbox_integrator{

    public:

        hyperbox_integrator(){
            _dalex_pts.set_name("integrator_dalex_pts");
            _dalex_chisq.set_name("integrator_dalex_chisq");
        }

        ~hyperbox_integrator(){}

        void set_pts(array_2d<double> &pts, array_1d<double> &xx){
             _dalex_pts.reset_preserving_room();
             _dalex_chisq.reset_preserving_room();
             int i;
             for(i=0;i<pts.get_rows();i++){
                 _dalex_pts.add_row(pts(i));
                 _dalex_chisq.add(xx.get_data(i));
             }
        }

        void add_pt(array_1d<double> &pt, double xx){
            _dalex_pts.add_row(pt);
            _dalex_chisq.add(xx);
        }

    private:
        array_2d<double> _dalex_pts;
        array_1d<double> _dalex_chisq;
};

void get_hyperbox_list(array_2d<double> &dalex_pts,
                       array_1d<double> &dalex_chisq,
                       double delta_chisq,
                       hyperbox_list &hb_list){

    double t_start=double(time(NULL));
    printf("starting with %d points\n",dalex_pts.get_rows());
    double pixel_factor=0.1;

    int i,j,k;
    int dim = dalex_pts.get_cols();
    array_1d<double> pt;
    pt.set_name("get_hyperbox_list_pt");

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

    double t_pre_pixel=double(time(NULL));
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
    double t_pixel=double(time(NULL))-t_pre_pixel;

    printf("n pixels %d\n",pixel_list.get_rows());

    hyperbox hb;
    hb_list.set_room(dalex_pts.get_rows());

    array_2d<double> box_pts;
    box_pts.set_name("box_pts");
    array_1d<double> box_min,box_max;
    box_min.set_name("box_min");
    box_max.set_name("box_max");

    int n_box_pts=0;
    int pt_dex;
    double t_pre_init=double(time(NULL));
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
    double t_init = double(time(NULL))-t_pre_init;

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
    double t_pre_split=double(time(NULL));
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
    printf("took %e splitting %e init %e\n",
    double(time(NULL))-t_start,
    double(time(NULL))-t_pre_split,
    t_init);
    printf("pixel %e\n",t_pixel);
    printf("\n");

}

int main(int iargc, char *argv[]){

    int i,j,k,dim;
    char in_name[letters];
    char out_name[letters];
    double delta_chisq=-1.0;
    in_name[0]=0;
    out_name[0]=0;
    dim=-1;
    int n_new_pts=0;
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
                case 'n':
                    i++;
                    n_new_pts=atoi(argv[i]);
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
    double chisq_min=2.0*exception_value;
    array_1d<double> xmin,xmax;
    xmin.set_name("xmin");
    xmax.set_name("xmax");
    while(fscanf(in_file,"%le",&xx)>0){
        pt.set(0,xx);
        if(xmin.get_dim()==0 || xx<xmin.get_data(0)){
            xmin.set(0,xx);
        }
        if(xmax.get_dim()==0 || xx>xmax.get_data(0)){
            xmax.set(0,xx);
        }
        for(i=1;i<dim;i++){
            fscanf(in_file,"%le",&xx);
            pt.set(i,xx);
            if(i>=xmin.get_dim() || xx<xmin.get_data(i)){
                xmin.set(i,xx);
            }
            if(i>=xmax.get_dim() || xx>xmax.get_data(i)){
                xmax.set(i,xx);
            }
        }
        dalex_pts.add_row(pt);
        fscanf(in_file,"%le",&xx);
        dalex_chisq.add(xx);
        if(xx<chisq_min){
            chisq_min=xx;
        }
        for(i=dim+1;i<n_cols;i++){
            fscanf(in_file,"%le",&xx);
        }
    }

    fclose(in_file);
    printf("read in data\n");

    kd_tree dalex_tree(dalex_pts,xmin,xmax);
    hyperbox_list hb_list;

    jellyBeanData *chifn;

    if(dim==12){
        chifn = new gaussianJellyBean12();
    }
    else if(dim==4){
        chifn = new integrableJellyBean();
    }
    else{
        printf("cannot do dim %d\n",dim);
        exit(1);
    }

    double factor;
    int pts_added;
    array_1d<double> dist;
    array_1d<int> neigh;
    dist.set_name("dist");
    neigh.set_name("neigh");
    array_1d<double> ln_posterior;
    array_1d<double> ln_vol_arr;
    double ln_vol;
    array_1d<double> sorted_ln_posterior;
    array_1d<int> ln_posterior_dex;
    double total_prob=0.0;
    array_1d<double> sorted_chisq;
    array_1d<int> chisq_dex;
    array_1d<double> posterior_chisq;
    double local_prob=0.0;
    int dex;

    array_1d<double> valid_vol,valid_vol_sorted;
    array_1d<int> valid_vol_dex;
    valid_vol.set_name("valid_vol");
    valid_vol_sorted.set_name("valid_vol_sorted");
    valid_vol_dex.set_name("valid_vol_dex");

    double vol_max;
    int vol_max_dex;
    double sgn;
    double max_valid_chisq;
    double max_valid_vol;

    FILE *out_file;

    int keep_going;
    int total_pts_added = 0;
    double t_start=double(time(NULL));

    double t_build_hyperbox=0.0;
    double t0;
    while(total_pts_added<n_new_pts){
        ln_posterior.reset_preserving_room();
        ln_vol_arr.reset_preserving_room();
        sorted_ln_posterior.reset_preserving_room();
        ln_posterior_dex.reset_preserving_room();
        sorted_chisq.reset_preserving_room();
        chisq_dex.reset_preserving_room();
        posterior_chisq.reset_preserving_room();
        hb_list.reset();

        t0=double(time(NULL));
        get_hyperbox_list(dalex_pts, dalex_chisq, delta_chisq, hb_list);
        t_build_hyperbox+=double(time(NULL))-t0;

        posterior_chisq.reset_preserving_room();
        for(i=0;i<hb_list.ct();i++){
            posterior_chisq.set(i,hb_list(i)->pts(0,dim));
        }
        ln_posterior.set_name("ln_posterior");
        ln_vol_arr.set_name("ln_vol_arr");
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

        sorted_ln_posterior.set_name("sorted_ln_posterior");
        ln_posterior_dex.set_name("ln_posterior_dex");
        for(i=0;i<ln_posterior.get_dim();i++){
            ln_posterior_dex.set(i,i);
        }
        sort(ln_posterior,sorted_ln_posterior,ln_posterior_dex);
        total_prob=0.0;
        for(i=0;i<ln_posterior.get_dim();i++){
            total_prob+=exp(sorted_ln_posterior.get_data(i));
        }


        sorted_chisq.set_name("sorted_chisq");
        chisq_dex.set_name("chisq_dex");
        for(i=0;i<posterior_chisq.get_dim();i++){
            chisq_dex.set(i,i);
        }

        sort(posterior_chisq,sorted_chisq,chisq_dex);

        //spock
        local_prob=0.0;
        vol_max_dex=-1;
        valid_vol.reset_preserving_room();
        valid_vol_sorted.reset_preserving_room();
        valid_vol_dex.reset_preserving_room();
        max_valid_chisq=-1.0;
        for(i=0;i<posterior_chisq.get_dim() && local_prob<0.95*total_prob;i++){
            dex=chisq_dex.get_data(i);
            if(posterior_chisq.get_data(dex)>max_valid_chisq){
                max_valid_chisq=posterior_chisq.get_data(dex);
            }
            local_prob+=exp(ln_posterior.get_data(dex));
            valid_vol.add(ln_vol_arr.get_data(dex));
            valid_vol_dex.add(dex);
            if(vol_max_dex<0 || ln_vol_arr.get_data(dex)>vol_max){
                vol_max_dex=i;
                vol_max=ln_vol_arr.get_data(dex);
            }
        }

        sort(valid_vol,valid_vol_sorted,valid_vol_dex);
        printf("max vol %e -- max chisq %e\n",
        valid_vol_sorted.get_data(valid_vol_dex.get_dim()-1),
        max_valid_chisq);
        j=0;
        xx=valid_vol_sorted.get_data(valid_vol_dex.get_dim()-1);
        for(i=0;i<valid_vol_sorted.get_dim();i++){
            if(fabs(valid_vol_sorted.get_data(i)-xx)<1.0e-10){
                j++;
            }
        }
        printf("degen %d\n",j);
        max_valid_vol=valid_vol_sorted.get_data(valid_vol_dex.get_dim()-1);

        pt.reset_preserving_room();
        pts_added=0;
        factor=0.25;
        keep_going=1;
        while(pts_added==0){
            for(k=valid_vol_dex.get_dim()-1;k>=0 && (pts_added==0 || keep_going==1);k--){
                dex=valid_vol_dex.get_data(k);
                printf("acting on %e %e\n",ln_vol_arr.get_data(dex),posterior_chisq.get_data(dex));
                for(i=0;i<dim;i++){
                    for(sgn=-1.0;sgn<2.0;sgn+=2.0){
                        for(j=0;j<dim;j++){
                            pt.set(j,0.5*(hb_list(dex)->max(j)+
                                          hb_list(dex)->min(j)));
                        }
                        pt.add_val(i,sgn*factor*(hb_list(dex)->max(i)-
                                                 hb_list(dex)->min(i)));


                        dalex_tree.nn_srch(pt,1,neigh,dist);
                        if(dist.get_data(0)>1.0e-20){
                            dalex_pts.add_row(pt);
                            dalex_tree.add(pt);
                            dalex_chisq.add(chifn[0](pt));
                            pts_added++;
                         }
                    }
                }
                keep_going=0;
                if(k>0){
                    if(fabs(valid_vol_sorted.get_data(k-1)-max_valid_vol)<0.1){
                        keep_going=1;
                    }
                }
            }
        }
        if(pts_added==0){
            printf("did not add any points\n");
            exit(1);
        }
        if(dalex_tree.get_pts()!=dalex_pts.get_rows()){
            printf("tree and array have different n %d %d\n",
            dalex_tree.get_pts(),dalex_pts.get_rows());
            exit(1);
        }
        total_pts_added+=pts_added;
        printf("ran %d in %e; build %e; estimate %e hours\n",
        total_pts_added,double(time(NULL))-t_start,
        t_build_hyperbox,
        n_new_pts*(double(time(NULL))-t_start)/(3600.0*total_pts_added));
    }

    sorted_chisq.set_name("sorted_chisq");
    chisq_dex.set_name("chisq_dex");
    for(i=0;i<posterior_chisq.get_dim();i++){
        chisq_dex.set(i,i);
    }

    sort(posterior_chisq,sorted_chisq,chisq_dex);

    local_prob=0.0;
    out_file=fopen(out_name,"w");
    for(i=0;i<posterior_chisq.get_dim();i++){
        dex=chisq_dex.get_data(i);
        local_prob+=exp(ln_posterior.get_data(dex));
        fprintf(out_file,"%e %e %e ",
                local_prob/total_prob,
                posterior_chisq.get_data(dex)+chisq_min,
                ln_vol_arr.get_data(dex));
        for(j=0;j<dim;j++){
            fprintf(out_file,"%e %e ",hb_list(dex)->min(j),hb_list(dex)->max(j));
        }
        fprintf(out_file,"\n");
    }
    fclose(out_file);

}
