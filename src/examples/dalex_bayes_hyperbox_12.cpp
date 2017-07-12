#include "hyperbox.h"
#include "kd.h"
#include "exampleLikelihoods.h"

void pixellate(const array_1d<double> &pt,
               const array_1d<double> &dx,
               array_1d<int> &px){

    int i,j;
    for(i=0;i<pt.get_dim();i++){
        j=int(pt.get_data(i)/dx.get_data(i));
        if(pt.get_data(i)-j*dx.get_data(i)>0.5*dx.get_data(i))j++;
        if(j*dx.get_data(i)-pt.get_data(i)>0.5*dx.get_data(i))j--;

        if(fabs(j*dx.get_data(i)-pt.get_data(i))>0.5*dx.get_data(i)){
            printf("pixellate failed\n");
            printf("input %e fit %e\n",pt.get_data(i),j*dx.get_data(i));
            printf("dx %e delta %e\n",dx.get_data(i), fabs(j*dx.get_data(i)-pt.get_data(i)));
            exit(1);
        }
        px.set(i,j);
    }

}

class hyperbox_integrator{

    public:

        hyperbox_list hb_list;

        hyperbox_integrator(){
            _dalex_pts.set_name("integrator_dalex_pts");
            _dalex_chisq.set_name("integrator_dalex_chisq");
            _dx.set_name("integrator_dx");
            _pixel_list.set_name("integrator_pixel_list");
            _pixel_mapping.set_name("integrator_pixel_mapping");
            _local_pixel.set_name("integrator_local_pixel");
        }

        ~hyperbox_integrator(){}

        int get_n_pts(){
            int sum=0;
            int i;
            for(i=0;i<_pixel_mapping.get_rows();i++){
                sum+=_pixel_mapping.get_cols(i);
            }
            if(sum!=_dalex_pts.get_rows()){
                printf("WARNING _dalex_pts %d; _pixel_mapping %d\n",
                _dalex_pts.get_rows(),sum);
                exit(1);
            }
            return _dalex_pts.get_rows();
        }

        double chisq_min(){
            return _chisq_min;
        }

        void _add_pixel(int pt_dex){
            pixellate(_dalex_pts(pt_dex),_dx,_local_pixel);
            int i,j;
            int pixel_dex=-1;
            for(i=0;i<_pixel_list.get_rows();i++){
                pixel_dex=i;
                for(j=0;j<_dalex_pts.get_cols();j++){
                    if(_pixel_list.get_data(i,j)!=_local_pixel.get_data(j)){
                        pixel_dex=-1;
                        break;
                    }
                }
                if(pixel_dex>=0){
                    break;
                }
            }

            int already_contained=0;
            if(pixel_dex<0){
                _pixel_list.add_row(_local_pixel);
                _pixel_mapping.set(_pixel_list.get_rows()-1,0,pt_dex);
                pixel_dex=_pixel_list.get_rows()-1;
            }
            else{
                already_contained=0;
                for(i=0;i<_pixel_mapping.get_cols(i);i++){
                    if(_pixel_mapping.get_data(pixel_dex,i)==pt_dex){
                        already_contained=1;
                        break;
                    }
                }
                if(already_contained==0){
                    _pixel_mapping.add(pixel_dex,pt_dex);
                }
            }

            double xmin,xmax;
            for(i=0;i<_dalex_pts.get_cols();i++){
                xmin=_pixel_list.get_data(pixel_dex,i)*_dx.get_data(i)-0.5*_dx.get_data(i);
                xmax=_pixel_list.get_data(pixel_dex,i)*_dx.get_data(i)+0.5*_dx.get_data(i);
                if(_dalex_pts.get_data(pt_dex,i)<xmin){
                    printf("pixelization failed on min: %e !> %e\n",
                    _dalex_pts.get_data(pt_dex,i),xmin);
                    exit(1);
                }
                if(_dalex_pts.get_data(pt_dex,i)>xmax){
                    printf("pixelization failed on max: %e !< %e\n",
                    _dalex_pts.get_data(pt_dex,i),xmax);
                    exit(1);
                }
            }

        }

        void set_pts(array_2d<double> &pts, array_1d<double> &xx, double delta_chisq){
             _dalex_pts.reset_preserving_room();
             _dalex_chisq.reset_preserving_room();
             _dx.reset_preserving_room();
             _pixel_list.reset_preserving_room();
             _pixel_mapping.reset_preserving_room();
             int i;
             _chisq_min=2.0*exception_value;
             for(i=0;i<pts.get_rows();i++){
                 _dalex_pts.add_row(pts(i));
                 _dalex_chisq.add(xx.get_data(i));
                 if(i==0 || xx.get_data(i)<_chisq_min){
                     _chisq_min=xx.get_data(i);
                 }
             }

             int j;
             array_1d<double> good_xmin, good_xmax;
             good_xmin.set_name("integrator_set_pts_good_xmin");
             good_xmax.set_name("integrator_set_pts_good_xmax");
             for(i=0;i<_dalex_pts.get_rows();i++){
                 if(_dalex_chisq.get_data(i)<=_chisq_min+delta_chisq){
                     for(j=0;j<_dalex_pts.get_cols();j++){
                         if(j>=good_xmin.get_dim() || _dalex_pts.get_data(i,j)<good_xmin.get_data(j)){
                             good_xmin.set(j,_dalex_pts.get_data(i,j));
                         }
                         if(j>=good_xmax.get_dim() || _dalex_pts.get_data(i,j)>good_xmax.get_data(j)){
                             good_xmax.set(j,_dalex_pts.get_data(i,j));
                         }
                     }
                 }
             }

             for(i=0;i<_dalex_pts.get_cols();i++){
                 _dx.set(i,0.1*(good_xmax.get_data(i)-good_xmin.get_data(i)));
             }

             for(i=0;i<_dalex_pts.get_rows();i++){
                 _add_pixel(i);
             }
             get_n_pts();
             printf("done with initialization\n");
        }

        void add_pt(array_1d<double> &pt, double xx, int box_dex){
            _dalex_pts.add_row(pt);
            _dalex_chisq.add(xx);
            if(_dalex_pts.get_rows()!=_dalex_chisq.get_dim()){
                printf("WARNING array points %d chisq values %d\n",
                _dalex_pts.get_rows(),_dalex_chisq.get_dim());
                exit(1);
            }
            _add_pixel(_dalex_pts.get_rows()-1);
            array_1d<double> dummy_pt;
            dummy_pt.set_name("dummy_pt");
            int i;
            for(i=0;i<pt.get_dim();i++){
                dummy_pt.set(i,pt.get_data(i));
            }
            dummy_pt.set(pt.get_dim(),xx);
            hb_list(box_dex)->add_point(dummy_pt);
        }

        void create_hyperboxes(){
            hb_list.reset();
            hb_list.set_room(_dalex_pts.get_rows());
            int n_box_pts=0;
            int pt_dex;
            double t_pre_init=double(time(NULL));
            array_1d<double> pt;
            pt.set_name("integrator_create_hyperboxes_pt");
            array_2d<double> box_pts;
            box_pts.set_name("integrator_create_hyperboxes_box_pts");
            array_1d<double> box_min,box_max;
            box_min.set_name("integrator_create_hyperboxes_box_min");
            box_max.set_name("integrator_create_hyperboxes_box_max");
            pt.reset_preserving_room();
            hyperbox hb;
            int i,j,k;
            int dim;
            dim=_dalex_pts.get_cols();
            for(i=0;i<_pixel_list.get_rows();i++){
                box_pts.reset_preserving_room();
                for(j=0;j<dim;j++){
                    box_min.set(j,_pixel_list.get_data(i,j)*_dx.get_data(j)-0.5*_dx.get_data(j));
                    box_max.set(j,_pixel_list.get_data(i,j)*_dx.get_data(j)+0.5*_dx.get_data(j));
                    if(box_min.get_data(j)>box_max.get_data(j)){
                       printf("setting min/max backwards\n");
                       printf("%e  %e\n",box_min.get_data(j),box_max.get_data(j));
                       printf("%d %e\n",_pixel_list.get_data(i,j),_dx.get_data(j));
                       exit(1);
                    }
                }
                for(j=0;j<_pixel_mapping.get_cols(i);j++){
                    pt_dex = _pixel_mapping.get_data(i,j);
                    for(k=0;k<dim;k++){
                        pt.set(k,_dalex_pts.get_data(pt_dex,k));
                        if(pt.get_data(k)>box_max.get_data(k) || pt.get_data(k)<box_min.get_data(k)){
                            printf("assigning point that is outside of box\n");
                            printf("failed: %e < %e < %e\n",
                            box_min.get_data(k),pt.get_data(k),box_max.get_data(k));
                            exit(1);
                        }
                    }
                    pt.set(dim,_dalex_chisq.get_data(pt_dex)-_chisq_min);
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

            printf("hyperboxes %d; points %d; array_pts %d\n",hb_list.ct(),n_box_pts,get_n_pts());

            for(i=0;i<hb_list.ct();i++){
                if(hb_list(i)->dim()!=dim){
                    printf("WARNING initial box has dim %d\n",
                    hb_list(i)->dim());
                    exit(1);
                 }
            }
        }

        void split_hyperboxes(){
            int i,j;
            hyperbox hb;
            int n_pts=0;
            for(i=0;i<hb_list.ct();i++){
                n_pts+=hb_list(i)->n_pts();
            }
            int room0=hb_list.room();
            if(n_pts>hb_list.room()){
                hb_list.set_room(room0+100000);
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
                if(hb_list(i)->n_pts()!=1){
                    printf("have box with %d pts\n",hb_list(i)->n_pts());
                    exit(1);
                }
            }

            n_pts=0;
            for(i=0;i<hb_list.ct();i++){
                n_pts+=hb_list(i)->n_pts();
            }
            if(n_pts!=_dalex_pts.get_rows()){
                printf("WARNING; pt mismatch\n");
                exit(1);
            }
        }

    private:
        array_2d<double> _dalex_pts;
        array_1d<double> _dalex_chisq,_dx;
        array_2d<int> _pixel_list;
        asymm_array_2d<int> _pixel_mapping;
        double _chisq_min;
        array_1d<int> _local_pixel;

};


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
    array_1d<double> dalex_chisq,xmin,xmax;
    array_2d<double> dalex_pts;
    double xx;
    dalex_chisq.set_name("dalex_chisq");
    dalex_pts.set_name("dalex_pts");
    xmin.set_name("xmin");
    xmax.set_name("xmax");

    printf("n_cols %d\n",n_cols);
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
        for(i=dim+1;i<n_cols;i++){
            fscanf(in_file,"%le",&xx);
        }
    }

    fclose(in_file);
    printf("read in data\n");

    kd_tree dalex_tree(dalex_pts,xmin,xmax);

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

    hyperbox_integrator hb_integrator;
    hb_integrator.set_pts(dalex_pts, dalex_chisq, delta_chisq);
    dalex_pts.reset();
    dalex_chisq.reset();

    t0=double(time(NULL));
    hb_integrator.create_hyperboxes();
    hb_integrator.split_hyperboxes();
    t_build_hyperbox+=double(time(NULL))-t0;

    array_2d<double> expanses;
    array_1d<double> local_expanse,sorted_expanse;
    array_1d<double> expanse_limit;
    array_1d<int> expanse_dex;
    expanses.set_name("expanses");
    local_expanse.set_name("local_expanse");
    sorted_expanse.set_name("sorted_expanse");
    expanse_dex.set_name("expanse_dex");
    expanse_limit.set_name("expanse_limit");
    int i_dim;
    int last_printed;
    int all_done=0;

    while(total_pts_added<n_new_pts){
        ln_posterior.reset_preserving_room();
        ln_vol_arr.reset_preserving_room();
        sorted_ln_posterior.reset_preserving_room();
        ln_posterior_dex.reset_preserving_room();
        sorted_chisq.reset_preserving_room();
        chisq_dex.reset_preserving_room();
        posterior_chisq.reset_preserving_room();

        t0=double(time(NULL));
        hb_integrator.split_hyperboxes();
        t_build_hyperbox+=double(time(NULL))-t0;

        posterior_chisq.reset_preserving_room();
        for(i=0;i<hb_integrator.hb_list.ct();i++){
            posterior_chisq.set(i,hb_integrator.hb_list(i)->pts(0,dim));
        }
        ln_posterior.set_name("ln_posterior");
        ln_vol_arr.set_name("ln_vol_arr");
        for(i=0;i<hb_integrator.hb_list.ct();i++){
            ln_posterior.set(i,-0.5*hb_integrator.hb_list(i)->pts(0,dim)+hb_integrator.hb_list(i)->ln_vol());
            ln_vol_arr.set(i,hb_integrator.hb_list(i)->ln_vol());
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
        expanses.reset_preserving_room();
        local_expanse.reset_preserving_room();
        expanses.set_cols(dim);
        for(i=0;i<posterior_chisq.get_dim() && local_prob<0.95*total_prob;i++){
            dex=chisq_dex.get_data(i);
            if(fabs(posterior_chisq.get_data(dex)-hb_integrator.hb_list(dex)->pts(0,dim))>1.0e-10){
                printf("WARNING; incorreclty mapping chisq\n");
                exit(1);
            }
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
            for(j=0;j<dim;j++){
                local_expanse.set(j,hb_integrator.hb_list(dex)->max(j)-hb_integrator.hb_list(dex)->min(j));
            }
            expanses.add_row(local_expanse);
        }

        expanse_limit.reset_preserving_room();
        for(i_dim=0;i_dim<dim;i_dim++){
            local_expanse.reset_preserving_room();
            sorted_expanse.reset_preserving_room();
            expanse_dex.reset_preserving_room();
            for(i=0;i<expanses.get_rows();i++){
                local_expanse.add(expanses.get_data(i,i_dim));
                expanse_dex.set(i,i);
            }
            sort(local_expanse,sorted_expanse,expanse_dex);
            expanse_limit.set(i_dim,sorted_expanse.get_data(expanse_dex.get_dim()-dim));
        }

        sort(valid_vol,valid_vol_sorted,valid_vol_dex);
        printf("max vol %e -- max chisq %e -- min chisq %e\n",
        valid_vol_sorted.get_data(valid_vol_dex.get_dim()-1),
        max_valid_chisq,hb_integrator.chisq_min());
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
        last_printed=0;
        for(dex=0;dex<hb_integrator.hb_list.ct();dex++){
            if(dex%10000==0 || pts_added-last_printed>10000){
                printf("dex %d of %d pts_added %d\n",dex,hb_integrator.hb_list.ct(),pts_added);
                last_printed=pts_added;
            }
            if(hb_integrator.hb_list(dex)->pts(0,dim)>max_valid_chisq+0.2){
                continue;
            }

            for(i_dim=0;i_dim<dim;i_dim++){

                if(hb_integrator.hb_list(dex)->max(i_dim)-hb_integrator.hb_list(dex)->min(i_dim)<expanse_limit.get_data(i_dim)){
                    continue;
                }

                //printf("acting on %e %e\n",hb_integrator.hb_list(dex)->ln_vol(),hb_integrator.hb_list(dex)->pts(0,dim));
                for(sgn=-1.0;sgn<2.0;sgn+=2.0){
                    for(j=0;j<dim;j++){
                        pt.set(j,0.5*(hb_integrator.hb_list(dex)->max(j)+
                                      hb_integrator.hb_list(dex)->min(j)));
                    }
                    pt.add_val(i_dim,sgn*factor*(hb_integrator.hb_list(dex)->max(i_dim)-
                                             hb_integrator.hb_list(dex)->min(i_dim)));


                    dalex_tree.nn_srch(pt,1,neigh,dist);
                    if(dist.get_data(0)>1.0e-20){
                        xx = chifn[0](pt);
                        hb_integrator.add_pt(pt, xx, dex);
                        dalex_tree.add(pt);
                        pts_added++;
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
        if(dalex_tree.get_pts()!=hb_integrator.get_n_pts()){
            printf("tree and array have different n %d %d\n",
            dalex_tree.get_pts(),hb_integrator.get_n_pts());
            exit(1);
        }
        total_pts_added+=pts_added;
        printf("ran %d in %e; build %e; estimate %e hours\n",
        total_pts_added,double(time(NULL))-t_start,
        t_build_hyperbox,
        n_new_pts*(double(time(NULL))-t_start)/(3600.0*total_pts_added));
        printf("max_valid_chisq %e\n",max_valid_chisq);
    }

    //need to process one last time

        ln_posterior.reset_preserving_room();
        ln_vol_arr.reset_preserving_room();
        sorted_ln_posterior.reset_preserving_room();
        ln_posterior_dex.reset_preserving_room();
        sorted_chisq.reset_preserving_room();
        chisq_dex.reset_preserving_room();
        posterior_chisq.reset_preserving_room();

        t0=double(time(NULL));
        hb_integrator.split_hyperboxes();
        t_build_hyperbox+=double(time(NULL))-t0;

        posterior_chisq.reset_preserving_room();
        for(i=0;i<hb_integrator.hb_list.ct();i++){
            posterior_chisq.set(i,hb_integrator.hb_list(i)->pts(0,dim));
        }
        ln_posterior.set_name("ln_posterior");
        ln_vol_arr.set_name("ln_vol_arr");
        for(i=0;i<hb_integrator.hb_list.ct();i++){
            ln_posterior.set(i,-0.5*hb_integrator.hb_list(i)->pts(0,dim)+hb_integrator.hb_list(i)->ln_vol());
            ln_vol_arr.set(i,hb_integrator.hb_list(i)->ln_vol());
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


    local_prob=0.0;
    FILE *pt_file;
    char pt_name[letters];
    sprintf(pt_name,"%s_pts.sav",out_name);
    pt_file=fopen(pt_name,"w");
    out_file=fopen(out_name,"w");
    for(i=0;i<posterior_chisq.get_dim();i++){
        dex=chisq_dex.get_data(i);
        local_prob+=exp(ln_posterior.get_data(dex));
        fprintf(out_file,"%e %e %e ",
                local_prob/total_prob,
                posterior_chisq.get_data(dex)+hb_integrator.chisq_min(),
                ln_vol_arr.get_data(dex));
        fprintf(pt_file,"%e %e %e ",
                local_prob/total_prob,
                posterior_chisq.get_data(dex)+hb_integrator.chisq_min(),
                ln_vol_arr.get_data(dex));
        for(j=0;j<dim;j++){
            fprintf(out_file,"%e %e ",hb_integrator.hb_list(dex)->min(j),hb_integrator.hb_list(dex)->max(j));
            fprintf(pt_file,"%e ",hb_integrator.hb_list(dex)->pts(0,j));
        }
        fprintf(out_file,"\n");
        fprintf(pt_file,"\n");
    }
    fclose(out_file);

}
