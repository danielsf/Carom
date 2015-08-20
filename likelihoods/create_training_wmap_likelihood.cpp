#include "kd.h"

class elliptical_model{

    public:

        elliptical_model(array_2d<double> *dd, array_1d<double> *ff, double chisq_target){
            _target=chisq_target;
            _data=dd;
            _fn=ff;
            _bases.set_name("bases");
            _bases.set_cols(_data->get_cols());
            _model_coeffs.set_name("coeffs");
            _basis_associates.set_name("associates");
            _basis_ddsq.set_name("ddsq");
            _basis_mm.set_name("mm");
            _basis_vv.set_name("vv");
            _basis_bb.set_name("bb");
            _min.set_name("min");
            _max.set_name("max");
            _dice=new Ran(42);

            _prepare();
        }

     ~elliptical_model(){
          delete _dice;
     }

    void find_bases();

    void output_tree(char *out_name){
        kd_tree tree(_data[0]);
        tree.write_to_file(out_name);
    }

    void output_model(char *out_name){
        FILE *output;
        output=fopen(out_name,"w");

        int i,j;
        fprintf(output,"%d",_data->get_cols());
        for(i=0;i<_data->get_cols();i++){
            fprintf(output,"%e\n",_model_coeffs.get_data(i));
        }
        for(i=0;i<_data->get_cols();i++){
            for(j=0;j<_data->get_cols();j++){
                fprintf(output,"%e ",_bases.get_data(i,j));
            }
            fprintf(output,"\n");
        }
        for(i=0;i<_fn->get_dim();i++){
            fprintf(output,"%e\n",_fn->get_data(i));
        }
        fclose(output);
    }

    private:

        array_2d<double> *_data;
        array_2d<double> _bases;
        array_1d<double> _model_coeffs;
        array_1d<double> _min,_max;
        array_1d<double> *_fn;
        array_1d<int> _basis_associates;

        array_2d<double> _basis_ddsq;
        array_1d<double> _basis_mm,_basis_vv,_basis_bb;

        int _mindex;
        double _target;
        Ran *_dice;

        void _prepare();
        double basis_error(array_2d<double>&, array_1d<double>&);
        void perturb_bases(int,array_1d<double>&,array_2d<double>&);
        void validate_bases(array_2d<double>&, char*);

        double distance(int i1, int i2){
            double ans=0.0;
            int i;
            for(i=0;i<_data->get_cols();i++){
                ans+=power((_data->get_data(i1,i)-_data->get_data(i2,i))/(_max.get_data(i)-_min.get_data(i)),2);
            }
            return sqrt(ans);
        }


        void project_to_bases(array_1d<double> &in, array_1d<double> &out){
             int i;

            if(_bases.get_rows()!=_data->get_cols()){
                for(i=0;i<_data->get_cols();i++){
                     out.set(i,in.get_data(i));
                }
                return;
            }

            int j;
            for(i=0;i<_data->get_cols();i++){
                out.set(i,0.0);
                for(j=0;j<_data->get_cols();j++){
                    out.add_val(i,in.get_data(j)*_bases.get_data(i,j));
                }
            }
        }
};


void elliptical_model::_prepare(){

    _mindex=-1;
    int i,j;
    double min_val;
    for(i=0;i<_fn->get_dim();i++){
        if(_mindex<0 || _fn->get_data(i)<min_val){
            min_val=_fn->get_data(i);
            _mindex=i;
        }
    }
    for(i=0;i<_data->get_cols();i++){
        for(j=0;j<_data->get_cols();j++){
            if(i==j){
                _bases.set(i,j,1.0);
            }
            else{
                _bases.set(i,j,0.0);
            }
        }
    }

    for(i=0;i<_data->get_cols();i++){
        _min.set(i,2.0*exception_value);
        _max.set(i,-2.0*exception_value);
    }

    for(i=0;i<_data->get_rows();i++){
        if(_fn->get_data(i)<=_target){
            for(j=0;j<_data->get_cols();j++){
                if(_data->get_data(i,j)<_min.get_data(j)){
                    _min.set(j,_data->get_data(i,j));
                }
                if(_data->get_data(i,j)>_max.get_data(j)){
                    _max.set(j,_data->get_data(i,j));
                }
            }
        }
    }

    double tol=0.1*(_target-_fn->get_data(_mindex));
    double midpt=0.5*(_target+_fn->get_data(_mindex));
    double dd,ddmin;
    for(i=0;i<_fn->get_dim();i++){
        if(fabs(_fn->get_data(i)-_target)<1.0e-2 || fabs(_fn->get_data(i)-midpt)<tol){
            ddmin=2.0*exception_value;
            for(j=0;j<_basis_associates.get_dim();j++){
                dd=distance(_basis_associates.get_data(j),i);
                if(dd<ddmin){
                    ddmin=dd;
                }
            }
            if(ddmin>0.1){
                _basis_associates.add(i);
            }
        }
    }
    printf("_basis_associates %d\n",_basis_associates.get_dim());
}


double elliptical_model::basis_error(array_2d<double> &trial_bases, array_1d<double> &trial_model){

    if(_basis_associates.get_dim()<=0){
        printf("WARNING cannot calculate basis error there are only %d associates\n",
        _basis_associates.get_dim());

        exit(1);
    }

    trial_model.zero();
    if(_basis_ddsq.get_rows()>_basis_associates.get_dim()){
        _basis_ddsq.reset();
    }

    if(_basis_ddsq.get_cols()!=_data->get_cols()){
        _basis_ddsq.set_cols(_data->get_cols());
    }

    if(_basis_mm.get_dim()!=_data->get_cols()*_data->get_cols()){
        _basis_mm.set_dim(_data->get_cols()*_data->get_cols());
    }

    if(_basis_bb.get_dim()!=_data->get_cols()){
        _basis_bb.set_dim(_data->get_cols());
    }

    if(_basis_vv.get_dim()!=_data->get_cols()){
        _basis_vv.set_dim(_data->get_cols());
    }

    _basis_mm.zero();
    _basis_bb.zero();
    _basis_vv.zero();
    _basis_ddsq.zero();

    int i,j,ix;
    double mu;
    for(ix=0;ix<_basis_associates.get_dim();ix++){
        for(i=0;i<_data->get_cols();i++){
            mu=0.0;
            for(j=0;j<_data->get_cols();j++){
                mu+=(_data->get_data(_basis_associates.get_data(ix),j)-_data->get_data(_mindex,j))*trial_bases.get_data(i,j);
            }
            _basis_ddsq.set(ix,i,mu*mu);
        }
    }

    for(i=0;i<_data->get_cols();i++){
        for(j=0;j<_basis_associates.get_dim();j++){
            _basis_bb.add_val(i,_basis_ddsq.get_data(j,i)*(_fn->get_data(_basis_associates.get_data(j))-_fn->get_data(_mindex)));
        }
    }

    int k;
    for(i=0;i<_data->get_cols();i++){
        for(j=i;j<_data->get_cols();j++){
            ix=i*_data->get_cols()+j;
            for(k=0;k<_basis_associates.get_dim();k++){
                _basis_mm.add_val(ix,_basis_ddsq.get_data(k,i)*_basis_ddsq.get_data(k,j));
            }
            if(j!=i){
                _basis_mm.set(j*_data->get_cols()+i,_basis_mm.get_data(ix));
            }
        }
    }

    try{
        naive_gaussian_solver(_basis_mm,_basis_bb,trial_model,_data->get_cols());
    }
    catch(int iex){
        printf("WARNING basis_error was no good\n");
        return 2.0*exception_value;
    }

    double error=0.0,chi_model;
    for(i=0;i<_basis_associates.get_dim();i++){
        chi_model=_fn->get_data(_mindex);
        for(j=0;j<_data->get_cols();j++){
            chi_model+=trial_model.get_data(j)*_basis_ddsq.get_data(i,j);
        }
        error+=power(_fn->get_data(_basis_associates.get_data(i))-chi_model,2);
    }

    return error/double(_basis_associates.get_dim());

}

void elliptical_model::perturb_bases(int idim, array_1d<double> &dx, array_2d<double> &bases_out){

    int i,j;
    bases_out.set_cols(_bases.get_cols());
    for(i=0;i<_data->get_cols();i++){
        for(j=0;j<_data->get_cols();j++){
            bases_out.set(i,j,_bases.get_data(i,j));
        }
    }

    for(i=0;i<_data->get_cols();i++){
        bases_out.add_val(idim,i,dx.get_data(i));
    }
    bases_out(idim)->normalize();

    int ix,jx;
    double mu;

    for(ix=idim+1;ix!=idim;){
        if(ix>=_data->get_cols()){
            ix=0;
        }

        for(jx=idim;jx!=ix;){
            if(ix>=_data->get_cols()){
                jx=0;
            }

            mu=0.0;
            for(i=0;i<_data->get_cols();i++){
                mu+=bases_out.get_data(ix,i)*bases_out.get_data(jx,i);
            }
            for(i=0;i<_data->get_cols();i++){
                bases_out.subtract_val(ix,i,mu*bases_out.get_data(jx,i));
            }

            if(jx<_data->get_cols()-1)jx++;
            else jx=0;
        }

        bases_out(ix)->normalize();

        if(ix<_data->get_cols()-1)ix++;
        else ix=0;

    }

    validate_bases(bases_out,"node_perturb_bases");

}

void elliptical_model::validate_bases(array_2d<double> &bases, char *whereami){
    int ix,i,jx;
    double mu;
    /////////////////testing
    for(ix=0;ix<_data->get_cols();ix++){
        bases(ix)->normalize();
        mu=0.0;
        for(i=0;i<_data->get_cols();i++){
            mu+=bases.get_data(ix,i)*bases.get_data(ix,i);
        }
        if(fabs(mu-1.0)>1.0e-6){
            printf("WARNING in %s, square norm %e\n",whereami,mu);
            exit(1);
        }

        for(jx=ix+1;jx<_data->get_cols();jx++){
            mu=0.0;
            for(i=0;i<_data->get_cols();i++){
                mu+=bases.get_data(ix,i)*bases.get_data(jx,i);
            }

            if(fabs(mu)>1.0e-6){
                printf("WARNING in %s, dot product %e\n",whereami,mu);
                exit(1);
            }
        }
    }

}

void elliptical_model::find_bases(){

    int i,j;
    int iFound;

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

    stdev=0.1/sqrt(double(_data->get_cols()));
    stdevlim=1.0e-5/sqrt(double(_data->get_cols()));
    max_abort=_data->get_cols()*100;

    error=basis_error(_bases,_model_coeffs);
    error0=error;
    error1=error;
    errorBest=error;
    aborted=0;
    total_aborted=0;
    changed_bases=0;
    ct=0;

    printf("error0 %e %d\n",error0,_basis_associates.get_dim());

    while(stdev>stdevlim && aborted<max_abort){
        ct++;
        idim=-1;
        while(idim>=_data->get_cols() || idim<0){
            idim=_dice->int32()%_data->get_cols();
        }

        for(i=0;i<_data->get_cols();i++){
            dx.set(i,normal_deviate(_dice,0.0,stdev));
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
            for(i=0;i<_data->get_cols();i++){
                _model_coeffs.set(i,trial_model.get_data(i));
                for(j=0;j<_data->get_cols();j++){
                    _bases.set(i,j,trial_bases.get_data(i,j));
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
            printf("    ct %d error %e from %e min %e\n",ct,errorBest,error0,_fn->get_data(_mindex));
        }
    }
    printf("    ct %d error %e from %e min %e\n",ct,errorBest,error0,_fn->get_data(_mindex));

}


int main(int iargc, char *argv[]){

    char in_name[letters];
    in_name[0]=0;

    char out_name[letters];
    out_name[0]=0;

    char tree_name[letters];
    tree_name[0]=0;

    int i,j;
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

                case 't':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        tree_name[j]=argv[i][j];
                    }
                    tree_name[j]=0;
                break;

                case 'o':
                    i++;
                    for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                        out_name[j]=argv[i][j];
                    }
                    out_name[j]=0;
                break;

            }
        }
    }

    if(in_name[0]==0){
        printf("WARNING did not set in_name\n");
        exit(1);
    }

    if(out_name[0]==0){
        printf("WARNING did not set out_name\n");
        exit(1);
    }

    if(tree_name[0]==0){
        printf("WARNING did not set tree_name\n");
        exit(1);
    }


    char word[letters];
    FILE *input;
    array_2d<double> data;
    array_1d<double> chisq;
    int dim=6;
    double mu,target;

    target=1283.1;
    data.set_cols(dim);

    input=fopen(in_name,"r");
    for(i=0;i<dim+5;i++){
        fscanf(input,"%s",word);
    }

    for(i=0;fscanf(input,"%le",&mu)>0;i++){
        data.set(i,0,mu);
        for(j=1;j<dim;j++){
            fscanf(input,"%le",&mu);
            data.set(i,j,mu);
        }
        fscanf(input,"%le",&mu);
        chisq.add(mu);
        for(j=0;j<3;j++){
             fscanf(input,"%le",&mu);
        }
    }
    fclose(input);

    printf("pts %d\n",data.get_rows());
    elliptical_model model(&data,&chisq,target);
    model.find_bases();
    model.output_model(out_name);
    model.output_tree(tree_name);


}
