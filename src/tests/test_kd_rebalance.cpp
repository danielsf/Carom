#include "kd.h"

int main(int iargc, char *argv[]){

    Ran chaos(1823);

    char data_name[letters];
    sprintf(data_name,"output/test_data/test_output.txt");
    int dim=33;
    char word[letters];
    FILE *input_file;

    array_2d<double> data;
    data.set_name("data");
    array_1d<double> vector;
    vector.set_name("vector");

    int i,j;
    input_file = fopen(data_name, "r");
    for(i=0;i<dim+3;i++){
        fscanf(input_file,"%s",word);
    }
    printf("last word %s\n",word);

    array_1d<double> min,max;
    min.set_name("min");
    max.set_name("max");

    double mu;
    for(i=0;i<2*dim;i++){
        for(j=0;j<dim;j++){
            fscanf(input_file,"%le",&mu);
            vector.set(j,mu);
            if(j<=min.get_dim() || mu<min.get_data(j)){
                min.set(j,mu);
            }
            if(j<=max.get_dim() || mu>max.get_data(j)){
                max.set(j,mu);
            }
        }
        data.add_row(vector);
        for(j=0;j<2;j++){
            fscanf(input_file,"%le",&mu);
        }
    }

    kd_tree kd_test(data,min,max);
    printf("built tree\n");

    while(fscanf(input_file,"%le",&mu)>0){
        //printf("pts %d -- %e\n",kd_test.get_pts(),mu);
        vector.set(0,mu);
        for(i=1;i<dim;i++){
            fscanf(input_file,"%le",&mu);
            vector.set(i,mu);
        }
        kd_test.add(vector);
        for(i=0;i<2;i++){
           fscanf(input_file,"%le",&mu);
        }
        if(kd_test.get_pts()%10000==0){
            printf("    %d pts loaded\n",kd_test.get_pts());
        }
        //printf("pts %d -- %e\n",kd_test.get_pts(),mu);

    }
    printf("left loop\n");

    printf("pts %d\n",kd_test.get_pts());

    array_1d<double> sorted_vals;
    sorted_vals.set_name("sorted_vals");
    array_1d<int> sorted_dexes;
    sorted_dexes.set_name("sorted_dexes");

    min.reset();
    max.reset();
    vector.reset();
    for(i=0;i<dim;i++){
        sorted_vals.reset_preserving_room();
        sorted_dexes.reset_preserving_room();
        vector.reset_preserving_room();
        for(j=0;j<kd_test.get_pts();j++){
            vector.set(j,kd_test.get_pt(j,i));
            sorted_dexes.set(j,j);
        }
        sort(vector,sorted_vals,sorted_dexes);
        j=sorted_dexes.get_dim();
        min.set(i,sorted_vals.get_data(j/3));
        max.set(i,sorted_vals.get_data(2*j/3));
        printf("    range for %d %e %e\n",
               i,min.get_data(i),max.get_data(i));
    }


    array_2d<double> random_vectors;
    random_vectors.set_name("random_vectors");
    int n_samples = 1000;
    random_vectors.set_dim(n_samples, dim);
    for(i=0;i<n_samples;i++){
        for(j=0;j<dim;j++){
           mu=min.get_data(j)+chaos.doub()*(max.get_data(j)-min.get_data(j));
           random_vectors.set(i,j,mu);
        }
    }

    printf("created random data\n");

    array_1d<int> neigh;
    neigh.set_name("neigh");
    array_1d<double> dd;
    dd.set_name("dd");
    double t_start = double(time(NULL));
    for(i=0;i<n_samples;i++){
        kd_test.nn_srch(random_vectors(i), 1, neigh, dd);
        //printf("    got nn  %d\n",i);
    }
    double t_unbalanced = double(time(NULL))-t_start;
    printf("unabalanced search took %e\n",t_unbalanced);
    t_start = double(time(NULL));
    kd_test.rebalance();
    double t_to_balance = double(time(NULL))-t_start;
    printf("rebalancing took %e\n",t_to_balance);
    t_start = double(time(NULL));
    for(i=0;i<n_samples;i++){
        kd_test.nn_srch(random_vectors(i), 5, neigh, dd);
    }
    double t_balanced = double(time(NULL))-t_start;
    printf("\n\nt_unbalanced %e\nt_to_balance %e\nt_balanced %e\n",
    t_unbalanced, t_to_balance, t_balanced);

}
