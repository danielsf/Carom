#include "containers.h"
#include "goto_tools.h"

#ifndef HYPERBOX_H
#define HYPERBOX_H

class hyperbox{

    public:

        hyperbox(){
            _pts.set_name("hyperbox_pts");
            _max.set_name("max");
            _min.set_name("min");
        }

        ~hyperbox(){}

        void copy(const hyperbox &h1){
            build(h1._pts, h1._min, h1._max);
        }

        const int dim(){return _pts.get_cols();}

        const int n_pts(){return _pts.get_rows();}

        const double min(int ii){return _min.get_data(ii);}

        const double max(int ii){return _max.get_data(ii);}

        void build(const array_2d<double> &pts,
                   const array_1d<double> &max_in,
                   const array_1d<double> &min_in){

            if(pts.get_cols()!=max_in.get_dim() ||
               pts.get_cols()!=min_in.get_dim()){

                printf("dimensions don't line up in hyperbox::build\n");
                printf("pts %d max %d min %d\n",
                pts.get_cols(),max_in.get_dim(),min_in.get_dim());

                exit(1);
            }

            _pts.reset_preserving_room();
            _min.reset_preserving_room();
            _max.reset_preserving_room();
            int i;
            for(i=0;i<pts.get_rows();i++){
                _pts.add_row(pts(i));
            }
            for(i=0;i<dim();i++){
                _min.set(i,min_in.get_data(i));
                _max.set(i,max_in.get_data(i));
                if(_min.get_data(i)>_max.get_data(i)){
                    printf("min/max backwards in hyperbox::build\n");
                    exit(1);
                }
            }
        }

        void _split_on_val(array_2d<double>&,
                           array_1d<double>&,
                           array_1d<double>&,
                           array_2d<double>&,
                           array_1d<double>&,
                           array_1d<double>&,
                           int, double);

        void split(array_2d<double>&,
                   array_1d<double>&,
                   array_1d<double>&,
                   array_2d<double>&,
                   array_1d<double>&,
                   array_1d<double>&);

    private:
        array_2d<double> _pts;
        array_1d<double> _max,_min;

};

#endif
