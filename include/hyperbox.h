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

        const int dim(){return _min.get_dim();}

        const int n_pts(){return _pts.get_rows();}

        const double ln_vol(){return _ln_vol;}

        const double min(int ii){return _min.get_data(ii);}

        const double max(int ii){return _max.get_data(ii);}

        const double pts(int ii, int jj){
            return _pts.get_data(ii,jj);
        }

        void build(const array_2d<double> &pts_in,
                   const array_1d<double> &min_in,
                   const array_1d<double> &max_in){

            if(pts_in.get_cols()<max_in.get_dim() ||
               pts_in.get_cols()<min_in.get_dim() ||
               min_in.get_dim()<max_in.get_dim()){

                printf("dimensions don't line up in hyperbox::build\n");
                printf("pts %d max %d min %d\n",
                pts_in.get_cols(),max_in.get_dim(),min_in.get_dim());

                exit(1);
            }

            _pts.reset_preserving_room();
            _min.reset_preserving_room();
            _max.reset_preserving_room();
            int i;
            for(i=0;i<pts_in.get_rows();i++){
                _pts.add_row(pts_in(i));
            }
            for(i=0;i<min_in.get_dim();i++){
                _min.set(i,min_in.get_data(i));
                _max.set(i,max_in.get_data(i));
                if(_min.get_data(i)>_max.get_data(i)){
                    printf("min/max backwards in hyperbox::build\n");
                    printf("%e %e\n",_min.get_data(i),_max.get_data(i));
                    exit(1);
                }
            }

            _ln_vol=0.0;
            for(i=0;i<dim();i++){
                if(_max.get_data(i)-_min.get_data(i)>1.0e-60){
                    _ln_vol+=log(_max.get_data(i)-_min.get_data(i));
                }
                else{
                    _ln_vol=-69.0*dim();
                }
            }
        }

        void add_point(array_1d<double> &pt){
            int i;
            for(i=0;i<dim();i++){
               if(pt.get_data(i)<_min.get_data(i)){
                   printf("cannot add this point; violates min\n");
                   exit(1);
               }
               if(pt.get_data(i)>_max.get_data(i)){
                   printf("cannot add this point; violates max\n");
                   exit(1);
               }
            }
            _pts.add_row(pt);
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
        double _ln_vol;

        double _n_metric(int,double,array_2d<double>&);

};


class hyperbox_list{

    public:

        hyperbox_list(){
            _ct=0;
            _room=0;
            _hyperbox_list=NULL;
        }

        ~hyperbox_list(){
            if(_hyperbox_list != NULL){
                delete [] _hyperbox_list;
            }
        }

        int ct(){return _ct;}

        int room(){return _room;}

        void set_room(int rr){
            if(_room>=rr){
                return;
            }
            int i;
            hyperbox *buffer;
            if(_ct>0){
                if(_hyperbox_list==NULL){
                    printf("_ct %d but _hyperbox_list is NULL\n",_ct);
                    exit(1);
                }
                buffer=new hyperbox[_ct];
                for(i=0;i<_ct;i++){
                    buffer[i].copy(_hyperbox_list[i]);
                }
                delete [] _hyperbox_list;
                _hyperbox_list = new hyperbox[rr];
                for(i=0;i<_ct;i++){
                    _hyperbox_list[i].copy(buffer[i]);
                }
                delete [] buffer;
            }
            else{
                if(_hyperbox_list!=NULL){
                    delete [] _hyperbox_list;
                }
                _hyperbox_list = new hyperbox[rr];
            }
            _room=rr;
        }

        void add(hyperbox &h_in){
            hyperbox *buffer;
            int i;
            int d_room=3;
            if(_room==0){
                if(_hyperbox_list!=NULL){
                    printf("somehow, room is 0 but _hyperbox_list is not NULL\n");
                    exit(1);
                }
                _hyperbox_list = new hyperbox[3];
                _room=3;
            }
            else if(_ct==_room && _room>0){
                buffer = new hyperbox[_room];
                for(i=0;i<_room;i++){
                    buffer[i].copy(_hyperbox_list[i]);
                }
                delete [] _hyperbox_list;
                _hyperbox_list = new hyperbox[_room+d_room];
                for(i=0;i<_room;i++){
                    _hyperbox_list[i].copy(buffer[i]);
                }
                delete [] buffer;
                _room += d_room;
            }
            _hyperbox_list[_ct].copy(h_in);
            _ct++;
        }

        hyperbox* operator()(int ii){
            if(ii>=_ct){
                printf("askint for hyperbox %d; only have %d\n",ii,_ct);
                exit(1);
            }
            return &_hyperbox_list[ii];
        }

        void reset(){
            _ct=0;
        }

    private:
        int _ct;
        int _room;
        hyperbox *_hyperbox_list;
};

#endif
