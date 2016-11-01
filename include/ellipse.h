#ifndef ELLIPSE_H
#define ELLIPSE_H

#include "containers.h"
#include "goto_tools.h"

class ellipse{

    public:

        ellipse(){
            _bases.set_name("ellipse_bases");
            _radii.set_name("ellipse_radii");
            _center.set_name("ellipse_center");
        }

        ~ellipse(){}
        void build(const array_2d<double>&);

        int get_dim(){return _bases.get_rows();}
        int contains(const array_1d<double>&);
        double bases(int i,int j){return _bases.get_data(i,j);}
        double center(int i){return _center.get_data(i);}
        double radii(int i){return _radii.get_data(i);}
        int dim(){return _center.get_dim();}
        void copy(ellipse&);

    private:
        array_2d<double> _bases;
        array_1d<double> _radii,_center;

        void _set_radii(const array_2d<double>&);
};


class ellipse_list{
    public:
        ellipse_list(){
            _ellipse_list=NULL;
            _ct=0;
        }

        ~ellipse_list(){
            if(_ellipse_list!=NULL){
                delete [] _ellipse_list;
            }
        }

        void add(ellipse&);
        int ct(){return _ct;}
        void reset(){
            _ct=0;
            delete [] _ellipse_list;
            _ellipse_list=NULL;
        }

        ellipse* operator()(int ii){
            if(ii<0 || ii>=_ct){
                printf("No such ellipse: %d\n",ii);
                exit(1);
            }
            return &_ellipse_list[ii];
        }

    private:
        int _ct;
        ellipse *_ellipse_list;
};

#endif
