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
        void build(array_2d<double>&);

        int get_dim(){return _bases.get_rows();}
        int contains(const array_1d<double>&);
        double bases(int i,int j){return _bases.get_data(i,j);}
        double center(int i){return _center.get_data(i);}
        double radii(int i){return _radii.get_data(i);}

    private:
        array_2d<double> _bases;
        array_1d<double> _radii,_center;

        void _set_radii(array_2d<double>&);

};

#endif
