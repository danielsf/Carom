#include <vector>
#include <iostream>
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/FunctionMinimum.h"


using namespace ROOT::Minuit2;


class test_fn : public FCNBase{

/*private:
    std::vector<double> theMeasurements;
    std::vector<double> thePositions;
    std::vector<double> theMVariances;
    double theErrorDef;*/

public:

    double operator()(const std::vector<double> &pp) const{
        int i;
        double ans=0.0;
        for(i=0;i<6;i++){
            ans += (pp[i]-i*2.1)*(pp[i]-i*2.1);
        }
        std::cout<<ans<<" "<<pp[3]<<std::endl;
        return ans;
    }

    virtual double Up() const{return 1.0;}
    /*std::vector<double> measurements() const{return theMeasurements;}
    std::vector<double> positions() const {return thePositions;}
    std::vector<double> variances() const{return theMVariances;}
    void setErrorDef(double def){theErrorDef=def;}*/

};

int main(int iargc, char *argv[]){

    MnUserParameters pp;
    int i;
    char word[100];
    for(i=0;i<6;i++){
        sprintf(word,"pp%d",i);
        pp.Add(word,7.0,3.0);
    }
    test_fn my_fn;
    MnMigrad test_migrad(my_fn, pp);
    FunctionMinimum min = test_migrad();
    std::vector<double> pp_out = min.UserParameters().Params();
    for(i=0;i<6;i++){
        std::cout<<pp_out[i]<<std::endl;
    }

}
