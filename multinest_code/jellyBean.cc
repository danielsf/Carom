#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include "multinest.h"
#include "jellyBean.h"
#include "exampleLikelihoods.h"

double *_global_min,*_global_max;
double _global_chimin;
int _global_called,_global_call_limit;
chisquared *chifn;
//FILE *caromOutput;
FILE *caromMinOutput;

/******************************************** loglikelihood routine ****************************************************/

// Now an example, sample an egg box likelihood

// Input arguments
// ndim 						= dimensionality (total number of free parameters) of the problem
// npars 						= total number of free plus derived parameters
// context						void pointer, any additional information
//
// Input/Output arguments
// Cube[npars] 						= on entry has the ndim parameters in unit-hypercube
//	 						on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//
// Output arguments
// lnew 						= loglikelihood

void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
	double chi = 1.0;
	int i;
        array_1d<double> jellyBeanParams;
        jellyBeanParams.set_name("caller_jellyBeanParams");
	for(i = 0; i < ndim; i++)
	{
		double x = _global_min[i]+(_global_max[i]-_global_min[i])*Cube[i];
		jellyBeanParams.set(i,x);
		Cube[i] = x;
	}

        /*for(i=0;i<ndim;i++){
            fprintf(caromOutput,"%e ",jellyBeanParams.get_data(i));
        }*/

        double chisq=chifn[0](jellyBeanParams);
        if(chisq<_global_chimin){
            _global_chimin=chisq;
            fprintf(caromMinOutput,"%d %e\n",_global_called,_global_chimin);
        }

        //fprintf(caromOutput,"%e 0 0 0\n",chisq);

	lnew = -0.5*chisq;
        _global_called++;
        /*if(_global_called>_global_call_limit){
            printf("aborting min %e calls %d\n",_global_chimin,_global_called);
            exit(1);
        }*/
}

/***********************************************************************************************************************/




/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
//
// Arguments:
//
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value from the default (non-INS) mode
// INSlogZ						= log evidence value from the INS mode
// logZerr						= error on log evidence value
// context						void pointer, any additional information

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context)
{
	// convert the 2D Fortran arrays to C++ arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	
	int i, j;
	
	double postdist[nSamples][nPar + 2];
	for( i = 0; i < nPar + 2; i++ )
		for( j = 0; j < nSamples; j++ )
			postdist[j][i] = posterior[0][i * nSamples + j];
	
	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[nlive][nPar + 1];
	for( i = 0; i < nPar + 1; i++ )
		for( j = 0; j < nlive; j++ )
			pLivePts[j][i] = physLive[0][i * nlive + j];
}

/***********************************************************************************************************************/




/************************************************** Main program *******************************************************/



int main(int argc, char *argv[])
{
	// set the MultiNest sampling parameters
	
	
	int IS = 0;					// do Nested Importance Sampling?
	
	int mmodal = 0;					// do mode separation?
	
	int ceff = 0;					// run in constant efficiency mode?
	
	int nlive = 50;				// number of live points
	
	double efr = 0.8;				// set the required efficiency
	
	double tol = 0.5;				// tol, defines the stopping criteria
	
	int ndims = 4;					// dimensionality (no. of free parameters)

	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	
        double width = 1.0;

        int ii;
        int raw_ndims;
        char root_root[100];
        double delta_space=0.0;
        for(ii=1;ii<argc;ii++){
            if(argv[ii][0]=='-'){
                switch(argv[ii][1]){
                    case 'x':
                        ii++;
                        delta_space=atof(argv[ii]);
                        break;
                    case 'd':
                        ii++;
                        ndims=atoi(argv[ii]);
                        raw_ndims=atoi(argv[ii]);
                        if(ndims==-4){
                            ndims=4;
                            chifn = new integrableJellyBean;
                            sprintf(root_root,"integrableJellyBean");
                        }
                        else if(ndims==4){
                            chifn = new gaussianJellyBean4;
                            sprintf(root_root,"gaussianJellyBean");
                            printf("not going to let you do 4d non integrable\n");
                            exit(1);
                        }
                        else if(ndims==12){
                            chifn = new gaussianJellyBean12;
                            sprintf(root_root,"gaussianJellyBean");
                        }
                        else if(ndims==-12){
                            ndims*=-1;
                            chifn = new nonGaussianLump12;
                            sprintf(root_root,"nonGaussianLump");
                        }
                        else if(ndims==16){
                            chifn = new gaussianJellyBean16;
                            sprintf(root_root,"gaussianJellyBean");
                        }
                        else if(ndims==24){
                            chifn = new gaussianJellyBean24;
                            sprintf(root_root,"gaussianJellyBean");
                        }
                        else if(ndims==8){
                            chifn = new nonGaussianLump8;
                            sprintf(root_root, "nonGaussianLump");
                        }
                        else{
                            printf("Do not know how to handle ndims %d\n",ndims);
                            exit(1);
                        }

                        break;
                    case 't':
                        ii++;
                        tol=atof(argv[ii]);
                        break;
                    case 's':
                        ii++;
                        seed=atoi(argv[ii]);
                        break;
                    case 'n':
                        ii++;
                        nlive=atoi(argv[ii]);
                        break;
                    case 'w':
                        ii++;
                        width=atof(argv[ii]);
                        break;
                    case 'h':
                        printf("d = dim\nt = tol\ns = seed\n");
                        printf("n = nlive\n");
                        exit(1);
                        break;
                }
            }
        }

	int nPar = ndims;					// total no. of parameters including free & derived parameters
	
	int nClsPar = ndims;				// no. of parameters to do mode separation on
	
	int updInt = 1000;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++) pWrap[i] = 0;
	
	char root[100];			// root for output files



	
	int fb = 1;					// need feedback on standard output?
	
	int resume = 1;					// resume from a previous job?
	
	int outfile = 1;				// write output files?
	
	int initMPI = 1;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization
	
	double logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest
	
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	
	void *context = 0;				// not required by MultiNest, any additional information user wants to pass

	
	// calling MultiNest

        _global_chimin=2.0*exception_value;
        _global_called=0;
        _global_call_limit=30000000;
        _global_min = new double[ndims];
        _global_max = new double[ndims];

        for(ii=0;ii<ndims;ii++){
            _global_min[ii]=-40.0;
            _global_max[ii]=40.0;
        }
        //_global_max[3]=80.0;

        if(ndims==12 && raw_ndims>0){
            printf("setting 12d min max\n");
            _global_min[0]=-25.0;
            _global_max[0]=17.0;

            _global_min[1]=-1.0;
            _global_max[1]=72.0;

            _global_min[2]=-20.0;
            _global_max[2]=65.0;

            _global_min[3]=3.0;
            _global_max[3]=54.0;

            _global_min[4]=-23.0;
            _global_max[4]=9.0;

            _global_min[5]=-28.0;
            _global_max[5]=12.0;

            _global_min[6]=-10.0;
            _global_max[6]=46.0;

            _global_min[7]=-27.0;
            _global_max[7]=35.0;

            _global_min[8]=-25.0;
            _global_max[8]=7.0;

            _global_min[9]=-19.0;
            _global_max[9]=23.0;

            _global_min[10]=-10.0;
            _global_max[10]=22.0;

            _global_min[11]=-12.0;
            _global_max[11]=35.0;
            for(ii=0;ii<ndims;ii++){
                _global_max[ii]-=5.0;
                _global_min[ii]+=5.0;
            }
        }
        else if(ndims==8){
            _global_min[7]=-20.0;
            _global_max[7]=50.0;
        }
        else if(ndims==12 && raw_ndims<0){
            _global_min[0]=-45.0;
            _global_max[0]=29.0;

            _global_min[1]=-7.0;
            _global_max[1]=40.0;

            _global_min[2]=-20.0;
            _global_max[2]=24.0;

            _global_min[3]=-15.0;
            _global_max[3]=60.0;

            _global_min[4]=-26.0;
            _global_max[4]=7.0;

            _global_min[5]=-28.0;
            _global_max[5]=17.0;

            _global_min[6]=-32.0;
            _global_max[6]=29.0;

            _global_min[7]=-35.0;
            _global_max[7]=60.0;

            _global_min[8]=-28.0;
            _global_max[8]=7.0;

            _global_min[9]=-40.0;
            _global_max[9]=40.0;

            _global_min[10]=-19.0;
            _global_max[10]=23.0;

            _global_min[11]=-20.0;
            _global_max[11]=35.0;

            for(ii=0;ii<ndims;ii++){
                _global_min[ii]+=5.0;
                _global_max[ii]-=5.0;
            }
        }
        else if(ndims==16 && raw_ndims>0){
            _global_min[0]=-72.0;
            _global_max[0]=-10.0;
            _global_min[1]=-47.0;
            _global_max[1]=11.0;
            _global_min[2]=-54.0;
            _global_max[2]=47.0;
            _global_min[3]=-31.0;
            _global_max[3]=74.0;
            _global_min[4]=-97.0;
            _global_max[4]=4.0;
            _global_min[5]=-27.0;
            _global_max[5]=52.0;
            _global_min[6]=-60.0;
            _global_max[6]=20.0;
            _global_min[7]=-130.0;
            _global_max[7]=0.0;
            _global_min[8]=-47.0;
            _global_max[8]=47.0;
            _global_min[9]=-126.0;
            _global_max[9]=7.0;
            _global_min[10]=-85.0;
            _global_max[10]=23.0;
            _global_min[11]=-57.0;
            _global_max[11]=36.0;
            _global_min[12]=-52.0;
            _global_max[12]=30.0;
            _global_min[13]=0.0;
            _global_max[13]=108.0;
            _global_min[14]=-71.0;
            _global_max[14]=71.0;
            _global_min[15]=-6.0;
            _global_max[15]=61.0;
        }
        /*else if(ndims==4 && raw_ndims<0){
            _global_min[0]=-5.0;
            _global_max[0]=16.0;
            _global_min[1]=-4.0;
            _global_max[1]=17.0;
            _global_min[2]=-15.0;
            _global_max[2]=-1.0;
            _global_min[3]=-20.0;
            _global_max[3]=-1.0;
        }*/

        for(ii=0;ii<ndims;ii++){
            _global_min[ii]-=delta_space;
            _global_max[ii]+=delta_space;
        }

        sprintf(root,"chains/%s_d%d_s%d_n%d_x%.2e_t%.2e",
        root_root,ndims,seed,nlive,delta_space,tol);
        //chifn=new jellyBean(ndims,1.0,20.0);
        //chifn=new jellyBeanData(ndims,1,width,100,0.4,0.4,0.02,20.0);

        /*char caromName[200];
        sprintf(caromName,"chains/%s_d%d_s%d_n%d_t%.2e_carom.sav",root_root,ndims,seed,nlive,tol);*/
        char caromMinName[200];
        sprintf(caromMinName,"chains/%s_d%d_s%d_n%d_x%.2e_t%.2e_carom_min.sav",
        root_root,ndims,seed,nlive,delta_space,tol);

        /*caromOutput=fopen(caromName,"w");
        fprintf(caromOutput,"# ");
        for(ii=0;ii<ndims;ii++){
            fprintf(caromOutput,"p%d ",ii);
        }
        fprintf(caromOutput,"chisq mu sig ling\n");*/

        caromMinOutput=fopen(caromMinName, "w");

	nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,
	logZero, maxiter, LogLike, dumper, context);

        //fclose(caromOutput);
        fclose(caromMinOutput);

        FILE *bounds;
        char bounds_name[500];
        sprintf(bounds_name,"chains/%s_d%d_s%d_n%d_x%.2e_t%.2e_carom_bounds.txt",
        root_root,ndims,seed,nlive,delta_space,tol);
        bounds = fopen(bounds_name, "w");

        fprintf(bounds,"delta_space %e\n",delta_space);
        fprintf(bounds,"total calls %d\n",chifn->get_called());
        fprintf(bounds,"chimin %e\n",_global_chimin);
        fprintf(bounds,"tol was %e\n",tol);
        fprintf(bounds,"chimin %e\n",_global_chimin);
        fprintf(bounds,"ndim %d\n", raw_ndims);
        for(ii=0;ii<ndims;ii++){
            fprintf(bounds,"%d min %e max %e\n",ii,_global_min[ii],_global_max[ii]);
        }
        fclose(bounds);

}

/***********************************************************************************************************************/
