This directory contains the source code necessary to get MultiNest version 3.9
to run with the likelihood functions defined used to test Dalex.  To run it, you
will need to copy the Makefile and jellyBean.cc* into a directory in the
MultiNest root directory, and then modify the main MultiNest Makefile to be
aware of the new directory (you should be able to intuit how to do this from the
MultiNest Makefile as downloaded).  The call signature for the executable will
be something like

./jellyBean -d D -n N -s S -x X

where

D = the dimensionality of the parameter space.  Different values of D correspond
to different likelihood functions

    D = 4 will give the 4-dimensional curved likelihood from the Dalex paper. 
    In this case, the output files will have names like
    chains/integrableJellyBean_d%d_s%d_n%d_t1.00e-03.txt
    
    D = 12 will give the 12-dimensional curved likelihood from the Dalex paper.
    In this case, the output files will have names like
    chains/gaussianJellyBean_d%d_s%d_n%d_t1.00e-03.txt
    
    D = -12 will give the 12-dimensional ellipsoidal likelihood from the Dalex
    paper.  In this case, the outptu files will have names like
    chains/nonGaussianLump_d%d_s%s_n%d_t1.00e-03.txt

N = the number of live points

S = the seed for the random number generator

X = a buffer added to the parameter space bounds.  MultiNest requires that you
specify hard bounds in parameter space.  Each of the pre-defined likelihood
functions above come with bounds that should contain the full 95% confidence
limit.  Specifying an X will expand these bounds to go from min-X to max+X, in
case you want to experiment with how expanding the range of parameter space
affects convergence.
*the file is called 'jellyBean.cc' because we were trying to make our curved
likelihood function take on the shape of a jelly bean.  It is up to the reader
whether or not we succeeded.
