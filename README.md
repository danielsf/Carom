# Dalex
Using cost function minimization to directly find confidence limits in
high-dimensional parameter spaces.

## Dependencies

To build and run Dalex, you will need the ARPACK, LAPACK, and BLAS linear
algebra libraries.  We downloaded those libraries from

www.caam.rice.edu/software/ARPACK
www.netlib.org/lapack (for both LAPACK and BLAS)

Once built and installed, you will need to modify the LAPACK_LIB, BLAS_LIB,and
ARPACK_LIB variables in the Makefile to point to the correct libraries.  You may
also need to modify FORTRAN_LIB to allow C++ to link to the libraries (which are
built with Fortran).

## Running examples

The likelihood functions from the paper are made available in easy to build and
run examples.  `Make curved_4d` will produce an executable `bin/curved_4d` that
runs Dalex on the 4-dimensional curved likelihood function from the paper.
`Make curved_12d` produces an example based on the curved 12-dimensional
likelihood function.  `Make ellipse_12d` produces an example based on the
ellipsoidal 12-dimensional likelihood function.  The call signatures for these
functions are all of the form

```
./bin/curved_4d -s S -n N -o out_name -t timing_name
```

where `S` is an integer seeding the random number generator, `N` is the number
of likelihood evaluations to make with Dalex, `out_name` is the name of the file
where the output data will be written and `timing_name` is the name of a file to
which information summarizing the wall-time efficiency of Dalex will be written.
