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
where the output data will be written, and `timing_name` is the name of a file to
which information summarizing the wall-time efficiency of Dalex will be written.

## Output

Running Dalex produces two output files.  The file specified by `out_name` above
contains a list of all of the points samples by Dalex.  The columns in that file
are
```
# p1 p2 ... pD chisq log
```
where `p*` are the parameter space coordinates of the points, `chi^2` is the
value of chi-squared at those points, and `log` is an integer denoting what type
of search was used to sample that particular point (these integers are defined at
the top of `include/chisq_wrapper.h` for the curious; they should not have any
bearing on uses of Dalex; contact the author at scott.f.daniel@gmail.com if you
are interested in using this information, as it is not currently
well-documented).

The file specified by `timing_file` lists information useful for keeping track
of how computationally efficient Dalex is being.  The columns of this file are:

- The number of points Dalex has evaluated and stored in its history.
- The number of times the chi-squared function has been called.  This could be
different from the first column if points have been called at which the
chi-squared function returned an exception but caught it gracefully, in which
case Dalex will not include the point in its history.  Also note thhat this
column will only be non-zero if the chi-squared function's `operator()`
increments the `_called` member variable (see below, where we explain how to
implement your own chi-squared function).
- The current minimum discovered value of chi-squared.
- The current target value of chi-squared on the confidence limit.
- The total wall-clock time elapsed since Dalex began.
- The average wall-clock time per call to chi-squared spent since the last time
timing data was written out.
- The total wall-clock time spent in the chi-squared function.
- The average wall-clock time spent in the chi-squared function per call to
chi-squared over the whole Dalex run.
- The wall-clock time overhead per call to chi-squared all of the chi-squared
evaluations since the last time timing data was written.
