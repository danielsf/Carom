CAROM_HOME = /Users/noldor/physics/nonparametric_cosmology/code_sfd/GitDirectory/Carom/

LAPACK_LIB =-L/Users/noldor/physics/lapack-3.1.1/ -llapack
BLAS_LIB =-L/Users/noldor/physics/BLAS/ -lblas

ARPACK_LIB = -L/Users/noldor/physics/ARPACK/ -larpack

FORTRAN_LIB = -L/opt/local/lib/ -lf95 -lgfortran

CAMB_LIB = -L/Users/noldor/physics/CAMB_110419/camb/ -lcamb
CAMB_INC = -I/Users/noldor/physics/CAMB_110419/camb/

CFITSIO_LIB = -L/Users/noldor/physics/cfitsio/ -lcfitsio

WMAP_LIB = -L/Users/noldor/physics/WMAP7likelihood/ -lwmap7 -lpthread
WMAP_INC = -I/Users/noldor/physics/WMAP7likelihood/

R_PATH = /Users/noldor/physics/lib/

LIBRARIES = $(LAPACK_LIB) $(BLAS_LIB) $(ARPACK_LIB) $(FORTRAN_LIB)

INCLUDE = -I$(CAROM_HOME)include/

WMAP_LIBRARIES = $(LIBRARIES) $(WMAP_LIB) $(CAMB_LIB) $(CFITSIO_LIB)

WMAP_INCLUDE = $(INCLUDE) $(CAMB_INC) $(WMAP_INC)

#do not use these compilers with omp
gg = g++ -Wno-write-strings -O3
ff = gfortran -O3

#use these compilers if you are going to use omp
#gg = /opt/local/bin/g++-mp-4.8 -Wno-write-strings -O3 -fopenmp -DUSE_OPENMP -g
#ff = /opt/local/bin/gfortran-mp-4.8 -O3 -g

object/containers.o: src/containers.cpp include/containers.h
	$(gg) -c -o object/containers.o src/containers.cpp

object/goto_tools.o: include/goto_tools.h src/goto_tools.cpp object/containers.o
	$(gg) -c -o object/goto_tools.o src/goto_tools.cpp $(INCLUDE)

test_containers: object/containers.o src/test_containers.cpp object/goto_tools.o
	$(gg) -o bin/test_containers src/test_containers.cpp object/containers.o \
        object/goto_tools.o $(INCLUDE) $(LIBRARIES)

all:
	make test_containers
	make test_kd
	make test_eigen
	make test_box
	make ellipse
	make aps_extract
	make s_control
	make s_curve
	make s_curve_analysis
	make s_curve_mcmc_analysis
	make s_curve_multinest_analysis
	make coverage
clean:
	rm *.o test_containers test_kd test_box test_eigen ellipse \
	aps_extract s_curve s_control s_curve_analysis \
	s_curve_mcmc_analysis s_curve_mcmc_analysis coverage
