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
gg = g++ -Wno-write-strings -O3 $(INCLUDE)
ff = gfortran -O3

#use these compilers if you are going to use omp
#gg = /opt/local/bin/g++-mp-4.8 -Wno-write-strings -O3 -fopenmp -DUSE_OPENMP -g
#ff = /opt/local/bin/gfortran-mp-4.8 -O3 -g

object/containers.o: src/utils/containers.cpp include/containers.h
	$(gg) -c -o object/containers.o src/utils/containers.cpp

object/goto_tools.o: include/goto_tools.h src/utils/goto_tools.cpp object/containers.o
	$(gg) -c -o object/goto_tools.o src/utils/goto_tools.cpp

test_containers: object/containers.o src/tests/test_containers.cpp object/goto_tools.o
	$(gg) -o bin/test_containers src/tests/test_containers.cpp object/containers.o \
        object/goto_tools.o $(LIBRARIES)

object/kd.o: src/utils/kd.cpp include/kd.h object/containers.o object/goto_tools.o
	$(gg) -c -o object/kd.o src/utils/kd.cpp

test_kd: src/tests/test_kd.cpp object/kd.o
	$(gg) -o bin/test_kd src/tests/test_kd.cpp object/containers.o object/goto_tools.o \
	object/kd.o $(LIBRARIES)

object/eigen_wrapper.o: src/utils/eigen_wrapper.cpp include/eigen_wrapper.h object/goto_tools.o
	$(gg) -c -o object/eigen_wrapper.o src/utils/eigen_wrapper.cpp

test_eigen: src/tests/test_eigen.cpp object/eigen_wrapper.o
	$(gg) -o bin/test_eigen src/tests/test_eigen.cpp object/containers.o object/goto_tools.o \
	object/eigen_wrapper.o $(LIBRARIES)

object/chisq.o: src/utils/chisq.cpp include/chisq.h object/goto_tools.o object/kd.o
	$(gg) -c -o object/chisq.o src/utils/chisq.cpp

object/aps_extractor.o: src/analysis/aps_extractor.cpp include/aps_extractor.h object/goto_tools.o
	$(gg) -c -o object/aps_extractor.o src/analysis/aps_extractor.cpp

s_curve_analysis: src/analysis/s_curve_analyzer.cpp object/chisq.o object/aps_extractor.o
	$(gg) -o bin/s_curve_analysis src/analysis/s_curve_analyzer.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/aps_extractor.o object/chisq.o \
	$(LIBRARIES)

object/wrappers.o: src/utils/wrappers.cpp include/wrappers.h object/chisq.o object/kd.o
	$(gg) -c -o object/wrappers.o src/utils/wrappers.cpp

all:
	make test_containers
	make test_kd
	make test_eigen
	make s_curve_analysis

clean:
	rm object/*.o bin/*
