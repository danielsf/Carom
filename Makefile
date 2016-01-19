CAROM_HOME = /Users/danielsf/physics/Carom/

LAPACK_LIB =-L/Users/danielsf/physics/lapack-3.5.0/ -llapack
BLAS_LIB =-L/Users/danielsf/physics/lapack-3.5.0/ -lrefblas

ARPACK_LIB = -L/Users/danielsf/physics/ARPACK/ -larpack_MACOSX

FORTRAN_LIB = -lgfortran

CAMB_LIB = -L/Users/danielsf/physics/CAMB_110419/camb/ -lcamb
CAMB_INCLUDE = -I/Users/danielsf/physics/CAMB_110419/camb/

CFITSIO_LIB = -L/Users/danielsf/physics/cfitsio/ -lcfitsio

WMAP_LIB = -L/Users/danielsf/physics/WMAP7likelihood/ -lwmap7 -lpthread
WMAP_INCLUDE = -I/Users/danielsf/physics/WMAP7likelihood/

R_PATH = /Users/danielsf/physics/lib/

LIBRARIES = $(LAPACK_LIB) $(BLAS_LIB) $(ARPACK_LIB) $(FORTRAN_LIB)

INCLUDE = -I$(CAROM_HOME)include/

WMAP_LIBRARIES = $(WMAP_LIB) $(CAMB_LIB) $(CFITSIO_LIB)

#do not use these compilers with omp
gg = g++ -Wno-write-strings -O3 $(INCLUDE)
ff = gfortran -O3

#use these compilers if you are going to use omp
#gg = /opt/local/bin/g++-mp-4.8 -Wno-write-strings -O3 -fopenmp -DUSE_OPENMP -g
#ff = /opt/local/bin/gfortran-mp-4.8 -O3 -g

object/containers.o: src/utils/containers.cpp include/containers.h
	$(gg) -c -o object/containers.o src/utils/containers.cpp

object/wrappers.o: src/utils/wrappers.cpp include/wrappers.h object/containers.o
	$(gg) -c -o object/wrappers.o src/utils/wrappers.cpp

object/goto_tools.o: include/goto_tools.h src/utils/goto_tools.cpp object/wrappers.o
	$(gg) -c -o object/goto_tools.o src/utils/goto_tools.cpp

test_containers: object/containers.o src/tests/test_containers.cpp object/goto_tools.o
	$(gg) -o bin/test_containers src/tests/test_containers.cpp object/containers.o \
        object/goto_tools.o $(LIBRARIES)

test_grid: object/containers.o src/tests/test_grid.cpp object/goto_tools.o
	$(gg) -o bin/test_grid src/tests/test_grid.cpp object/containers.o \
        object/goto_tools.o $(LIBRARIES)

object/kd.o: src/utils/kd.cpp include/kd.h object/containers.o object/goto_tools.o
	$(gg) -c -o object/kd.o src/utils/kd.cpp

object/kde.o: include/kde.h src/utils/kde.cpp object/goto_tools.o object/kd.o
	$(gg) -c -o object/kde.o src/utils/kde.cpp

test_kd: src/tests/test_kd.cpp object/kd.o
	$(gg) -o bin/test_kd src/tests/test_kd.cpp object/containers.o object/goto_tools.o \
	object/kd.o $(LIBRARIES)

object/eigen_wrapper.o: src/utils/eigen_wrapper.cpp include/eigen_wrapper.h object/goto_tools.o
	$(gg) -c -o object/eigen_wrapper.o src/utils/eigen_wrapper.cpp

test_eigen: src/tests/test_eigen.cpp object/eigen_wrapper.o
	$(gg) -o bin/test_eigen src/tests/test_eigen.cpp object/containers.o object/goto_tools.o \
	object/eigen_wrapper.o $(LIBRARIES)

object/chisq.o: src/utils/chisq.cpp include/chisq.h object/goto_tools.o object/kd.o object/wrappers.o
	$(gg) -c -o object/chisq.o src/utils/chisq.cpp

object/chain.o: src/mcmc/chain.cpp include/mcmc/chain.h object/containers.o \
object/goto_tools.o object/kde.o
	$(gg) -c -o object/chain.o src/mcmc/chain.cpp

object/mcmc.o: src/mcmc/mcmc.cpp include/mcmc/mcmc.h object/chain.o object/chisq.o \
object/simplex.o object/eigen_wrapper.o
	$(gg) -c -o object/mcmc.o src/mcmc/mcmc.cpp

object/camb_wrapper_wmap.o: likelihoods/camb_wrapper_wmap.F90
	$(ff) -c -o object/camb_wrapper_wmap.o likelihoods/camb_wrapper_wmap.F90 $(CAMB_INCLUDE)

object/wmap_wrapper.o: likelihoods/wmap_wrapper.F90
	$(ff) -c -o object/wmap_wrapper.o likelihoods/wmap_wrapper.F90 $(WMAP_INCLUDE)

object/wmap_likelihood_function.o: object/chisq.o object/wmap_wrapper.o object/camb_wrapper_wmap.o \
likelihoods/wmap_likelihood_function.cpp include/wmap_likelihood_function.h
	$(gg) -c -o object/wmap_likelihood_function.o likelihoods/wmap_likelihood_function.cpp \
	$(WMAP_INCLUDE) $(CAMB_INCLUDE)

object/aps_extractor.o: src/analysis/aps_extractor.cpp include/aps_extractor.h object/kd.o object/goto_tools.o
	$(gg) -c -o object/aps_extractor.o src/analysis/aps_extractor.cpp

analysis: src/analysis/generic_analyzer.cpp object/aps_extractor.o
	$(gg) -o bin/analysis src/analysis/generic_analyzer.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/aps_extractor.o \
	$(LIBRARIES)

contours: src/analysis/chisquared_contours.cpp object/kd.o
	$(gg) -o bin/contours src/analysis/chisquared_contours.cpp \
	object/containers.o object/goto_tools.o object/kd.o \
	$(LIBRARIES)

object/chisq_wrapper.o: src/utils/chisq_wrapper.cpp include/chisq_wrapper.h object/wrappers.o object/chisq.o object/kd.o
	$(gg) -c -o object/chisq_wrapper.o src/utils/chisq_wrapper.cpp

object/dchi_simplex.o: src/carom/dchi_simplex.cpp include/dchi_simplex.h object/chisq_wrapper.o
	$(gg) -c -o object/dchi_simplex.o src/carom/dchi_simplex.cpp

object/jellyBean.o: src/utils/jellyBean.cpp include/jellyBean.h object/chisq.o
	$(gg) -c -o object/jellyBean.o src/utils/jellyBean.cpp

object/simplex.o: src/utils/simplex.cpp include/simplex.h object/wrappers.o object/chisq_wrapper.o
	$(gg) -c -o object/simplex.o src/utils/simplex.cpp

object/node.o: src/carom/node.cpp include/node.h object/wrappers.o object/chisq_wrapper.o object/eigen_wrapper.o \
object/simplex.o object/dchi_simplex.o include/search_types.h
	$(gg) -c -o object/node.o src/carom/node.cpp

object/carom.o: src/carom/carom.cpp include/carom.h \
object/simplex.o object/node.o object/eigen_wrapper.o include/search_types.h
	$(gg) -c -o object/carom.o src/carom/carom.cpp

object/control_integrator.o: src/controls/control_integrator.cpp \
include/controls/control_integrator.h object/simplex.o object/kd.o
	$(gg) -c -o object/control_integrator.o src/controls/control_integrator.cpp

jellyBean_control: object/control_integrator.o src/controls/jellyBean_control.cpp \
object/jellyBean.o include/exampleLikelihoods.h
	$(gg) -o bin/jellyBean_control src/controls/jellyBean_control.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/wrappers.o \
	object/chisq.o object/jellyBean.o object/control_integrator.o \
	object/simplex.o \
	$(LIBRARIES)

jellyBean_test: src/examples/jellyBean_example.cpp object/carom.o \
object/jellyBean.o include/exampleLikelihoods.h
	$(gg) -o bin/jellyBean_test src/examples/jellyBean_example.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/wrappers.o object/chisq_wrapper.o object/eigen_wrapper.o object/simplex.o \
	object/dchi_simplex.o object/node.o object/carom.o object/jellyBean.o \
	$(LIBRARIES)

d24_test: src/examples/test_d24_chisq.cpp object/jellyBean.o \
include/exampleLikelihoods.h
	$(gg) -o bin/d24_test src/examples/test_d24_chisq.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/wrappers.o object/chisq_wrapper.o object/eigen_wrapper.o object/simplex.o \
	object/jellyBean.o \
	$(LIBRARIES)

ellipse_test: src/examples/ellipse_example.cpp object/carom.o \
object/jellyBean.o
	$(gg) -o bin/ellipse_test src/examples/ellipse_example.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/wrappers.o object/chisq_wrapper.o object/eigen_wrapper.o object/simplex.o \
	object/node.o object/carom.o object/jellyBean.o \
	$(LIBRARIES)


jellyBean_bayesianControl: src/controls/jellyBeanBayesianControl.cpp object/jellyBean.o
	$(gg) -o bin/jellyBean_bayesianControl src/controls/jellyBeanBayesianControl.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o object/jellyBean.o \
	object/wrappers.o \
	$(LIBRARIES)

jellyBean_frequentistControl: src/controls/jellyBeanFrequentistControl.cpp object/jellyBean.o
	$(gg) -o bin/jellyBean_frequentistControl src/controls/jellyBeanFrequentistControl.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o object/jellyBean.o \
	object/wrappers.o \
	$(LIBRARIES)

jellyBeanMCMC: src/examples/jellyBean_mcmc_example.cpp object/mcmc.o \
object/wmap_likelihood_function.o object/jellyBean.o \
include/exampleLikelihoods.h
	$(gg) -o bin/jellyBeanMCMC src/examples/jellyBean_mcmc_example.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/wrappers.o object/eigen_wrapper.o object/simplex.o \
	object/chain.o object/kde.o object/mcmc.o object/jellyBean.o \
	$(LIBRARIES)

ellipseMCMC: src/examples/ellipse_mcmc_example.cpp object/mcmc.o \
object/wmap_likelihood_function.o object/jellyBean.o
	$(gg) -o bin/ellipseMCMC src/examples/ellipse_mcmc_example.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/wrappers.o object/eigen_wrapper.o object/simplex.o \
	object/chain.o object/kde.o object/mcmc.o object/jellyBean.o \
	$(LIBRARIES)

wmap7: src/examples/wmap7_example.cpp object/carom.o \
object/wmap_likelihood_function.o
	$(gg) -o bin/wmap7 src/examples/wmap7_example.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/camb_wrapper_wmap.o object/wmap_wrapper.o \
	object/wmap_likelihood_function.o \
	object/wrappers.o object/chisq_wrapper.o object/eigen_wrapper.o object/simplex.o \
	object/dchi_simplex.o object/node.o object/carom.o \
	$(WMAP_LIBRARIES) \
	$(WMAP_INCLUDE) $(CAMB_INCLUDE) $(LIBRARIES)

wmap7mcmc: src/examples/wmap7_mcmc_example.cpp object/mcmc.o \
object/wmap_likelihood_function.o
	$(gg) -o bin/wmap7mcmc src/examples/wmap7_mcmc_example.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/camb_wrapper_wmap.o object/wmap_wrapper.o \
	object/wmap_likelihood_function.o \
	object/wrappers.o object/eigen_wrapper.o object/simplex.o \
	object/chain.o object/kde.o object/mcmc.o \
	$(WMAP_LIBRARIES) \
	$(WMAP_INCLUDE) $(CAMB_INCLUDE) $(LIBRARIES)

wmap7reader: src/examples/wmap7_mcmc_reader_example.cpp object/chain.o
	$(gg) -o bin/wmap7reader src/examples/wmap7_mcmc_reader_example.cpp \
	object/containers.o object/goto_tools.o object/chain.o \
	object/kde.o object/kd.o \
	$(LIBRARIES)

MCMCreader: src/examples/generic_mcmc_reader.cpp object/chain.o
	$(gg) -o bin/MCMCreader src/examples/generic_mcmc_reader.cpp \
	object/containers.o object/goto_tools.o object/chain.o \
	object/kde.o object/kd.o \
	$(LIBRARIES)

multiNestReader: src/examples/multiNest_reader.cpp object/chain.o
	$(gg) -o bin/multiNestReader src/examples/multiNest_reader.cpp \
	object/containers.o object/goto_tools.o object/chain.o \
	object/kde.o object/kd.o \
	$(LIBRARIES)

multiNestFullReader: src/examples/multiNest_fullD_reader.cpp object/goto_tools.o
	$(gg) -o bin/multiNestFullReader src/examples/multiNest_fullD_reader.cpp \
	object/containers.o object/goto_tools.o \
	$(LIBRARIES)

wmap7_2d: src/examples/wmap7_2d_example.cpp object/carom.o \
object/wmap_likelihood_function.o
	$(gg) -o bin/wmap7_2d src/examples/wmap7_2d_example.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/camb_wrapper_wmap.o object/wmap_wrapper.o \
	object/wmap_likelihood_function.o \
	object/wrappers.o object/chisq_wrapper.o object/eigen_wrapper.o object/simplex.o \
	object/node.o object/carom.o \
	$(WMAP_LIBRARIES) \
	$(WMAP_INCLUDE) $(CAMB_INCLUDE) $(LIBRARIES)

wmap7_2d_control: src/examples/wmap7_2d_control.cpp \
object/wmap_likelihood_function.o
	$(gg) -o bin/wmap7_2d_control src/examples/wmap7_2d_control.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/camb_wrapper_wmap.o object/wmap_wrapper.o \
	object/wmap_likelihood_function.o \
	$(WMAP_LIBRARIES) \
	$(WMAP_INCLUDE) $(CAMB_INCLUDE) $(LIBRARIES)

all:
	make test_containers
	make test_kd
	make test_eigen

clean:
	rm object/*.o bin/*
