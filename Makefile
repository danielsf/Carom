DALEX_HOME = /Users/danielsf/physics/Carom/

LAPACK_LIB =-L/Users/danielsf/physics/lapack-3.5.0/ -llapack
BLAS_LIB =-L/Users/danielsf/physics/lapack-3.5.0/ -lrefblas

ARPACK_LIB = -L/Users/danielsf/physics/ARPACK/ -larpack_MACOSX

FORTRAN_LIB = -lgfortran

R_PATH = /Users/danielsf/physics/lib/

LIBRARIES = $(LAPACK_LIB) $(BLAS_LIB) $(ARPACK_LIB) $(FORTRAN_LIB)

INCLUDE = -I$(DALEX_HOME)include/

# compilers
cc = g++ -Wno-write-strings -O3 $(INCLUDE)

object/containers.o: src/utils/containers.cpp include/containers.h
	$(cc) -c -o object/containers.o src/utils/containers.cpp

object/wrappers.o: src/utils/wrappers.cpp include/wrappers.h object/containers.o
	$(cc) -c -o object/wrappers.o src/utils/wrappers.cpp

object/goto_tools.o: include/goto_tools.h src/utils/goto_tools.cpp object/wrappers.o
	$(cc) -c -o object/goto_tools.o src/utils/goto_tools.cpp

object/hyperbox.o: include/hyperbox.h src/utils/hyperbox.cpp object/goto_tools.o \
object/containers.o
	$(cc) -c -o object/hyperbox.o src/utils/hyperbox.cpp

test_containers: object/containers.o src/tests/test_containers.cpp object/goto_tools.o
	$(cc) -o bin/test_containers src/tests/test_containers.cpp object/containers.o \
        object/goto_tools.o $(LIBRARIES)

object/ellipse.o: object/goto_tools.o object/containers.o \
include/ellipse.h src/utils/ellipse.cpp
	$(cc) -c -o object/ellipse.o src/utils/ellipse.cpp

test_ellipse_creator: object/ellipse.o src/tests/test_ellipse_creator.cpp
	$(cc) -o bin/test_ellipse_creator src/tests/test_ellipse_creator.cpp \
	object/containers.o object/goto_tools.o object/ellipse.o

object/kd.o: src/utils/kd.cpp include/kd.h object/containers.o object/goto_tools.o
	$(cc) -c -o object/kd.o src/utils/kd.cpp

test_kd: src/tests/test_kd.cpp object/kd.o
	$(cc) -o bin/test_kd src/tests/test_kd.cpp object/containers.o object/goto_tools.o \
	object/kd.o $(LIBRARIES)

object/eigen_wrapper.o: src/utils/eigen_wrapper.cpp include/eigen_wrapper.h object/goto_tools.o
	$(cc) -c -o object/eigen_wrapper.o src/utils/eigen_wrapper.cpp

test_eigen: src/tests/test_eigen.cpp object/eigen_wrapper.o
	$(cc) -o bin/test_eigen src/tests/test_eigen.cpp object/containers.o object/goto_tools.o \
	object/eigen_wrapper.o $(LIBRARIES)

object/chisq.o: src/utils/chisq.cpp include/chisq.h object/goto_tools.o object/kd.o object/wrappers.o
	$(cc) -c -o object/chisq.o src/utils/chisq.cpp

object/aps_extractor.o: src/analysis/aps_extractor.cpp include/aps_extractor.h object/kd.o object/goto_tools.o
	$(cc) -c -o object/aps_extractor.o src/analysis/aps_extractor.cpp

analysis: src/analysis/generic_analyzer.cpp object/aps_extractor.o
	$(cc) -o bin/analysis src/analysis/generic_analyzer.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/aps_extractor.o \
	$(LIBRARIES)

contours: src/analysis/chisquared_contours.cpp object/kd.o
	$(cc) -o bin/contours src/analysis/chisquared_contours.cpp \
	object/containers.o object/goto_tools.o object/kd.o \
	$(LIBRARIES)

object/chisq_wrapper.o: src/utils/chisq_wrapper.cpp include/chisq_wrapper.h object/wrappers.o object/chisq.o object/kd.o
	$(cc) -c -o object/chisq_wrapper.o src/utils/chisq_wrapper.cpp

object/cost_fn.o: src/dalex/cost_fn.cpp include/cost_fn.h object/chisq_wrapper.o \
object/ellipse.o
	$(cc) -c -o object/cost_fn.o src/dalex/cost_fn.cpp

object/jellyBean.o: src/utils/jellyBean.cpp include/jellyBean.h object/chisq.o
	$(cc) -c -o object/jellyBean.o src/utils/jellyBean.cpp

object/simplex.o: src/utils/simplex.cpp include/simplex.h object/wrappers.o object/chisq_wrapper.o
	$(cc) -c -o object/simplex.o src/utils/simplex.cpp

object/control_integrator.o: src/controls/control_integrator.cpp \
include/controls/control_integrator.h object/simplex.o object/kd.o
	$(cc) -c -o object/control_integrator.o src/controls/control_integrator.cpp

curved_4d_truth: object/control_integrator.o src/controls/curved_4d_truth.cpp \
object/jellyBean.o include/exampleLikelihoods.h
	$(cc) -o bin/curved_4d_truth src/controls/curved_4d_truth.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/wrappers.o \
	object/chisq.o object/jellyBean.o object/control_integrator.o \
	object/simplex.o \
	$(LIBRARIES)

object/explorers.o: src/dalex/explorers.cpp include/explorers.h \
object/chisq_wrapper.o object/cost_fn.o
	$(cc) -c -o object/explorers.o src/dalex/explorers.cpp

object/dalex.o: src/dalex/dalex.cpp include/dalex.h object/containers.o \
object/goto_tools.o object/chisq_wrapper.o object/simplex.o object/explorers.o \
object/ellipse.o
	$(cc) -c -o object/dalex.o src/dalex/dalex.cpp

object/dalex_initializer.o: src/dalex/dalex_initializer.cpp include/dalex_initializer.h \
object/chisq_wrapper.o object/simplex.o
	$(cc) -c -o object/dalex_initializer.o src/dalex/dalex_initializer.cpp

test_dalex_init: src/tests/test_dalex_init.cpp object/dalex_initializer.o \
include/exampleLikelihoods.h object/jellyBean.o
	$(cc) -o bin/test_dalex_init src/tests/test_dalex_init.cpp \
	object/containers.o object/goto_tools.o \
	object/wrappers.o object/dalex_initializer.o object/simplex.o \
	object/chisq_wrapper.o object/kd.o object/jellyBean.o \
	object/chisq.o

object/dalex_driver.o: src/dalex/dalex_driver.cpp include/dalex_driver.h \
object/simplex.o object/cost_fn.o object/eigen_wrapper.o \
object/dalex.o object/dalex_initializer.o object/ellipse.o
	$(cc) -c -o object/dalex_driver.o src/dalex/dalex_driver.cpp

basic_example: src/examples/basic_example.cpp object/dalex_driver.o \
object/eigen_wrapper.o object/ellipse.o
	$(cc) -o bin/basic_example src/examples/basic_example.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/wrappers.o object/chisq_wrapper.o object/eigen_wrapper.o object/simplex.o \
	object/cost_fn.o object/dalex_driver.o object/ellipse.o \
        object/dalex.o object/dalex_initializer.o object/explorers.o \
	$(LIBRARIES)

curved_4d: src/examples/curved_4d_example.cpp object/dalex_driver.o \
object/jellyBean.o include/exampleLikelihoods.h object/eigen_wrapper.o \
object/ellipse.o
	$(cc) -o bin/curved_4d src/examples/curved_4d_example.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/wrappers.o object/chisq_wrapper.o object/eigen_wrapper.o object/simplex.o \
	object/cost_fn.o object/dalex_driver.o object/jellyBean.o object/ellipse.o \
        object/dalex.o object/dalex_initializer.o object/explorers.o \
	$(LIBRARIES)

curved_12d: src/examples/curved_12d_example.cpp object/dalex_driver.o \
object/jellyBean.o include/exampleLikelihoods.h object/eigen_wrapper.o \
object/ellipse.o
	$(cc) -o bin/curved_12d src/examples/curved_12d_example.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/wrappers.o object/chisq_wrapper.o object/eigen_wrapper.o object/simplex.o \
	object/cost_fn.o object/dalex_driver.o object/jellyBean.o object/ellipse.o \
        object/dalex.o object/dalex_initializer.o object/explorers.o \
	$(LIBRARIES)

ellipse_12d: src/examples/ellipse_12d_example.cpp object/dalex_driver.o \
object/jellyBean.o include/exampleLikelihoods.h object/eigen_wrapper.o \
object/ellipse.o
	$(cc) -o bin/ellipse_12d src/examples/ellipse_12d_example.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/wrappers.o object/chisq_wrapper.o object/eigen_wrapper.o object/simplex.o \
	object/cost_fn.o object/dalex_driver.o object/jellyBean.o object/ellipse.o \
        object/dalex.o object/dalex_initializer.o object/explorers.o \
	$(LIBRARIES)


test_opt: src/examples/test_opt.cpp object/dalex_driver.o \
object/jellyBean.o include/exampleLikelihoods.h object/mcmc.o object/eigen_wrapper.o
	$(cc) -o bin/test_opt src/examples/test_opt.cpp \
	object/containers.o object/goto_tools.o object/kd.o object/chisq.o \
	object/wrappers.o object/chisq_wrapper.o object/eigen_wrapper.o object/simplex.o \
	object/cost_fn.o object/dalex_driver.o object/jellyBean.o object/mcmc.o object/chain.o \
        object/kde.o object/gp_lin.o object/dalex.o \
	object/dalex_initializer.o \
	$(LIBRARIES)

preburner: src/preburner/preburner.cpp object/ellipse.o
	$(cc) -o bin/preburner src/preburner/preburner.cpp \
	object/containers.o object/goto_tools.o object/ellipse.o

object/mcmc.o: include/mcmc.h src/mcmc/mcmc.cpp object/containers.o \
object/goto_tools.o object/eigen_wrapper.o object/wrappers.o
	$(cc) -c -o object/mcmc.o src/mcmc/mcmc.cpp

mcmc_example: src/examples/mcmc_curved_12d_example.cpp object/mcmc.o \
object/chisq.o object/jellyBean.o
	$(cc) -o bin/mcmc_example src/examples/mcmc_curved_12d_example.cpp \
	object/containers.o object/goto_tools.o object/wrappers.o \
	object/eigen_wrapper.o object/chisq.o object/jellyBean.o \
	object/mcmc.o \
	$(LIBRARIES)

ellipse_sampler: src/examples/elliptical_sampling_example.cpp \
object/ellipse.o object/chisq.o object/jellyBean.o
	$(cc) -o bin/ellipse_sampler src/examples/elliptical_sampling_example.cpp \
	object/containers.o object/ellipse.o object/goto_tools.o \
	object/wrappers.o object/chisq.o object/jellyBean.o \
	$(LIBRARIES)

bayes: src/examples/dalex_bayes.cpp object/kd.o include/exampleLikelihoods.h \
object/jellyBean.o
	$(cc) -o bin/bayes src/examples/dalex_bayes.cpp \
	object/containers.o object/goto_tools.o object/kd.o \
	object/wrappers.o object/chisq.o object/jellyBean.o \
	$(LIBRARIES)

hyper_bayes: src/examples/dalex_bayes_hyperbox.cpp object/hyperbox.o
	$(cc) -o bin/hyper_bayes src/examples/dalex_bayes_hyperbox.cpp \
	object/containers.o object/goto_tools.o object/hyperbox.o \
	$(LIBRARIES)

hyper_bayes_12: src/examples/dalex_bayes_hyperbox_12.cpp object/hyperbox.o \
object/jellyBean.o include/exampleLikelihoods.h object/kd.o
	$(cc) -o bin/hyper_bayes_12 src/examples/dalex_bayes_hyperbox_12.cpp \
	object/containers.o object/goto_tools.o object/hyperbox.o \
	object/wrappers.o object/chisq.o object/jellyBean.o \
	object/kd.o \
	$(LIBRARIES)


all:
	make test_containers
	make test_kd
	make test_eigen

clean:
	rm object/*.o bin/*
