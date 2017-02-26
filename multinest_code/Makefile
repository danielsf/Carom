NESTLIBDIR = ../
LIBS =  -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB) -lc++
APSDIR = /Users/danielsf/physics/Carom/
 
OBJFILES = jellyBean.o \
        $(APSDIR)object/containers.o \
        $(APSDIR)object/goto_tools.o \
        $(APSDIR)object/kd.o \
        $(APSDIR)object/chisq.o \
        $(APSDIR)object/jellyBean.o \
        $(APSDIR)object/wrappers.o

all: jellyBean

%.o: %.cc
	$(CXX) $(CFLAGS) -c $*.cc
 
jellyBean: $(OBJFILES)
	$(FC) $(FFLAGS) -I$(APSDIR)include/ -o ../jellyBean $(OBJFILES) $(LIBS)

clean:
	rm -f *.o *.mod ../jellyBean
