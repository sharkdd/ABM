OS := $(shell uname)

## This should be updated to point at your adevs installation
##DEVSPATH=$(HOME)/Temp/adevs-2.5
LIBPATH=/usr/local/include

CXX = gcc
INCFLAGS = -I$(LIBPATH)/adves-2.6/include -I$(LIBPATH) -I/Users/tata/Documents/code/IBM/src


## Note that the GNU scientific library is required for random variate
## generation
LDFLAGS = -L$(LIBPATH)/gsl -lgsl -ladevs -lm

## ## Compile flags: ##
##
## If GERMS partnership formation is desired, then the -DGERMS flag
## must be set, and the Makefile must specify that the executable
## depends on PshipClassic.o rather than Partnership.o
##
## To replicate the KZN ODE, the IBM must be compiled _without_ -DGERMS 
##
## CODEFLAGS += -DHISTORY         ## store event histories (discarded otherwise to reduce memory requirements)
## CODEFLAGS += -DPARTNER_HISTORY ## with HISTORY, extends histories to include all partnership formation and dissolution
## CODEFLAGS += -DGERMS        ## GERMS partnership module used (PshipClassic otherwise)
CODEFLAGS += -DTYPE_TWO_AGE    ## aging simulated as a poisson process (deterministic otherwise)
CODEFLAGS += -DGEOMETRIC       ## partnership formation is the geometric mean of partners' preferences (arithmetic mean otherwise)

## DEBUGFLAGS = -Wall -g
DEBUGFLAGS = -Wall -O2
ARCHFLAGS = -D$(OS)
CXXFLAGS += $(DEBUGFLAGS) $(CODEFLAGS) $(ARCHFLAGS) $(INCFLAGS)

simulate: simulate.cpp PshipClassic.o Bookkeeper.o Person.o ModelTypes.o Demographer.o Progression.o Transmission.o Circumcision.o Behavior.o Treatment.o Resistance.o

.o: $@.cpp $@.H
	$(CXX) $(CXXFLAGS) -stdlib=libstdc++ -c $@ $^ 

.cpp: $@.cpp
	$(CXX) $(CXXFLAGS) -stdlib=libstdc++ -o $@ $^ $(LDFLAGS) 

clean:
	rm -rf *.dSYM *.o *~
