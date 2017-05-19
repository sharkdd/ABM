OS := $(shell uname)

CXX = g++
INCFLAGS = -I.
LDFLAGS =  -lgsl -lgslcblas -lm

## Compile flags:
## CD4_500 split early chronic infection (CD4>350) into CD4>500 and 350<CD4<500
## CODEFLAGS += -DCD4_500

## DEBUGFLAGS = -Wall -g
DEBUGFLAGS = -Wall -O2
ARCHFLAGS = -D$(OS)
CXXFLAGS += $(DEBUGFLAGS) $(CODEFLAGS) $(ARCHFLAGS) $(INCFLAGS)

simulate: simulate.cpp Model.o attributes.o

.o: $@.cpp $@.H
	$(CXX) $(CXXFLAGS) -c $@ $^

.cpp: $@.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -rf *.dSYM *.o *~
