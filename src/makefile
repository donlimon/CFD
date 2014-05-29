IDIR =../include
EIGENDIR =/usr/include/eigen3
CXX=g++
CPPFLAGS=-I$(IDIR) -I$(EIGENDIR)
BINARY=fracstep

ODIR=obj
LDIR =../lib
BDIR=../bin

LIBS=

_DEPS = settings.hpp matinc.hpp matrixgen.hpp io.hpp taylorgreen.hpp testing.hpp vortex.hpp structs.hpp
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o matrixgen.o io.o taylorgreen.o testing.o vortex.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

debug: CPPFLAGS += -g
debug: $(OBJ)
	$(CXX) -o $(BDIR)/$(BINARY)-DEBUG $^ $(CPPFLAGS) $(LIBS)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CPPFLAGS)

release: CPPFLAGS += -O3 -DNDEBUG -ftree-vectorize -msse2 -fopenmp
release: $(OBJ)
	$(CXX) -o $(BDIR)/$(BINARY) $^ $(CPPFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
