IDIR =../include
EIGENDIR =/usr/include/eigen3
CXX=g++
CPPFLAGS=-I$(IDIR) -I$(EIGENDIR)
BINARY=fracstep

ODIR=obj
LDIR =../lib
BDIR=../bin

LIBS=

_DEPS = settings.hpp matrixgen.hpp
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o matrixgen.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

debug: CPPFLAGS += -g
debug: $(OBJ)
	$(CXX) -o $(BDIR)/$(BINARY)-DEBUG $^ $(CPPFLAGS) $(LIBS)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CPPFLAGS)

release: $(OBJ)
	$(CXX) -o $(BDIR)/$(BINARY) $^ $(CPPFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
