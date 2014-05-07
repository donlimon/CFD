IDIR =../include
CXX=g++
CPPFLAGS=-I$(IDIR)

ODIR=obj
LDIR =../lib
BDIR=../bin

LIBS=

_DEPS =
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CPPFLAGS)

debug: $(OBJ)
	$(CXX) -o $(BDIR)/fracstep-DEBUG $^ -g $(CPPFLAGS) $(LIBS)
	
release: $(OBJ)
	$(CXX) -o $(BDIR)/fracstep $^ $(CPPFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
