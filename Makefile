.SILENT:

.SUFFIXES: .cxx

.PHONY: library

CXXOBJECTS=examesh.o refinePart.o UMesh.o GeomUtils.o 

LIBOBJECTS=TetDivider.o PyrDivider.o PrismDivider.o HexDivider.o CellDivider.o \
BdryTriDivider.o BdryQuadDivider.o

OBJECTS=$(CXXOBJECTS) $(LIBOBJECTS)
OPT_DEBUG=-O3
CPPFLAGS=-I/home/cfog/Research/Projects/ExaMesh -I/home/cfog/GMGW1/src
CXX_COMPILE=g++ -Wall -Wextra -fPIC $(OPT_DEBUG) $(CPPFLAGS) $(EXTRAFLAGS)
CXX_LINK=g++ -fPIC $(EXTRAFLAGS) $(OPT_DEBUG)
THISDIR=/home/cfog/Research/Projects/ExaMesh/src
MESHIOLIB=-L/home/cfog/GMGW1/src -Wl,-rpath=/home/cfog/GMGW1/src -lMeshIO
EXAMESHLIB=-L$(THISDIR) -lexamesh -Wl,-rpath=$(THISDIR)
LDFLAGS=$(EXAMESHLIB) $(MESHIOLIB)

.cxx.o:
	echo Compiling $*.cxx
	-rm -f $*.o
	-$(CXX_COMPILE) -c $*.cxx || \
	(echo Compile of $*.cxx failed!	Command was: ; \
	echo $(CXX_COMPILE) -c $*.cxx)

all default: $(OBJECTS) library test-exa examesh
	./test-exa

library: $(LIBOBJECTS)
	-$(CXX_COMPILE) -shared -o libexamesh.so $(EXTRAFLAGS) $(OBJECTS)

examesh: $(CXXOBJECTS)
	echo Linking examesh
	-$(CXX_LINK) -o examesh $(CXXOBJECTS) $(LDFLAGS) 
	echo " "

test-exa.o: test-exa.cxx
	echo Compiling $*.cxx
	$(CXX_COMPILE) -DBOOST_TEST_DYN_LINK -Wall -g -c test-exa.cxx 

test-exa: test-exa.o libexamesh.so
	$(CXX_LINK) -o test-exa test-exa.o $(LDFLAGS) -lboost_unit_test_framework

clean:
	rm -f *.o examesh

# DO NOT DELETE

BdryQuadDivider.o: BdryQuadDivider.h CellDivider.h examesh.h UMesh.h
BdryTriDivider.o: BdryTriDivider.h CellDivider.h examesh.h UMesh.h
CellDivider.o: examesh.h CellDivider.h UMesh.h
HexDivider.o: HexDivider.h examesh.h CellDivider.h UMesh.h
PrismDivider.o: PrismDivider.h examesh.h CellDivider.h UMesh.h
PyrDivider.o: GeomUtils.h PyrDivider.h examesh.h CellDivider.h UMesh.h
refinePart.o: HexDivider.h examesh.h CellDivider.h UMesh.h PrismDivider.h
refinePart.o: PyrDivider.h TetDivider.h BdryTriDivider.h BdryQuadDivider.h
test-exa.o: examesh.h UMesh.h
TetDivider.o: GeomUtils.h TetDivider.h examesh.h CellDivider.h UMesh.h
UMesh.o: examesh.h UMesh.h
