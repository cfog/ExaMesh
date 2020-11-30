#  Copyright 2019 by Carl Ollivier-Gooch.  The University of British
#  Columbia disclaims all copyright interest in the software ExaMesh.//
#
#  This file is part of ExaMesh.
#
#  ExaMesh is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  ExaMesh is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with ExaMesh.  If not, see <https://www.gnu.org/licenses/>.

.SILENT:

.SUFFIXES: .cxx

.PHONY: library

CXXOBJECTS=refine.o 

LIBOBJECTS=TetDivider.o PyrDivider.o PrismDivider.o HexDivider.o CellDivider.o \
BdryTriDivider.o BdryQuadDivider.o refinePart.o ExaMesh.o UMesh.o CubicMesh.o GeomUtils.o \
LagrangeMapping.o LengthScaleMapping.o UniformMapping.o \
LagrangeCubicTet.o LagrangeCubicPyr.o LagrangeCubicPrism.o LagrangeCubicHex.o \
Part.o partition.o

OBJECTS=$(CXXOBJECTS) $(LIBOBJECTS)
DEBUG=-g
OPT=-O3 -DNDEBUG -g
OPT_DEBUG=$(OPT) 
CPPFLAGS=-I/home/cfog/Research/Projects/ExaMesh/src -I/home/cfog/Research/External/GMGW/src
CXX_COMPILE=g++ -Wall -Wextra -fPIC $(OPT_DEBUG) $(CPPFLAGS) $(EXTRAFLAGS) -fopenmp
CXX_LINK=g++ -fPIC $(EXTRAFLAGS) $(OPT_DEBUG) -fopenmp
THISDIR=/home/cfog/Research/Projects/ExaMesh/src
MESHIOLIB=-L/home/cfog/Research/External/GMGW/src -Wl,-rpath=/home/cfog/Research/External/GMGW/src -lMeshIO
EXAMESHLIB=-L$(THISDIR) -lexamesh -Wl,-rpath=$(THISDIR)
CGNSLIB=-lcgns -lhdf5_serial -lsz
LDFLAGS=$(EXAMESHLIB) $(MESHIOLIB) $(CGNSLIB)
LIBRARY=libexamesh.so

.cxx.o:
	echo Compiling $*.cxx
	-rm -f $*.o
	-$(CXX_COMPILE) -c $*.cxx || \
	(echo Compile of $*.cxx failed!	Command was: ; \
	echo $(CXX_COMPILE) -c $*.cxx)

all default: $(OBJECTS) library test-exa refine
	./test-exa

library $(LIBRARY) : $(LIBOBJECTS)
	-$(CXX_COMPILE) -shared -o libexamesh.so $(EXTRAFLAGS) $(LIBOBJECTS)

refine: $(CXXOBJECTS) $(LIBRARY)
	echo Linking refine
	-$(CXX_LINK) -o refine $(CXXOBJECTS) $(LDFLAGS) 
	echo " "

test-exa.o: test-exa.cxx 
	echo Compiling $*.cxx
	$(CXX_COMPILE) -DBOOST_TEST_DYN_LINK -Wall -g -c test-exa.cxx 

test-exa: test-exa.o $(LIBRARY)
	$(CXX_LINK) -o test-exa test-exa.o $(LDFLAGS) -lboost_unit_test_framework

clean:
	rm -f *.o refine test-exa libexamesh.so *.gcda *.gcno

depend:
	makedepend $(CPPFLAGS) -Y *.cxx 2> /dev/null

cover:
	lcov --capture --directory . --output coverage/lcov.raw.info
	lcov --remove coverage/lcov.raw.info "/usr*" "*GMGW1*" -o coverage/lcov.info
	genhtml --demangle-cpp --output-directory coverage coverage/lcov.info
	firefox coverage/index.html

# DO NOT DELETE

BdryQuadDivider.o: BdryQuadDivider.h CellDivider.h ExaMesh.h Mapping.h
BdryQuadDivider.o: exa-defs.h Part.h UMesh.h CubicMesh.h
BdryTriDivider.o: BdryTriDivider.h CellDivider.h ExaMesh.h Mapping.h
BdryTriDivider.o: exa-defs.h Part.h UMesh.h CubicMesh.h
CellDivider.o: ExaMesh.h Mapping.h exa-defs.h Part.h GeomUtils.h
CellDivider.o: CellDivider.h UMesh.h CubicMesh.h
CubicMesh.o: CubicMesh.h ExaMesh.h Mapping.h exa-defs.h Part.h UMesh.h
ExaMesh.o: ExaMesh.h Mapping.h exa-defs.h Part.h GeomUtils.h UMesh.h
ExaMesh.o: CubicMesh.h
HexDivider.o: HexDivider.h CellDivider.h ExaMesh.h Mapping.h exa-defs.h
HexDivider.o: Part.h UMesh.h CubicMesh.h
LagrangeCubicHex.o: Mapping.h exa-defs.h
LagrangeCubicPrism.o: Mapping.h exa-defs.h
LagrangeCubicPyr.o: Mapping.h exa-defs.h
LagrangeCubicTet.o: Mapping.h exa-defs.h
LagrangeMapping.o: ExaMesh.h Mapping.h exa-defs.h Part.h
LengthScaleMapping.o: ExaMesh.h Mapping.h exa-defs.h Part.h
Part.o: exa-defs.h Part.h
PrismDivider.o: PrismDivider.h CellDivider.h ExaMesh.h Mapping.h exa-defs.h
PrismDivider.o: Part.h UMesh.h CubicMesh.h
PyrDivider.o: GeomUtils.h PyrDivider.h CellDivider.h ExaMesh.h Mapping.h
PyrDivider.o: exa-defs.h Part.h UMesh.h CubicMesh.h
TetDivider.o: GeomUtils.h TetDivider.h CellDivider.h ExaMesh.h Mapping.h
TetDivider.o: exa-defs.h Part.h UMesh.h CubicMesh.h
UMesh.o: ExaMesh.h Mapping.h exa-defs.h Part.h UMesh.h CubicMesh.h
UMesh.o: /home/cfog/Research/External/GMGW/src/GMGW_unstr.hxx
UMesh.o: /home/cfog/Research/External/GMGW/src/config.h
UMesh.o: /home/cfog/Research/External/GMGW/src/GMGW_FileWrapper.hxx
UniformMapping.o: ExaMesh.h Mapping.h exa-defs.h Part.h
partition.o: ExaMesh.h Mapping.h exa-defs.h Part.h
refine.o: ExaMesh.h Mapping.h exa-defs.h Part.h CubicMesh.h UMesh.h
refinePart.o: ExaMesh.h Mapping.h exa-defs.h Part.h HexDivider.h
refinePart.o: CellDivider.h UMesh.h CubicMesh.h PrismDivider.h PyrDivider.h
refinePart.o: TetDivider.h BdryTriDivider.h BdryQuadDivider.h
test-exa.o: ExaMesh.h Mapping.h exa-defs.h Part.h UMesh.h CubicMesh.h
test-exa.o: TetDivider.h CellDivider.h
