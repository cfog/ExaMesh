//  Copyright 2019 by Carl Ollivier-Gooch.  The University of British
//  Columbia disclaims all copyright interest in the software ExaMesh.//
//
//  This file is part of ExaMesh.
//
//  ExaMesh is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as
//  published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  ExaMesh is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with ExaMesh.  If not, see <https://www.gnu.org/licenses/>.

/*
 * test-exa.cxx
 *
 *  Created on: Jul. 5, 2019
 *      Author: cfog
 */

#define BOOST_TEST_MODULE test-exa
#include <boost/test/unit_test.hpp>
#include "ExaMesh.h"
#include "UMesh.h"
#include "CubicMesh.h"
#include "TetDivider.h"
#include "PyrDivider.h"
#include "PrismDivider.h"
#include "HexDivider.h"
#include "Mapping.h"

#define DO_SUBDIVISION_TESTS

#ifdef DO_SUBDIVISION_TESTS

static void checkExpectedSize(const UMesh &UM) {
	BOOST_CHECK_EQUAL(UM.maxNVerts(), UM.numVerts());
	BOOST_CHECK_EQUAL(UM.maxNBdryTris(), UM.numBdryTris());
	BOOST_CHECK_EQUAL(UM.maxNBdryQuads(), UM.numBdryQuads());
	BOOST_CHECK_EQUAL(UM.maxNTets(), UM.numTets());
	BOOST_CHECK_EQUAL(UM.maxNPyrs(), UM.numPyramids());
	BOOST_CHECK_EQUAL(UM.maxNPrisms(), UM.numPrisms());
	BOOST_CHECK_EQUAL(UM.maxNHexes(), UM.numHexes());
}
#endif

void SetArtificialIntVertQuad(QuadFaceVerts &quad, const emInt nDivs){
	// Numbering from 1 to .. 
	emInt k=1; 	
	for (int jj = 0; jj <= nDivs ; jj++) {
	 	for (int ii = 0; ii <= nDivs; ii++) {
			quad.setIntVertInd(ii,jj,k);
			k++; 
		}
	}

}
void setExpectedMapping (const emInt rotation, 
std::unordered_map<emInt,emInt> &map){
	switch (rotation) {
  		case -1:
   		 	map[1]=1; 
			map[2]=5;
			map[3]=9;
			map[4]=13; 
			map[5]=2;
			map[6]=6;
			map[7]=10;
			map[8]=14; 
			map[9]=3;
			map[10]=7;
			map[11]=11;
			map[12]=15;
			map[13]=4; 
			map[14]=8;
			map[15]=12;
			map[16]=16;
   	 	break;
  		case -2:
		   	map[1]=4; 
			map[2]=3;
			map[3]=2;
			map[4]=1; 
			map[5]=8;
			map[6]=7;
			map[7]=6;
			map[8]=5; 
			map[9]=12;
			map[10]=11;
			map[11]=10;
			map[12]=9;
			map[13]=16; 
			map[14]=15;
			map[15]=14;
			map[16]=13;
    		
    	break;
  		case -3:
			map[1]=16; 
			map[2]=12;
			map[3]=8;
			map[4]=4; 
			map[5]=15;
			map[6]=11;
			map[7]=7;
			map[8]=3; 
			map[9]=14;
			map[10]=10;
			map[11]=6;
			map[12]=2;
			map[13]=13; 
			map[14]=9;
			map[15]=5;
			map[16]=1;
    		
    	break;
  		case -4:
			map[1]=13; 
			map[2]=14;
			map[3]=15;
			map[4]=16; 
			map[5]=9;
			map[6]=10;
			map[7]=11;
			map[8]=12; 
			map[9]=5;
			map[10]=6;
			map[11]=7;
			map[12]=8;
			map[13]=1; 
			map[14]=2;
			map[15]=3;
			map[16]=4;
    		
    	break;
  		case 2:
			map[1]=13; 
			map[2]=9;
			map[3]=5;
			map[4]=1; 
			map[5]=14;
			map[6]=10;
			map[7]=6;
			map[8]=2; 
			map[9]=15;
			map[10]=11;
			map[11]=7;
			map[12]=3;
			map[13]=16; 
			map[14]=12;
			map[15]=8;
			map[16]=4 ;
    		
    	break;
  		case 3:
			map[1]=16; 
			map[2]=15;
			map[3]=14;
			map[4]=13; 
			map[5]=12;
			map[6]=11;
			map[7]=10;
			map[8]=9; 
			map[9]=8;
			map[10]=7;
			map[11]=6;
			map[12]=5;
			map[13]=4; 
			map[14]=3;
			map[15]=2;
			map[16]=1;
    		
    	break;
  		case 4:
			map[1]=4; 
			map[2]=8;
			map[3]=12;
			map[4]=16; 
			map[5]=3;
			map[6]=7;
			map[7]=11;
			map[8]=15; 
			map[9]=2;
			map[10]=6;
			map[11]=10;
			map[12]=14;
			map[13]=1; 
			map[14]=5;
			map[15]=9;
			map[16]=13;
    		
    	break;
	}
}
void setArbitraryTriDataForTesting (const emInt i ,emInt (&local)[3], emInt (&global)[3], 
emInt (&remote)[3], emInt &nDivs, emInt &partId, emInt &remoteId, emInt &type, 
emInt &elemInd, bool &globalCompare){
	local[0] = i;
    local[1] = 2 * i;
    local[2] = 3 * i;
    
    global[0] = 2 * i;
    global[1] = 3 * i;
    global[2] = 4 * i;
    
    remote[0] = i+1;
    remote[1] = i + 2;
    remote[2] = i + 3;

	nDivs= i ; 
	partId= 2*i; 
	remoteId= 3*i; 
	type= 4*i; 
	elemInd= 5*i;
	

	if(i%2==0){
		globalCompare=true; 
	}else{
		globalCompare=false;
	}
}
void setArbitraryQuadDataForTesting (const emInt i ,emInt (&local)[4], emInt (&global)[4], 
emInt (&remote)[4], emInt &nDivs, emInt &partId, emInt &remoteId, emInt &type, 
emInt &elemInd, bool &globalCompare){
	local[0] = i;
    local[1] = 2 * i;
    local[2] = 3 * i;
	local[3]=  4 * i; 
    
    global[0] = 2 * i;
    global[1] = 3 * i;
    global[2] = 4 * i;
	global[3]=  5*i ; 
    
    remote[0] = i+1;
    remote[1] = i + 2;
    remote[2] = i + 3;
	remote[3]=  i + 4 ; 

	nDivs= i ; 
	partId= 2*i; 
	remoteId= 3*i; 
	type= 4*i; 
	elemInd= 5*i;
	

	if(i%2==0){
		globalCompare=true; 
	}else{
		globalCompare=false;
	}
}
void setArbitrary_CellPartData_ForTesting (const emInt i ,double (&coords)[3],
 emInt &index, emInt &cellType){

    
    coords[0] = i+1;
    coords[1] = i + 2;
    coords[2] = i + 3;

	index= 10*i; 
	cellType= 20*i; 
}
void setArbitrary_Part_DataForTesting (const emInt i , double &xmin, double &xmax, 
double &ymin, double &ymax, double &zmin, double &zmax,
emInt &first, emInt &last, emInt &parts ){
	xmin= 1.2*i; 
	xmax= 2.4*i ; 

	ymin= 3.2*i; 
	ymax= 1.2+i; 

	zmin= 2.3+i; 
	zmax= 5.2*i; 

	first= i; 
	last=  i+1; 
	parts= i+2; 
}
struct MixedMeshFixture {
	UMesh *pUM_In, *pUM_Out;
	exa_map<Edge, EdgeVerts> vertsOnEdges;
	exa_set<TriFaceVerts> vertsOnTris;
	exa_set<QuadFaceVerts> vertsOnQuads;
	MixedMeshFixture() {
		pUM_In = new UMesh(11, 11, 6, 6, 1, 1, 1, 1);
		double coords[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 },
				{ 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, -1 }, { 1, 0, -1 },
				{ 1, 1, -1 }, { 0, 1, -1 }, { 0, -1, 0 }, { 0, -1, -1 } };
		emInt triVerts[][3] = { { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 },
				{ 0, 9, 4 }, { 9, 1, 4 }, { 10, 6, 5 } };
		emInt quadVerts[][4] = { { 6, 7, 2, 1 }, { 7, 8, 3, 2 }, { 8, 5, 0, 3 },
				{ 10, 6, 1, 9 }, { 5, 10, 9, 0 }, { 5, 6, 7, 8 } };
		emInt tetVerts[4] = { 9, 1, 0, 4 };
		emInt pyrVerts[5] = { 0, 1, 2, 3, 4 };
		emInt prismVerts[6] = { 10, 6, 5, 9, 1, 0 };
		emInt hexVerts[8] = { 5, 6, 7, 8, 0, 1, 2, 3 };

		for (int ii = 0; ii < 11; ii++) {
			pUM_In->addVert(coords[ii]);
		}
		for (int ii = 0; ii < 6; ii++) {
			pUM_In->addBdryTri(triVerts[ii]);
			pUM_In->addBdryQuad(quadVerts[ii]);
		}
		pUM_In->addTet(tetVerts);
		pUM_In->addPyramid(pyrVerts);
		pUM_In->addPrism(prismVerts);
		pUM_In->addHex(hexVerts);

		MeshSize MSOut = pUM_In->computeFineMeshSize(4);
		pUM_Out = new UMesh(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris,
				MSOut.nBdryQuads, MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms,
				MSOut.nHexes);
		for (int ii = 0; ii < 11; ii++) {
			pUM_Out->addVert(coords[ii]);
		}
	}
	~MixedMeshFixture() {
		delete pUM_In;
	}
};

void makeLengthScaleUniform(const UMesh* pUM_In) {
	// Check the uniform length scale cases.
	for (int ii = 0; ii < pUM_In->numVerts(); ii++) {
		pUM_In->setLengthScale(ii, 1.);
	}
}
void setPrescribedLengthScale(const UMesh* pUM_In) {
	// A simple analytic length scale.
	for (int ii = 0; ii < pUM_In->numVerts(); ii++) {
		double len = 1 + 0.1*(pUM_In->getX(ii) + pUM_In->getY(ii) + pUM_In->getZ(ii));
		pUM_In->setLengthScale(ii, len);
	}
}

// The following test suite confirms correctness of the parametric-to-
// physical space mapping for cubic Lagrange element shapes.
BOOST_AUTO_TEST_SUITE(LagrangeCubicMappingTest)
static void cubicXYZ(const double uvw[3], double xyz[3]) {
	const double &u = uvw[0];
	const double &v = uvw[1];
	const double &w = uvw[2];
	xyz[0] = u * u * u + 2 * u * u * v + 4 * u * v * v - 30 * u * v * w
			+ 5 * u * u + 6 * u - 3 * u * v + 3;
	xyz[1] = v * v * v + 2 * v * v * w + 4 * v * w * w - 54 * u * v * w
			+ 5 * v * v + 6 * v - 4 * v * w + 7;
	xyz[2] = w * w * w + 2 * w * w * u + 4 * w * u * u - 72 * u * v * w
			+ 5 * w * w + 7 * w - 5 * w * u + 10.5;
}

BOOST_AUTO_TEST_CASE(CubicTet) {
	printf("Testing cubic tet mapping\n");
	CubicMesh CM(0, 0, 0, 0, 0, 0, 0, 0);
	LagrangeCubicTetMapping LCTM(&CM);

	// So the nodal values are just these functions evaluated at nodes, where for the
	// canonical tet, we have:
	double uvw[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, // Corner nodes
			{ 1 / 3., 0, 0 }, { 2 / 3., 0, 0 },  // Nodes between 0 and 1
			{ 2 / 3., 1 / 3., 0 }, { 1 / 3., 2 / 3., 0 }, // Nodes between 1 and 2
			{ 0, 2 / 3., 0 }, { 0, 1 / 3., 0 }, // Nodes between 2 and 0
			{ 0, 0, 1 / 3. }, { 0, 0, 2 / 3. }, // Nodes between 0 and 3
			{ 2 / 3., 0, 1 / 3. }, { 1 / 3., 0, 2 / 3. }, // Nodes between 1 and 3
			{ 0, 2 / 3., 1 / 3. }, { 0, 1 / 3., 2 / 3. }, // Nodes between 2 and 3
			{ 1 / 3., 1 / 3., 0 }, // Node on 012
			{ 1 / 3., 0, 1 / 3. }, // Node on 013
			{ 1 / 3., 1 / 3., 1 / 3. }, // Node on 123
			{ 0, 1 / 3., 1 / 3. } }; // Node on 203
	double xyz[20][3] = { 0 };

	// Now the nodal values, by simple functional evaluation.
	for (int ii = 0; ii < 20; ii++) {
		cubicXYZ(uvw[ii], xyz[ii]);
	}
	LCTM.setNodalValues(xyz);
	LCTM.setModalValues();

	double testUVW[] = { 1. / M_PI, 1. / M_E, 0.5 * (1 - 1. / M_PI - 1. / M_E) };
	double LCTMxyz[3], funcxyz[3];
	LCTM.computeTransformedCoords(testUVW, LCTMxyz);
	cubicXYZ(testUVW, funcxyz);
	BOOST_CHECK_CLOSE(LCTMxyz[0], funcxyz[0], 1.e-8);
	BOOST_CHECK_CLOSE(LCTMxyz[1], funcxyz[1], 1.e-8);
	BOOST_CHECK_CLOSE(LCTMxyz[2], funcxyz[2], 1.e-8);

	for (int ii = 0; ii < 20; ii++) {
		LCTM.computeTransformedCoords(uvw[ii], LCTMxyz);
//			printf("%d\n", ii);
		BOOST_CHECK_CLOSE(LCTMxyz[0], xyz[ii][0], 1.e-8);
		BOOST_CHECK_CLOSE(LCTMxyz[1], xyz[ii][1], 1.e-8);
		BOOST_CHECK_CLOSE(LCTMxyz[2], xyz[ii][2], 1.e-8);
	}

}

BOOST_AUTO_TEST_CASE(CubicPyramid) {
	printf("Testing cubic pyramid mapping\n");
	CubicMesh CM(0, 0, 0, 0, 0, 0, 0, 0);
	LagrangeCubicPyramidMapping LCPM(&CM);

	// So the nodal values are just these functions evaluated at nodes, where for the
	// canonical pyramid, we have:
	double uvw[][3] = { { -1, -1, 0 }, { 1, -1, 0 }, { 1, 1, 0 }, { -1, 1, 0 },
			{ 0, 0, 1 }, // Corner nodes
			{ -1 / 3., -1, 0 },
			{ 1 / 3., -1, 0 }, // Nodes between 0 and 1
			{ 1, -1 / 3., 0 },
			{ 1, 1 / 3., 0 }, // Nodes between 1 and 2
			{ 1 / 3., 1, 0 },
			{ -1 / 3., 1, 0 }, // Nodes between 2 and 3
			{ -1, 1 / 3., 0 },
			{ -1, -1 / 3., 0 }, // Nodes between 3 and 0
			{ -2 / 3., -2 / 3., 1 / 3. },
			{ -1 / 3., -1 / 3., 2 / 3. }, // Nodes between 0 and 4
			{ 2 / 3., -2 / 3., 1 / 3. },
			{ 1 / 3., -1 / 3., 2 / 3. }, // Nodes between 1 and 4
			{ 2 / 3., 2 / 3., 1 / 3. },
			{ 1 / 3., 1 / 3., 2 / 3. }, // Nodes between 2 and 4
			{ -2 / 3., 2 / 3., 1 / 3. },
			{ -1 / 3., 1 / 3., 2 / 3. }, // Nodes between 3 and 4
			{ -1 / 3., -1 / 3., 0 }, { 1 / 3., -1 / 3., 0 },
			{ 1 / 3., 1 / 3., 0 }, { -1 / 3., 1 / 3., 0 }, // Nodes on the base
			{ 0, -2 / 3., 1 / 3. }, // Node on 014
			{ 2 / 3., 0, 1 / 3. }, // Node on 124
			{ 0, 2 / 3., 1 / 3. }, // Node on 234
			{ -2 / 3., 0, 1 / 3. }, // Node on 304
			{ 0, 0, 1 / 3. } }; // Node in the interior

	double uvwToMap[30][3];

	double xyz[30][3] = { 0 };

	// Now the nodal values, by simple functional evaluation.
	for (int ii = 0; ii < 30; ii++) {
		double w = uvw[ii][2];
		uvwToMap[ii][0] = (uvw[ii][0] + (1 - w)) / 2;
		uvwToMap[ii][1] = (uvw[ii][1] + (1 - w)) / 2;
		uvwToMap[ii][2] = w;
		cubicXYZ(uvw[ii], xyz[ii]);
	}
	LCPM.setNodalValues(xyz);
	LCPM.setModalValues();

	double testUVW[] = { 1. / M_PI, 1. / M_E, 2 * (1 - 1. / M_PI - 1. / M_E) };
	double LCPMxyz[3], funcxyz[3];
	cubicXYZ(testUVW, funcxyz);
	testUVW[0] = (testUVW[0] + (1 - testUVW[2])) / 2;
	testUVW[1] = (testUVW[1] + (1 - testUVW[2])) / 2;
	LCPM.computeTransformedCoords(testUVW, LCPMxyz);
	BOOST_CHECK_CLOSE(LCPMxyz[0], funcxyz[0], 1.e-8);
	BOOST_CHECK_CLOSE(LCPMxyz[1], funcxyz[1], 1.e-8);
	BOOST_CHECK_CLOSE(LCPMxyz[2], funcxyz[2], 1.e-8);

	for (int ii = 0; ii < 30; ii++) {
		LCPM.computeTransformedCoords(uvwToMap[ii], LCPMxyz);
//			printf("%d\n", ii);
		BOOST_CHECK_CLOSE(LCPMxyz[0], xyz[ii][0], 1.e-8);
		BOOST_CHECK_CLOSE(LCPMxyz[1], xyz[ii][1], 1.e-8);
		BOOST_CHECK_CLOSE(LCPMxyz[2], xyz[ii][2], 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(CubicPrism) {
	printf("Testing cubic prism mapping\n");
	CubicMesh CM(0, 0, 0, 0, 0, 0, 0, 0);
	LagrangeCubicPrismMapping LCPM(&CM);

	// So the nodal values are just these functions evaluated at nodes, where for the
	// canonical prism, we have:
	double uvw[][3] = { { 0, 0, 0 }, // 0
			{ 1, 0, 0 }, // 1
			{ 0, 1, 0 }, // 2
			{ 0, 0, 1 }, // 3
			{ 1, 0, 1 }, // 4
			{ 0, 1, 1 }, // 5
			{ 1. / 3, 0, 0 }, // 6
			{ 2. / 3, 0, 0 }, // 7
			{ 2. / 3, 1. / 3, 0 }, // 8
			{ 1. / 3, 2. / 3, 0 }, // 9
			{ 0, 2. / 3, 0 }, // 10
			{ 0, 1. / 3, 0 }, // 11
			{ 0, 0, 1. / 3 }, // 12
			{ 0, 0, 2. / 3 }, // 13
			{ 1, 0, 1. / 3 }, // 14
			{ 1, 0, 2. / 3 }, // 15
			{ 0, 1, 1. / 3 }, // 16
			{ 0, 1, 2. / 3 }, // 17
			{ 1. / 3, 0, 1 }, // 18
			{ 2. / 3, 0, 1 }, // 19
			{ 2. / 3, 1. / 3, 1 }, // 20
			{ 1. / 3, 2. / 3, 1 }, // 21
			{ 0, 2. / 3, 1 }, // 22
			{ 0, 1. / 3, 1 }, // 23
			{ 1. / 3, 1. / 3, 0 }, // 24
			{ 1. / 3, 0, 1. / 3 }, // 25
			{ 2. / 3, 0, 1. / 3 }, // 26
			{ 2. / 3, 0, 2. / 3 }, // 27
			{ 1. / 3, 0, 2. / 3 }, // 28
			{ 2. / 3, 1. / 3, 1. / 3 }, // 29
			{ 1. / 3, 2. / 3, 1. / 3 }, // 30
			{ 1. / 3, 2. / 3, 2. / 3 }, // 31
			{ 2. / 3, 1. / 3, 2. / 3 }, // 32
			{ 0, 2. / 3, 1. / 3 }, // 33
			{ 0, 1. / 3, 1. / 3 }, // 34
			{ 0, 1. / 3, 2. / 3 }, // 35
			{ 0, 2. / 3, 2. / 3 }, // 36
			{ 1. / 3, 1. / 3, 1 }, // 37
			{ 1. / 3, 1. / 3, 1. / 3 }, // 38
			{ 1. / 3, 1. / 3, 2. / 3 } // 39
	};

	double xyz[40][3] = { 0 };

	// Now the nodal values, by simple functional evaluation.
	for (int ii = 0; ii < 40; ii++) {
		cubicXYZ(uvw[ii], xyz[ii]);
	}
	LCPM.setNodalValues(xyz);
	LCPM.setModalValues();

	double testUVW[] = { 1. / M_PI, 1. / M_E, 0.5 * (1 - 1. / M_PI - 1. / M_E) };
	double LCPMxyz[3], funcxyz[3];
	cubicXYZ(testUVW, funcxyz);
	LCPM.computeTransformedCoords(testUVW, LCPMxyz);
	BOOST_CHECK_CLOSE(LCPMxyz[0], funcxyz[0], 1.e-8);
	BOOST_CHECK_CLOSE(LCPMxyz[1], funcxyz[1], 1.e-8);
	BOOST_CHECK_CLOSE(LCPMxyz[2], funcxyz[2], 1.e-8);

	for (int ii = 0; ii < 40; ii++) {
		LCPM.computeTransformedCoords(uvw[ii], LCPMxyz);
//			printf("%d\n", ii);
		BOOST_CHECK_CLOSE(LCPMxyz[0], xyz[ii][0], 1.e-8);
		BOOST_CHECK_CLOSE(LCPMxyz[1], xyz[ii][1], 1.e-8);
		BOOST_CHECK_CLOSE(LCPMxyz[2], xyz[ii][2], 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(CubicHex) {
	printf("Testing cubic hex mapping\n");
	CubicMesh CM(0, 0, 0, 0, 0, 0, 0, 0);
	LagrangeCubicHexMapping LCHM(&CM);

	// So the nodal values are just these functions evaluated at nodes, where for the
	// canonical pyramid, we have:
	double uvw[][3] = {
	// Bottom layer
			{ 0, 0, 0 }, // 0
			{ 1, 0, 0 }, // 1
			{ 1, 1, 0 }, // 2
			{ 0, 1, 0 }, // 3
			{ 0, 0, 1 }, // 4
			{ 1, 0, 1 }, // 5
			{ 1, 1, 1 }, // 6
			{ 0, 1, 1 }, // 7
			{ 1. / 3, 0, 0 }, // 8
			{ 2. / 3, 0, 0 }, // 9
			{ 1, 1. / 3, 0 }, // 10
			{ 1, 2. / 3, 0 }, // 11
			{ 2. / 3, 1, 0 }, // 12
			{ 1. / 3, 1, 0 }, // 13
			{ 0, 2. / 3, 0 }, // 14
			{ 0, 1. / 3, 0 }, // 15
			{ 0, 0, 1. / 3 }, // 16
			{ 0, 0, 2. / 3 }, // 17
			{ 1, 0, 1. / 3 }, // 18
			{ 1, 0, 2. / 3 }, // 19
			{ 1, 1, 1. / 3 }, // 20
			{ 1, 1, 2. / 3 }, // 21
			{ 0, 1, 1. / 3 }, // 22
			{ 0, 1, 2. / 3 }, // 23
			{ 1. / 3, 0, 1 }, // 24
			{ 2. / 3, 0, 1 }, // 25
			{ 1, 1. / 3, 1 }, // 26
			{ 1, 2. / 3, 1 }, // 27
			{ 2. / 3, 1, 1 }, // 28
			{ 1. / 3, 1, 1 }, // 29
			{ 0, 2. / 3, 1 }, // 30
			{ 0, 1. / 3, 1 }, // 31
			{ 1. / 3, 1. / 3, 0 }, // 32
			{ 2. / 3, 1. / 3, 0 }, // 33
			{ 2. / 3, 2. / 3, 0 }, // 34
			{ 1. / 3, 2. / 3, 0 }, // 35
			{ 1. / 3, 0, 1. / 3 }, // 36
			{ 2. / 3, 0, 1. / 3 }, // 37
			{ 2. / 3, 0, 2. / 3 }, // 38
			{ 1. / 3, 0, 2. / 3 }, // 39
			{ 1, 1. / 3, 1. / 3 }, // 40
			{ 1, 2. / 3, 1. / 3 }, // 41
			{ 1, 2. / 3, 2. / 3 }, // 42
			{ 1, 1. / 3, 2. / 3 }, // 43
			{ 2. / 3, 1, 1. / 3 }, // 44
			{ 1. / 3, 1, 1. / 3 }, // 45
			{ 1. / 3, 1, 2. / 3 }, // 46
			{ 2. / 3, 1, 2. / 3 }, // 47
			{ 0, 2. / 3, 1. / 3 }, // 48
			{ 0, 1. / 3, 1. / 3 }, // 49
			{ 0, 1. / 3, 2. / 3 }, // 50
			{ 0, 2. / 3, 2. / 3 }, // 51
			{ 1. / 3, 1. / 3, 1 }, // 52
			{ 2. / 3, 1. / 3, 1 }, // 53
			{ 2. / 3, 2. / 3, 1 }, // 54
			{ 1. / 3, 2. / 3, 1 }, // 55
			{ 1. / 3, 1. / 3, 1. / 3 }, // 56
			{ 2. / 3, 1. / 3, 1. / 3 }, // 57
			{ 2. / 3, 2. / 3, 1. / 3 }, // 58
			{ 1. / 3, 2. / 3, 1. / 3 }, // 59
			{ 1. / 3, 1. / 3, 2. / 3 }, // 60
			{ 2. / 3, 1. / 3, 2. / 3 }, // 61
			{ 2. / 3, 2. / 3, 2. / 3 }, // 62
			{ 1. / 3, 2. / 3, 2. / 3 } // 63
	};

	double xyz[64][3] = { 0 };

	// Now the nodal values, by simple functional evaluation.
	for (int ii = 0; ii < 64; ii++) {
		cubicXYZ(uvw[ii], xyz[ii]);
	}
	LCHM.setNodalValues(xyz);
	LCHM.setModalValues();

	double testUVW[] = { 1. / M_PI, 1. / M_E, 0.5 * (1 - 1. / M_PI - 1. / M_E) };
	double LCHMxyz[3], funcxyz[3];
	cubicXYZ(testUVW, funcxyz);
	LCHM.computeTransformedCoords(testUVW, LCHMxyz);
	BOOST_CHECK_CLOSE(LCHMxyz[0], funcxyz[0], 1.e-8);
	BOOST_CHECK_CLOSE(LCHMxyz[1], funcxyz[1], 1.e-8);
	BOOST_CHECK_CLOSE(LCHMxyz[2], funcxyz[2], 1.e-8);

	for (int ii = 0; ii < 64; ii++) {
		LCHM.computeTransformedCoords(uvw[ii], LCHMxyz);
//			printf("%d\n", ii);
		BOOST_CHECK_CLOSE(LCHMxyz[0], xyz[ii][0], 1.e-8);
		BOOST_CHECK_CLOSE(LCHMxyz[1], xyz[ii][1], 1.e-8);
		BOOST_CHECK_CLOSE(LCHMxyz[2], xyz[ii][2], 1.e-8);
	}
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(MappingTests, MixedMeshFixture)

BOOST_AUTO_TEST_CASE(EdgeMappingUniform)
 {
	// Check the uniform length scale cases.
	makeLengthScaleUniform(pUM_In);

	TetDivider TD(pUM_Out, pUM_In, 4);

	// Compute point locations and compare to the analytic
	// result
	const emInt *const thisTet = pUM_In->getTetConn(0);
	TD.setupCoordMapping(thisTet);
	TD.divideEdges(vertsOnEdges);

	BOOST_CHECK_EQUAL(vertsOnEdges.size(), 6);
	for (auto thisEdgeData : vertsOnEdges) {
		EdgeVerts EV = thisEdgeData.second;
		emInt startInd = EV.m_verts[0];
		emInt vertInd = EV.m_verts[1];
		emInt endInd = EV.m_verts[4];
		BOOST_CHECK_CLOSE(pUM_Out->getX(vertInd),
				0.75 * pUM_Out->getX(startInd) + 0.25 * pUM_Out->getX(endInd),
				1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getY(vertInd),
				0.75 * pUM_Out->getY(startInd) + 0.25 * pUM_Out->getY(endInd),
				1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getZ(vertInd),
				0.75 * pUM_Out->getZ(startInd) + 0.25 * pUM_Out->getZ(endInd),
				1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(TetFaceMappingUniform) {
	printf("Tet face uniform mapping test\n");
	// Check the uniform length scale cases.
	makeLengthScaleUniform(pUM_In);

	TetDivider TD(pUM_Out, pUM_In, 4);
	// Compute point locations and compare to the analytic
	// result
	const emInt *const thisTet = pUM_In->getTetConn(0);
	TD.setupCoordMapping(thisTet);
	TD.divideEdges(vertsOnEdges);
	TD.divideFaces(vertsOnTris, vertsOnQuads);

	BOOST_CHECK_EQUAL(vertsOnTris.size(), 4);
	BOOST_CHECK_EQUAL(vertsOnQuads.size(), 0);

	for (auto thisTriData : vertsOnTris) {
		emInt corners[] = { thisTriData.getCorner(0), thisTriData.getCorner(1),
				thisTriData.getCorner(2) };
		double coords[3][3];
		for (int ii = 0; ii < 3; ii++) {
			coords[ii][0] = pUM_Out->getX(corners[ii]);
			coords[ii][1] = pUM_Out->getY(corners[ii]);
			coords[ii][2] = pUM_Out->getZ(corners[ii]);
		}
		emInt vertInd = thisTriData.getIntVertInd(2, 1);
		BOOST_CHECK_CLOSE(pUM_Out->getX(vertInd),
				0.5 * coords[1][0] + 0.25 * coords[2][0] + 0.25 * coords[0][0],
				1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getY(vertInd),
				0.5 * coords[1][1] + 0.25 * coords[2][1] + 0.25 * coords[0][1],
				1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getZ(vertInd),
				0.5 * coords[1][2] + 0.25 * coords[2][2] + 0.25 * coords[0][2],
				1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(PyramidFaceMappingUniform) {
	printf("Pyramid face uniform mapping test\n");
	// Check the uniform length scale cases.
makeLengthScaleUniform(pUM_In);

	PyrDivider PD(pUM_Out, pUM_In, 4);
	// Compute point locations and compare to the analytic
	// result
	const emInt *const thisPyr = pUM_In->getPyrConn(0);
	PD.setupCoordMapping(thisPyr);
	PD.divideEdges(vertsOnEdges);
	PD.divideFaces(vertsOnTris, vertsOnQuads);

	BOOST_CHECK_EQUAL(vertsOnTris.size(), 4);
	BOOST_CHECK_EQUAL(vertsOnQuads.size(), 1);

	for (auto thisTriData : vertsOnTris) {
		emInt corners[] = { thisTriData.getCorner(0), thisTriData.getCorner(1),
				thisTriData.getCorner(2) };
		double coords[3][3];
		for (int ii = 0; ii < 3; ii++) {
			coords[ii][0] = pUM_Out->getX(corners[ii]);
			coords[ii][1] = pUM_Out->getY(corners[ii]);
			coords[ii][2] = pUM_Out->getZ(corners[ii]);
		}
		emInt vertInd = thisTriData.getIntVertInd(2,1);
		BOOST_CHECK_CLOSE(pUM_Out->getX(vertInd),
				0.5 * coords[1][0] + 0.25 * coords[2][0] + 0.25 * coords[0][0],
				1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getY(vertInd),
				0.5 * coords[1][1] + 0.25 * coords[2][1] + 0.25 * coords[0][1],
				1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getZ(vertInd),
				0.5 * coords[1][2] + 0.25 * coords[2][2] + 0.25 * coords[0][2],
				1.e-8);
	}

	for (auto thisQuadData : vertsOnQuads) {
		emInt corners[] = { thisQuadData.getCorner(0), thisQuadData.getCorner(1),
				thisQuadData.getCorner(2), thisQuadData.getCorner(3) };
		double coords[4][3];
		for (int ii = 0; ii < 4; ii++) {
			coords[ii][0] = pUM_Out->getX(corners[ii]);
			coords[ii][1] = pUM_Out->getY(corners[ii]);
			coords[ii][2] = pUM_Out->getZ(corners[ii]);
		}
		emInt vertInd = thisQuadData.getIntVertInd(2, 1);
		BOOST_CHECK_CLOSE(pUM_Out->getX(vertInd),
				0.375 * (coords[0][0] + coords[1][0])
						+ 0.125 * (coords[2][0] + coords[3][0]), 1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getY(vertInd),
				0.375 * (coords[0][1] + coords[1][1])
						+ 0.125 * (coords[2][1] + coords[3][1]), 1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getZ(vertInd),
				0.375 * (coords[0][2] + coords[1][2])
						+ 0.125 * (coords[2][2] + coords[3][2]), 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(PrismFaceMappingUniform) {
	printf("Prism face uniform mapping test\n");
	// Check the uniform length scale cases.
makeLengthScaleUniform(pUM_In);

	PrismDivider PD(pUM_Out, pUM_In, 4);
	// Compute point locations and compare to the analytic
	// result
	const emInt *const thisPrism = pUM_In->getPrismConn(0);
	PD.setupCoordMapping(thisPrism);
	PD.divideEdges(vertsOnEdges);
	PD.divideFaces(vertsOnTris, vertsOnQuads);

	BOOST_CHECK_EQUAL(vertsOnTris.size(), 2);
	BOOST_CHECK_EQUAL(vertsOnQuads.size(), 3);

	for (auto thisTriData : vertsOnTris) {
		emInt corners[] = { thisTriData.getCorner(0), thisTriData.getCorner(1),
				thisTriData.getCorner(2) };
		double coords[3][3];
		for (int ii = 0; ii < 3; ii++) {
			coords[ii][0] = pUM_Out->getX(corners[ii]);
			coords[ii][1] = pUM_Out->getY(corners[ii]);
			coords[ii][2] = pUM_Out->getZ(corners[ii]);
		}
		emInt vertInd = thisTriData.getIntVertInd(2, 1);
		BOOST_CHECK_CLOSE(pUM_Out->getX(vertInd),
				0.5 * coords[1][0] + 0.25 * coords[2][0] + 0.25 * coords[0][0],
				1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getY(vertInd),
				0.5 * coords[1][1] + 0.25 * coords[2][1] + 0.25 * coords[0][1],
				1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getZ(vertInd),
				0.5 * coords[1][2] + 0.25 * coords[2][2] + 0.25 * coords[0][2],
				1.e-8);
	}

	for (auto thisQuadData : vertsOnQuads) {
		emInt corners[] = { thisQuadData.getCorner(0), thisQuadData.getCorner(1),
				thisQuadData.getCorner(2), thisQuadData.getCorner(3) };
		double coords[4][3];
		for (int ii = 0; ii < 4; ii++) {
			coords[ii][0] = pUM_Out->getX(corners[ii]);
			coords[ii][1] = pUM_Out->getY(corners[ii]);
			coords[ii][2] = pUM_Out->getZ(corners[ii]);
		}
		emInt vertInd = thisQuadData.getIntVertInd(2, 1);
		BOOST_CHECK_CLOSE(pUM_Out->getX(vertInd),
				0.375 * (coords[0][0] + coords[1][0])
						+ 0.125 * (coords[2][0] + coords[3][0]), 1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getY(vertInd),
				0.375 * (coords[0][1] + coords[1][1])
						+ 0.125 * (coords[2][1] + coords[3][1]), 1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getZ(vertInd),
				0.375 * (coords[0][2] + coords[1][2])
						+ 0.125 * (coords[2][2] + coords[3][2]), 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(HexFaceMappingUniform) {
	printf("Hex face uniform mapping test\n");
	// Check the uniform length scale cases.
makeLengthScaleUniform(pUM_In);

	HexDivider HD(pUM_Out, pUM_In, 4);
	// Compute point locations and compare to the analytic
	// result
	const emInt *const thisHex = pUM_In->getHexConn(0);
	HD.setupCoordMapping(thisHex);
	HD.divideEdges(vertsOnEdges);
	HD.divideFaces(vertsOnTris, vertsOnQuads);

	BOOST_CHECK_EQUAL(vertsOnTris.size(), 0);
	BOOST_CHECK_EQUAL(vertsOnQuads.size(), 6);

	for (auto thisQuadData : vertsOnQuads) {
		emInt corners[] = { thisQuadData.getCorner(0), thisQuadData.getCorner(1),
				thisQuadData.getCorner(2), thisQuadData.getCorner(3) };
		double coords[4][3];
		for (int ii = 0; ii < 4; ii++) {
			coords[ii][0] = pUM_Out->getX(corners[ii]);
			coords[ii][1] = pUM_Out->getY(corners[ii]);
			coords[ii][2] = pUM_Out->getZ(corners[ii]);
		}
		emInt vertInd = thisQuadData.getIntVertInd(2, 1);
		BOOST_CHECK_CLOSE(pUM_Out->getX(vertInd),
				0.375 * (coords[0][0] + coords[1][0])
						+ 0.125 * (coords[2][0] + coords[3][0]), 1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getY(vertInd),
				0.375 * (coords[0][1] + coords[1][1])
						+ 0.125 * (coords[2][1] + coords[3][1]), 1.e-8);
		BOOST_CHECK_CLOSE(pUM_Out->getZ(vertInd),
				0.375 * (coords[0][2] + coords[1][2])
						+ 0.125 * (coords[2][2] + coords[3][2]), 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(TetMappingUniform) {
	printf("Tet uniform mapping\n");
	// Check the uniform length scale cases.
makeLengthScaleUniform(pUM_In);

	TetDivider TD(pUM_Out, pUM_In, 5);

	const emInt *const thisTet = pUM_In->getTetConn(0);
	TD.setupCoordMapping(thisTet);

	TD.createDivisionVerts(vertsOnEdges, vertsOnTris, vertsOnQuads);

	double uvw[4][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
	double xyz[4][3];
	TD.getPhysCoordsFromParamCoords(uvw[0], xyz[0]);
	TD.getPhysCoordsFromParamCoords(uvw[1], xyz[1]);
	TD.getPhysCoordsFromParamCoords(uvw[2], xyz[2]);
	TD.getPhysCoordsFromParamCoords(uvw[3], xyz[3]);
	BOOST_CHECK_CLOSE(xyz[0][0], 0, 1.e-8);
	BOOST_CHECK_CLOSE(xyz[0][1], -1, 1.e-8);
	BOOST_CHECK_CLOSE(xyz[0][2], 0, 1.e-8);

	BOOST_CHECK_CLOSE(xyz[1][0], 1, 1.e-8);
	BOOST_CHECK_CLOSE(xyz[1][1], 0, 1.e-8);
	BOOST_CHECK_CLOSE(xyz[1][2], 0, 1.e-8);

	BOOST_CHECK_CLOSE(xyz[2][0], 0, 1.e-8);
	BOOST_CHECK_CLOSE(xyz[2][1], 0, 1.e-8);
	BOOST_CHECK_CLOSE(xyz[2][2], 0, 1.e-8);

	BOOST_CHECK_CLOSE(xyz[3][0], 0, 1.e-8);
	BOOST_CHECK_CLOSE(xyz[3][1], 0, 1.e-8);
	BOOST_CHECK_CLOSE(xyz[3][2], 1, 1.e-8);

	// Check correctness of the uvw values.
	emInt iArray[] = { 1, 1, 1, 2 };
	emInt jArray[] = { 1, 1, 2, 1 };
	emInt kArray[] = { 2, 1, 1, 1 };
	double uArray[] = { 0.2, 0.2, 0.2, 0.4 };
	double vArray[] = { 0.2, 0.2, 0.4, 0.2 };
	double wArray[] = { 0.4, 0.2, 0.2, 0.2 };

	double xArray[] = { 0.2, 0.2, 0.2, 0.4 };
	double yArray[] = { -0.2, -0.4, -0.2, -0.2 };
	double zArray[] = { 0.4, 0.2, 0.2, 0.2 };

	// Now extract the uvw values and check correctness
	for (int ii = 0; ii < 4; ii++) {
		double my_uvw[3], my_xyz[3];
		TD.getParamCoords(iArray[ii], jArray[ii], kArray[ii], my_uvw);
		BOOST_CHECK_CLOSE(uArray[ii], my_uvw[0], 1.e-8);
		BOOST_CHECK_CLOSE(vArray[ii], my_uvw[1], 1.e-8);
		BOOST_CHECK_CLOSE(wArray[ii], my_uvw[2], 1.e-8);

		TD.getPhysCoordsFromParamCoords(my_uvw, my_xyz);
		BOOST_CHECK_CLOSE(xArray[ii], my_xyz[0], 1.e-8);
		BOOST_CHECK_CLOSE(yArray[ii], my_xyz[1], 1.e-8);
		BOOST_CHECK_CLOSE(zArray[ii], my_xyz[2], 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(PyramidMappingUniform) {
	printf("Pyramid uniform mapping\n");
	// Check the uniform length scale cases.
makeLengthScaleUniform(pUM_In);

	PyrDivider PD(pUM_Out, pUM_In, 5);

	const emInt *const thisPyr = pUM_In->getPyrConn(0);
	PD.setupCoordMapping(thisPyr);

	PD.createDivisionVerts(vertsOnEdges, vertsOnTris, vertsOnQuads);

	double uvw[5][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, {1, 1, 0}, { 0, 0, 1 } };
	double xyz[5][3];
	PD.getPhysCoordsFromParamCoords(uvw[0], xyz[0]);
	PD.getPhysCoordsFromParamCoords(uvw[1], xyz[1]);
	PD.getPhysCoordsFromParamCoords(uvw[2], xyz[2]);
	PD.getPhysCoordsFromParamCoords(uvw[3], xyz[3]);
	PD.getPhysCoordsFromParamCoords(uvw[4], xyz[4]);

	// Check correctness of the uvw values.
	emInt iArray[] = { 1, 1, 1, 2, 2 };
	emInt jArray[] = { 1, 1, 2, 1, 2 };
	emInt kArray[] = { 2, 1, 1, 1, 1 };
	double uArray[] = { 0.2, 0.2, 0.2, 0.4, 0.4 };
	double vArray[] = { 0.2, 0.2, 0.4, 0.2, 0.4 };
	double wArray[] = { 0.4, 0.2, 0.2, 0.2, 0.2 };

	double xArray[] = { 0.2, 0.2, 0.2, 0.4, 0.4 };
	double yArray[] = { 0.2, 0.2, 0.4, 0.2, 0.4 };
	double zArray[] = { 0.4, 0.2, 0.2, 0.2, 0.2 };

	// Now extract the uvw values and check correctness
	for (int ii = 0; ii < 4; ii++) {
		double my_uvw[3], my_xyz[3];
		PD.getParamCoords(iArray[ii], jArray[ii], kArray[ii], my_uvw);
		BOOST_CHECK_CLOSE(uArray[ii], my_uvw[0], 1.e-8);
		BOOST_CHECK_CLOSE(vArray[ii], my_uvw[1], 1.e-8);
		BOOST_CHECK_CLOSE(wArray[ii], my_uvw[2], 1.e-8);

		PD.getPhysCoordsFromParamCoords(my_uvw, my_xyz);
		BOOST_CHECK_CLOSE(xArray[ii], my_xyz[0], 1.e-8);
		BOOST_CHECK_CLOSE(yArray[ii], my_xyz[1], 1.e-8);
		BOOST_CHECK_CLOSE(zArray[ii], my_xyz[2], 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(PrismMappingUniform) {
	printf("Prism uniform mapping\n");
	// Check the uniform length scale cases.
makeLengthScaleUniform(pUM_In);

	PrismDivider PD(pUM_Out, pUM_In, 4);

	const emInt *const thisPrism = pUM_In->getPrismConn(0);
	PD.setupCoordMapping(thisPrism);

	PD.createDivisionVerts(vertsOnEdges, vertsOnTris, vertsOnQuads);

	double uvw[6][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 },
			{ 1, 0, 1 }, { 0, 1, 1 }};
	double xyz[6][3];
	PD.getPhysCoordsFromParamCoords(uvw[0], xyz[0]);
	PD.getPhysCoordsFromParamCoords(uvw[1], xyz[1]);
	PD.getPhysCoordsFromParamCoords(uvw[2], xyz[2]);
	PD.getPhysCoordsFromParamCoords(uvw[3], xyz[3]);
	PD.getPhysCoordsFromParamCoords(uvw[4], xyz[4]);
	PD.getPhysCoordsFromParamCoords(uvw[5], xyz[5]);

	// Check correctness of the uvw values.
	emInt iArray[] = { 1, 1, 2, 1, 1, 2, 1, 1, 2 };
	emInt jArray[] = { 1, 2, 1, 1, 2, 1, 1, 2, 1 };
	emInt kArray[] = { 1, 1, 1, 2, 2, 2, 3, 3, 3 };
	double uArray[] = { 0.25, 0.25, 0.5, 0.25, 0.25, 0.5, 0.25, 0.25, 0.5 };
	double vArray[] = { 0.25, 0.5, 0.25, 0.25, 0.5, 0.25, 0.25, 0.5, 0.25 };
	double wArray[] = { 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75 };

	double xArray[] = { 0.25, 0.25, 0.5, 0.25, 0.25, 0.5, 0.25, 0.25, 0.5 };
	double yArray[] = { -0.5, -0.25, -0.25, -0.5, -0.25, -0.25, -0.5, -0.25, -0.25 };
	double zArray[] = { -0.75, -0.75, -0.75, -0.5, -0.5, -0.5, -0.25, -0.25, -0.25 };

	// Now extract the uvw values and check correctness
	for (int ii = 0; ii < 9; ii++) {
		double my_uvw[3], my_xyz[3];
		PD.getParamCoords(iArray[ii], jArray[ii], kArray[ii], my_uvw);
		BOOST_CHECK_CLOSE(uArray[ii], my_uvw[0], 1.e-8);
		BOOST_CHECK_CLOSE(vArray[ii], my_uvw[1], 1.e-8);
		BOOST_CHECK_CLOSE(wArray[ii], my_uvw[2], 1.e-8);

		PD.getPhysCoordsFromParamCoords(my_uvw, my_xyz);
		BOOST_CHECK_CLOSE(xArray[ii], my_xyz[0], 1.e-8);
		BOOST_CHECK_CLOSE(yArray[ii], my_xyz[1], 1.e-8);
		BOOST_CHECK_CLOSE(zArray[ii], my_xyz[2], 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(HexMappingUniform) {
	printf("Hex uniform mapping\n");
	// Check the uniform length scale cases.
	makeLengthScaleUniform(pUM_In);

	HexDivider HD(pUM_Out, pUM_In, 3);

	const emInt *const thisHex = pUM_In->getHexConn(0);
	HD.setupCoordMapping(thisHex);

	HD.createDivisionVerts(vertsOnEdges, vertsOnTris, vertsOnQuads);

	double uvw[8][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 },
			{ 0, 0, 1 }, { 1, 0, 1 }, { 0, 1, 1 }, { 1, 1, 1 }	};
	double xyz[8][3];
	HD.getPhysCoordsFromParamCoords(uvw[0], xyz[0]);
	HD.getPhysCoordsFromParamCoords(uvw[1], xyz[1]);
	HD.getPhysCoordsFromParamCoords(uvw[2], xyz[2]);
	HD.getPhysCoordsFromParamCoords(uvw[3], xyz[3]);
	HD.getPhysCoordsFromParamCoords(uvw[4], xyz[4]);
	HD.getPhysCoordsFromParamCoords(uvw[5], xyz[5]);
	HD.getPhysCoordsFromParamCoords(uvw[6], xyz[6]);
	HD.getPhysCoordsFromParamCoords(uvw[7], xyz[7]);

	// Check correctness of the uvw values.
	emInt iArray[] = { 1, 1, 2, 2, 1, 1, 2, 2 };
	emInt jArray[] = { 1, 2, 2, 1, 1, 2, 2, 1 };
	emInt kArray[] = { 1, 1, 1, 1, 2, 2, 2, 2 };

	const double oneThird = 1./3;

	// Now extract the uvw values and check correctness
	for (int ii = 0; ii < 8; ii++) {
		double my_uvw[3], my_xyz[3];
		HD.getParamCoords(iArray[ii], jArray[ii], kArray[ii], my_uvw);
		BOOST_CHECK_CLOSE(iArray[ii]*oneThird, my_uvw[0], 1.e-8);
		BOOST_CHECK_CLOSE(jArray[ii]*oneThird, my_uvw[1], 1.e-8);
		BOOST_CHECK_CLOSE(kArray[ii]*oneThird, my_uvw[2], 1.e-8);

		HD.getPhysCoordsFromParamCoords(my_uvw, my_xyz);
		BOOST_CHECK_CLOSE(my_uvw[0], my_xyz[0], 1.e-8);
		BOOST_CHECK_CLOSE(my_uvw[1], my_xyz[1], 1.e-8);
		BOOST_CHECK_CLOSE(-1+my_uvw[2], my_xyz[2], 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(EdgeMappingNonuniformPrescribed) {
	printf("Edge non-uniform mapping\n");
	setPrescribedLengthScale(pUM_In);

	TetDivider TD(pUM_Out, pUM_In, 4);

	// Compute point locations and compare to the analytic
	// result
	const emInt *const thisTet = pUM_In->getTetConn(0);
	TD.setupCoordMapping(thisTet);
	TD.divideEdges(vertsOnEdges);

	BOOST_CHECK_EQUAL(vertsOnEdges.size(), 6);
	for (auto thisEdgeData : vertsOnEdges) {
		EdgeVerts EV = thisEdgeData.second;
		emInt startInd = EV.m_verts[0];
		emInt vertInd = EV.m_verts[1];
		emInt endInd = EV.m_verts[4];
		double startLenOrig = pUM_In->getLengthScale(startInd);
		double endLenOrig = pUM_In->getLengthScale(endInd);

		double startCoords[3], endCoords[3];
		pUM_In->getCoords(startInd, startCoords);
		pUM_In->getCoords(endInd, endCoords);
		double delta[] = {endCoords[0] - startCoords[0],
				endCoords[1] - startCoords[1], endCoords[2] - startCoords[2]};

		// Working out parametric coordinate along the edge, so (see Carl's
		// notes for June 16, 2020, or the journal article) x_A = 0, dx = 1,
		// and:
		//
		// u = startLen * xi + (3 - 2*startLen - endLen) * xi^2
		//     + (startLen + endLen - 2) * xi^3

		double startLen = sqrt(startLenOrig/endLenOrig);
		double endLen = 1. / startLen;

		double xi = 0.25;
		double u = startLen * xi + (3 - 2*startLen - endLen) * xi*xi
				   + (startLen + endLen - 2) * xi*xi*xi;

		// Adding 10 to avoid problems comparing machine zero to exactly zero.
		BOOST_CHECK_CLOSE(10+pUM_Out->getX(vertInd), 10+startCoords[0] + u*delta[0],
				1.e-10);
		BOOST_CHECK_CLOSE(10+pUM_Out->getY(vertInd), 10+startCoords[1] + u*delta[1],
				1.e-10);
		BOOST_CHECK_CLOSE(10+pUM_Out->getZ(vertInd), 10+startCoords[2] + u*delta[2],
				1.e-10);
	}
}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(NonuniformParamMapping)

BOOST_AUTO_TEST_CASE(FaceInteriorIntersectionTrivial) {
	// Trivial: coordinate aligned
	double uvL[] = {0, 0.25};
	double uvR[] = {1, 0.25};

	double uvB[] = {0.43, 0};
	double uvT[] = {0.43, 1};

	double uv[] = {0, 0};
	double uvExp[] = {0.43, 0.25};

	getFaceParametricIntersectionPoint(uvL, uvR, uvB, uvT, uv);

	BOOST_CHECK_CLOSE(uv[0], uvExp[0], 1.e-8);
	BOOST_CHECK_CLOSE(uv[1], uvExp[1], 1.e-8);
}

BOOST_AUTO_TEST_CASE(FaceInteriorIntersectionQuad)
{
	double uvL[] = {0, 0.207};
	double uvR[] = {1, 0.307};

	double uvB[] = {0.38, 0};
	double uvT[] = {0.58, 1};

	double uv[] = {0, 0};
	double uvExp[] = {0.43, 0.25};

	getFaceParametricIntersectionPoint(uvL, uvR, uvB, uvT, uv);

	BOOST_CHECK_CLOSE(uv[0], uvExp[0], 1.e-8);
	BOOST_CHECK_CLOSE(uv[1], uvExp[1], 1.e-8);
}

BOOST_AUTO_TEST_CASE(FaceInteriorIntersectionTri) {
	// Simulated triangle data, same result
	double uvL[] = {0, 0.36};
	double uvR[] = {0.86, 0.14};

	double uvB[] = {0.36, 0};
	double uvT[] = {0.5, 0.5};

	double uv[] = {0, 0};
	double uvExp[] = {0.43, 0.25};

	getFaceParametricIntersectionPoint(uvL, uvR, uvB, uvT, uv);

	BOOST_CHECK_CLOSE(uv[0], uvExp[0], 1.e-8);
	BOOST_CHECK_CLOSE(uv[1], uvExp[1], 1.e-8);
}

BOOST_AUTO_TEST_CASE(CellInteriorIntersectionTrivial) {
	// Trivial case: coordinate aligned, intersection
	double uvwXLo[] = {0, 0.25, 0.45};
	double uvwXHi[] = {1, 0.25, 0.45};
	double uvwYLo[] = {0.35, 0, 0.45};
	double uvwYHi[] = {0.35, 1, 0.45};
	double uvwZLo[] = {0.35, 0.25, 0};
	double uvwZHi[] = {0.35, 0.25, 1};
	double uvwExp[] = {0.35, 0.25, 0.45};
	double uvw[] = {0,0,0};
	getCellInteriorParametricIntersectionPoint(uvwXLo, uvwXHi, uvwYLo, uvwYHi, uvwZLo, uvwZHi, uvw);
	BOOST_CHECK_CLOSE(uvw[0], uvwExp[0], 1.e-8);
	BOOST_CHECK_CLOSE(uvw[1], uvwExp[1], 1.e-8);
	BOOST_CHECK_CLOSE(uvw[2], uvwExp[2], 1.e-8);
}

BOOST_AUTO_TEST_CASE(CellInteriorIntersectionHexExact) {
	// Harder but has the same (analytically exact) answer
	double uvwXLo[] = {0, 0.215, 0.485};
	double uvwXHi[] = {1, 0.315, 0.385};
	double uvwYLo[] = {0.325, 0, 0.475};
	double uvwYHi[] = {0.425, 1, 0.375};
	double uvwZLo[] = {0.305, 0.295, 0};
	double uvwZHi[] = {0.405, 0.195, 1};

	// s1 = 0.35  uvwX[] = {s1, 0.25+(s1-0.35)/10, 0.45-(s1-0.35)/10}
	// s2 = 0.25  uvwY[] = {0.35+(s2-0.25)/10, s2, 0.45-(s2-0.25)/10}
	// s3 = 0.45  uvwZ[] = {0.35+(s3-0.45)/10, 0.25-(s3-0.45)/10, s3}

	double uvwExp[] = {0.35, 0.25, 0.45};
	double uvw[] = {0,0,0};
	getCellInteriorParametricIntersectionPoint(uvwXLo, uvwXHi, uvwYLo, uvwYHi, uvwZLo, uvwZHi, uvw);
	BOOST_CHECK_CLOSE(uvw[0], uvwExp[0], 1.e-8);
	BOOST_CHECK_CLOSE(uvw[1], uvwExp[1], 1.e-8);
	BOOST_CHECK_CLOSE(uvw[2], uvwExp[2], 1.e-8);
}

BOOST_AUTO_TEST_CASE(CellInteriorIntersectionTetExact) {
	// Harder but has the same (analytically exact) answer
	double uvwXLo[] = {0, 0.235, 0.465};
	double uvwXHi[] = {0.3, 0.265, 0.435};
	double uvwYLo[] = {0.175, 0, 0.475};
	double uvwYHi[] = {0.135, 0.4, 0.435};
	double uvwZLo[] = {0.105, 0.295, 0};
	double uvwZHi[] = {0.165, 0.235, 0.6};

	// s1 = 0.15  uvwX[] = {s1, 0.25+(s1-0.15)/10, 0.45-(s1-0.15)/10}
	// s2 = 0.25  uvwY[] = {0.15+(s2-0.25)/10, s2, 0.45-(s2-0.25)/10}
	// s3 = 0.45  uvwZ[] = {0.15+(s3-0.45)/10, 0.25-(s3-0.45)/10, s3}

	double uvwExp[] = {0.15, 0.25, 0.45};
	double uvw[] = {0,0,0};
	getCellInteriorParametricIntersectionPoint(uvwXLo, uvwXHi, uvwYLo, uvwYHi, uvwZLo, uvwZHi, uvw);
	BOOST_CHECK_CLOSE(uvw[0], uvwExp[0], 1.e-8);
	BOOST_CHECK_CLOSE(uvw[1], uvwExp[1], 1.e-8);
	BOOST_CHECK_CLOSE(uvw[2], uvwExp[2], 1.e-8);
}

BOOST_AUTO_TEST_SUITE_END()


#ifdef DO_SUBDIVISION_TESTS
BOOST_AUTO_TEST_SUITE(SimpleMeshSubdivisionTests)

BOOST_AUTO_TEST_CASE(SizeTestSingleTetBy2) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 4;
	MSIn.nVerts = 4;
	MSIn.nBdryTris = 4;
	MSIn.nBdryQuads = 0;
	MSIn.nTets = 1;
	MSIn.nPyrs = 0;
	MSIn.nPrisms = 0;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 2, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 10);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 16);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 0);
	BOOST_CHECK_EQUAL(MSOut.nTets, 8);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 0);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 0);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 0);
}

BOOST_AUTO_TEST_CASE(SizeTestSingleTetBy3) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 4;
	MSIn.nVerts = 4;
	MSIn.nBdryTris = 4;
	MSIn.nBdryQuads = 0;
	MSIn.nTets = 1;
	MSIn.nPyrs = 0;
	MSIn.nPrisms = 0;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 3, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 20);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 36);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 0);
	BOOST_CHECK_EQUAL(MSOut.nTets, 27);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 0);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 0);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 0);
}

BOOST_AUTO_TEST_CASE(SizeTestSingleTetBy4) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 4;
	MSIn.nVerts = 4;
	MSIn.nBdryTris = 4;
	MSIn.nBdryQuads = 0;
	MSIn.nTets = 1;
	MSIn.nPyrs = 0;
	MSIn.nPrisms = 0;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 4, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 35);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 64);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 0);
	BOOST_CHECK_EQUAL(MSOut.nTets, 64);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 0);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 0);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 0);
}

BOOST_AUTO_TEST_CASE(SizeTestSinglePyrBy2) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 5;
	MSIn.nVerts = 5;
	MSIn.nBdryTris = 4;
	MSIn.nBdryQuads = 1;
	MSIn.nTets = 0;
	MSIn.nPyrs = 1;
	MSIn.nPrisms = 0;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 2, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 14);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 16);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 4);
	BOOST_CHECK_EQUAL(MSOut.nTets, 4);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 6);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 0);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 0);
}

BOOST_AUTO_TEST_CASE(SizeTestSinglePyrBy3) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 5;
	MSIn.nVerts = 5;
	MSIn.nBdryTris = 4;
	MSIn.nBdryQuads = 1;
	MSIn.nTets = 0;
	MSIn.nPyrs = 1;
	MSIn.nPrisms = 0;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 3, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 30);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 36);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 9);
	BOOST_CHECK_EQUAL(MSOut.nTets, 16);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 19);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 0);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 0);
}

BOOST_AUTO_TEST_CASE(SizeTestSinglePyrBy4) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 5;
	MSIn.nVerts = 5;
	MSIn.nBdryTris = 4;
	MSIn.nBdryQuads = 1;
	MSIn.nTets = 0;
	MSIn.nPyrs = 1;
	MSIn.nPrisms = 0;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 4, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 55);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 64);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 16);
	BOOST_CHECK_EQUAL(MSOut.nTets, 40);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 44);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 0);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 0);
}

BOOST_AUTO_TEST_CASE(SizeTestSinglePrismBy2) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 6;
	MSIn.nVerts = 6;
	MSIn.nBdryTris = 2;
	MSIn.nBdryQuads = 3;
	MSIn.nTets = 0;
	MSIn.nPyrs = 0;
	MSIn.nPrisms = 1;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 2, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 18);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 8);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 12);
	BOOST_CHECK_EQUAL(MSOut.nTets, 0);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 0);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 8);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 0);
}

BOOST_AUTO_TEST_CASE(SizeTestSinglePrismBy3) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 6;
	MSIn.nVerts = 6;
	MSIn.nBdryTris = 2;
	MSIn.nBdryQuads = 3;
	MSIn.nTets = 0;
	MSIn.nPyrs = 0;
	MSIn.nPrisms = 1;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 3, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 40);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 18);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 27);
	BOOST_CHECK_EQUAL(MSOut.nTets, 0);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 0);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 27);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 0);
}

BOOST_AUTO_TEST_CASE(SizeTestSinglePrismBy4) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 6;
	MSIn.nVerts = 6;
	MSIn.nBdryTris = 2;
	MSIn.nBdryQuads = 3;
	MSIn.nTets = 0;
	MSIn.nPyrs = 0;
	MSIn.nPrisms = 1;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 4, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 75);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 32);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 48);
	BOOST_CHECK_EQUAL(MSOut.nTets, 0);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 0);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 64);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 0);
}

BOOST_AUTO_TEST_CASE(SizeTestTwoPrismsBy5) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 8;
	MSIn.nVerts = 8;
	MSIn.nBdryTris = 4;
	MSIn.nBdryQuads = 4;
	MSIn.nTets = 0;
	MSIn.nPyrs = 0;
	MSIn.nPrisms = 2;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 5, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 216);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 100);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 100);
	BOOST_CHECK_EQUAL(MSOut.nTets, 0);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 0);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 250);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 0);
}

BOOST_AUTO_TEST_CASE(SizeTestSingleHexBy2) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 8;
	MSIn.nVerts = 8;
	MSIn.nBdryTris = 0;
	MSIn.nBdryQuads = 6;
	MSIn.nTets = 0;
	MSIn.nPyrs = 0;
	MSIn.nPrisms = 0;
	MSIn.nHexes = 1;
	computeMeshSize(MSIn, 2, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 27);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 0);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 24);
	BOOST_CHECK_EQUAL(MSOut.nTets, 0);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 0);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 0);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 8);
}

BOOST_AUTO_TEST_CASE(SizeTestSingleHexBy3) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 8;
	MSIn.nVerts = 8;
	MSIn.nBdryTris = 0;
	MSIn.nBdryQuads = 6;
	MSIn.nTets = 0;
	MSIn.nPyrs = 0;
	MSIn.nPrisms = 0;
	MSIn.nHexes = 1;
	computeMeshSize(MSIn, 3, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 64);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 0);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 54);
	BOOST_CHECK_EQUAL(MSOut.nTets, 0);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 0);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 0);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 27);
}

BOOST_AUTO_TEST_CASE(SizeTestSingleHexBy4) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 8;
	MSIn.nVerts = 8;
	MSIn.nBdryTris = 0;
	MSIn.nBdryQuads = 6;
	MSIn.nTets = 0;
	MSIn.nPyrs = 0;
	MSIn.nPrisms = 0;
	MSIn.nHexes = 1;
	computeMeshSize(MSIn, 4, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 125);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 0);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 96);
	BOOST_CHECK_EQUAL(MSOut.nTets, 0);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 0);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 0);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 64);
}

BOOST_AUTO_TEST_CASE(SizeTestMixedMeshBy5) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 182;
	MSIn.nVerts = 1093;
	MSIn.nBdryTris = 40;
	MSIn.nBdryQuads = 160;
	MSIn.nTets = 40;
	MSIn.nPyrs = 160;
	MSIn.nPrisms = 200;
	MSIn.nHexes = 800;
	computeMeshSize(MSIn, 5, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 122461);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 1000);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 4000);
	BOOST_CHECK_EQUAL(MSOut.nTets, 17800);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 13600);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 25000);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 100000);
}

BOOST_AUTO_TEST_CASE(SizeTestMixedMeshBy2) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 182;
	MSIn.nVerts = 1093;
	MSIn.nBdryTris = 40;
	MSIn.nBdryQuads = 160;
	MSIn.nTets = 40;
	MSIn.nPyrs = 160;
	MSIn.nPrisms = 200;
	MSIn.nHexes = 800;
	computeMeshSize(MSIn, 2, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 8125);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 160);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 640);
	BOOST_CHECK_EQUAL(MSOut.nTets, 960);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 960);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 1600);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 6400);
}

BOOST_AUTO_TEST_CASE(SizeTestMixedMeshBy3) {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 182;
	MSIn.nVerts = 1093;
	MSIn.nBdryTris = 40;
	MSIn.nBdryQuads = 160;
	MSIn.nTets = 40;
	MSIn.nPyrs = 160;
	MSIn.nPrisms = 200;
	MSIn.nHexes = 800;
	computeMeshSize(MSIn, 3, MSOut);
	BOOST_CHECK_EQUAL(MSOut.nVerts, 26857);
	BOOST_CHECK_EQUAL(MSOut.nBdryTris, 360);
	BOOST_CHECK_EQUAL(MSOut.nBdryQuads, 1440);
	BOOST_CHECK_EQUAL(MSOut.nTets, 3640);
	BOOST_CHECK_EQUAL(MSOut.nPyrs, 3040);
	BOOST_CHECK_EQUAL(MSOut.nPrisms, 5400);
	BOOST_CHECK_EQUAL(MSOut.nHexes, 21600);
}

BOOST_AUTO_TEST_CASE(SingleTetN2) {
	UMesh UM(4, 4, 4, 0, 1, 0, 0, 0);

	BOOST_CHECK_EQUAL(UM.maxNVerts(), 4);
	BOOST_CHECK_EQUAL(UM.maxNBdryTris(), 4);
	BOOST_CHECK_EQUAL(UM.maxNTets(), 1);

	double coords[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
	emInt triVerts[][3] = { { 0, 1, 2 }, { 0, 1, 3 }, { 1, 2, 3 }, { 2, 0, 3 } };
	emInt tetVerts[4] = { 0, 1, 2, 3 };

	UM.addVert(coords[0]);
	BOOST_CHECK_EQUAL(UM.numVerts(), 1);
	UM.addVert(coords[1]);
	UM.addVert(coords[2]);
	UM.addVert(coords[3]);
	BOOST_CHECK_EQUAL(UM.numVerts(), 4);
	BOOST_CHECK_EQUAL(UM.maxNVerts(), 4);

	UM.addBdryTri(triVerts[0]);
	BOOST_CHECK_EQUAL(UM.numBdryTris(), 1);
	UM.addBdryTri(triVerts[1]);
	UM.addBdryTri(triVerts[2]);
	UM.addBdryTri(triVerts[3]);
	BOOST_CHECK_EQUAL(UM.numBdryTris(), 4);
	BOOST_CHECK_EQUAL(UM.maxNBdryTris(), 4);

	UM.addTet(tetVerts);
	BOOST_CHECK_EQUAL(UM.numTets(), 1);
	BOOST_CHECK_EQUAL(UM.maxNTets(), 1);

	makeLengthScaleUniform(&UM);

	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 4;
	MSIn.nVerts = 4;
	MSIn.nBdryTris = 4;
	MSIn.nBdryQuads = 0;
	MSIn.nTets = 1;
	MSIn.nPyrs = 0;
	MSIn.nPrisms = 0;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 2, MSOut);

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris,
			MSOut.nBdryQuads, MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms,
			MSOut.nHexes);
	subdividePartMesh(&UM, &UMOut, 2);
	checkExpectedSize(UMOut);
}

BOOST_AUTO_TEST_CASE(SingleTetN5) {
	UMesh UM(4, 4, 4, 0, 1, 0, 0, 0);

	BOOST_CHECK_EQUAL(UM.maxNVerts(), 4);
	BOOST_CHECK_EQUAL(UM.maxNBdryTris(), 4);
	BOOST_CHECK_EQUAL(UM.maxNTets(), 1);

	double coords[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
	emInt triVerts[][3] = { { 0, 1, 2 }, { 1, 0, 3 }, { 2, 1, 3 }, { 0, 2, 3 } };
	emInt tetVerts[4] = { 0, 1, 2, 3 };

	UM.addVert(coords[0]);
	UM.addVert(coords[1]);
	UM.addVert(coords[2]);
	UM.addVert(coords[3]);

	UM.addBdryTri(triVerts[0]);
	UM.addBdryTri(triVerts[1]);
	UM.addBdryTri(triVerts[2]);
	UM.addBdryTri(triVerts[3]);

	UM.addTet(tetVerts);

	makeLengthScaleUniform(&UM);

	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 4;
	MSIn.nVerts = 4;
	MSIn.nBdryTris = 4;
	MSIn.nBdryQuads = 0;
	MSIn.nTets = 1;
	MSIn.nPyrs = 0;
	MSIn.nPrisms = 0;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 5, MSOut);

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris,
			MSOut.nBdryQuads, MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms,
			MSOut.nHexes);
	subdividePartMesh(&UM, &UMOut, 5);
	checkExpectedSize(UMOut);
}

BOOST_AUTO_TEST_CASE(SinglePyrN3) {
	UMesh UM(5, 5, 4, 1, 0, 1, 0, 0);
	double coords[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, {
			0, 0, 1 } };
	emInt triVerts[][3] = { { 0, 1, 4 }, { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 } };
	emInt quadVerts[4] = { 0, 1, 2, 3 };
	emInt pyrVerts[5] = { 0, 1, 2, 3, 4 };
	for (int ii = 0; ii < 5; ii++) {
		UM.addVert(coords[ii]);
	}
	for (int ii = 0; ii < 4; ii++) {
		UM.addBdryTri(triVerts[ii]);
	}
	UM.addBdryQuad(quadVerts);
	UM.addPyramid(pyrVerts);

	makeLengthScaleUniform(&UM);

	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 5;
	MSIn.nVerts = 5;
	MSIn.nBdryTris = 4;
	MSIn.nBdryQuads = 1;
	MSIn.nTets = 0;
	MSIn.nPyrs = 1;
	MSIn.nPrisms = 0;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 3, MSOut);

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris,
			MSOut.nBdryQuads, MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms,
			MSOut.nHexes);
	subdividePartMesh(&UM, &UMOut, 3);
	checkExpectedSize(UMOut);
}

BOOST_AUTO_TEST_CASE(SinglePrismN3) {
	UMesh UM(6, 6, 2, 3, 0, 0, 1, 0);
	double coords[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, {
			1, 0, 1 }, { 0, 1, 1 } };
	emInt triVerts[][3] = { { 0, 1, 2 }, { 3, 4, 5 } };
	emInt quadVerts[][4] = { { 0, 1, 4, 3 }, { 1, 2, 5, 4 }, { 2, 0, 3, 5 } };
	emInt prismVerts[6] = { 0, 1, 2, 3, 4, 5 };
	for (int ii = 0; ii < 6; ii++) {
		UM.addVert(coords[ii]);
	}
	for (int ii = 0; ii < 2; ii++) {
		UM.addBdryTri(triVerts[ii]);
	}
	for (int ii = 0; ii < 3; ii++) {
		UM.addBdryQuad(quadVerts[ii]);
	}
	UM.addPrism(prismVerts);

	makeLengthScaleUniform(&UM);

	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 6;
	MSIn.nVerts = 6;
	MSIn.nBdryTris = 2;
	MSIn.nBdryQuads = 3;
	MSIn.nTets = 0;
	MSIn.nPyrs = 0;
	MSIn.nPrisms = 1;
	MSIn.nHexes = 0;
	computeMeshSize(MSIn, 3, MSOut);

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris,
			MSOut.nBdryQuads, MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms,
			MSOut.nHexes);
	subdividePartMesh(&UM, &UMOut, 3);
	checkExpectedSize(UMOut);
}

BOOST_AUTO_TEST_CASE(MixedN2) {
	UMesh UM(11, 11, 6, 6, 1, 1, 1, 1);
	double coords[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, {
			0, 0, 1 }, { 0, 0, -1 }, { 1, 0, -1 }, { 1, 1, -1 }, { 0, 1, -1 }, {
			0, -1, 0 }, { 0, -1, -1 } };
	emInt triVerts[][3] = { { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 }, { 0, 9, 4 },
			{ 9, 1, 4 }, { 10, 6, 5 } };
	emInt quadVerts[][4] = { { 6, 7, 2, 1 }, { 7, 8, 3, 2 }, { 8, 5, 0, 3 }, {
			10, 6, 1, 9 }, { 5, 10, 9, 0 }, { 5, 6, 7, 8 } };
	emInt tetVerts[4] = { 9, 1, 0, 4 };
	emInt pyrVerts[5] = { 0, 1, 2, 3, 4 };
	emInt prismVerts[6] = { 10, 6, 5, 9, 1, 0 };
	emInt hexVerts[8] = { 5, 6, 7, 8, 0, 1, 2, 3 };

	for (int ii = 0; ii < 11; ii++) {
		UM.addVert(coords[ii]);
	}
	for (int ii = 0; ii < 6; ii++) {
		UM.addBdryTri(triVerts[ii]);
		UM.addBdryQuad(quadVerts[ii]);
	}
	UM.addTet(tetVerts);
	UM.addPyramid(pyrVerts);
	UM.addPrism(prismVerts);
	UM.addHex(hexVerts);

	makeLengthScaleUniform(&UM);

	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 11;
	MSIn.nVerts = 11;
	MSIn.nBdryTris = 6;
	MSIn.nBdryQuads = 6;
	MSIn.nTets = 1;
	MSIn.nPyrs = 1;
	MSIn.nPrisms = 1;
	MSIn.nHexes = 1;
	computeMeshSize(MSIn, 2, MSOut);

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris,
			MSOut.nBdryQuads, MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms,
			MSOut.nHexes);
	subdividePartMesh(&UM, &UMOut, 2);
	checkExpectedSize(UMOut);
}

BOOST_AUTO_TEST_CASE(MixedN3) {
	UMesh UM(11, 11, 6, 6, 1, 1, 1, 1);
	double coords[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, {
			0, 0, 1 }, { 0, 0, -1 }, { 1, 0, -1 }, { 1, 1, -1 }, { 0, 1, -1 }, {
			0, -1, 0 }, { 0, -1, -1 } };
	emInt triVerts[][3] = { { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 }, { 0, 9, 4 },
			{ 9, 1, 4 }, { 10, 6, 5 } };
	emInt quadVerts[][4] = { { 6, 7, 2, 1 }, { 7, 8, 3, 2 }, { 8, 5, 0, 3 }, {
			10, 6, 1, 9 }, { 5, 10, 9, 0 }, { 5, 6, 7, 8 } };
	emInt tetVerts[4] = { 9, 1, 0, 4 };
	emInt pyrVerts[5] = { 0, 1, 2, 3, 4 };
	emInt prismVerts[6] = { 10, 6, 5, 9, 1, 0 };
	emInt hexVerts[8] = { 5, 6, 7, 8, 0, 1, 2, 3 };

	for (int ii = 0; ii < 11; ii++) {
		UM.addVert(coords[ii]);
	}
	for (int ii = 0; ii < 6; ii++) {
		UM.addBdryTri(triVerts[ii]);
		UM.addBdryQuad(quadVerts[ii]);
	}
	UM.addTet(tetVerts);
	UM.addPyramid(pyrVerts);
	UM.addPrism(prismVerts);
	UM.addHex(hexVerts);

	makeLengthScaleUniform(&UM);

	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 11;
	MSIn.nVerts = 11;
	MSIn.nBdryTris = 6;
	MSIn.nBdryQuads = 6;
	MSIn.nTets = 1;
	MSIn.nPyrs = 1;
	MSIn.nPrisms = 1;
	MSIn.nHexes = 1;
	computeMeshSize(MSIn, 3, MSOut);

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris,
			MSOut.nBdryQuads, MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms,
			MSOut.nHexes);
	subdividePartMesh(&UM, &UMOut, 3);
	checkExpectedSize(UMOut);
	bool result = UMOut.writeVTKFile("/tmp/test-exa.vtk");
	BOOST_CHECK(result);
	result = UMOut.writeUGridFile("/tmp/test-exa.b8.ugrid");
	BOOST_CHECK(result);
}

BOOST_AUTO_TEST_CASE(MixedN5) {
	UMesh UM(11, 11, 6, 6, 1, 1, 1, 1);
	double coords[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, {
			0, 0, 1 }, { 0, 0, -1 }, { 1, 0, -1 }, { 1, 1, -1 }, { 0, 1, -1 }, {
			0, -1, 0 }, { 0, -1, -1 } };
	emInt triVerts[][3] = { { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 }, { 0, 9, 4 },
			{ 9, 1, 4 }, { 10, 6, 5 } };
	emInt quadVerts[][4] = { { 6, 7, 2, 1 }, { 7, 8, 3, 2 }, { 8, 5, 0, 3 }, {
			10, 6, 1, 9 }, { 5, 10, 9, 0 }, { 5, 6, 7, 8 } };
	emInt tetVerts[4] = { 9, 1, 0, 4 };
	emInt pyrVerts[5] = { 0, 1, 2, 3, 4 };
	emInt prismVerts[6] = { 10, 6, 5, 9, 1, 0 };
	emInt hexVerts[8] = { 5, 6, 7, 8, 0, 1, 2, 3 };

	for (int ii = 0; ii < 11; ii++) {
		UM.addVert(coords[ii]);
	}
	for (int ii = 0; ii < 6; ii++) {
		UM.addBdryTri(triVerts[ii]);
		UM.addBdryQuad(quadVerts[ii]);
	}
	UM.addTet(tetVerts);
	UM.addPyramid(pyrVerts);
	UM.addPrism(prismVerts);
	UM.addHex(hexVerts);

	makeLengthScaleUniform(&UM);

	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = 11;
	MSIn.nVerts = 11;
	MSIn.nBdryTris = 6;
	MSIn.nBdryQuads = 6;
	MSIn.nTets = 1;
	MSIn.nPyrs = 1;
	MSIn.nPrisms = 1;
	MSIn.nHexes = 1;
	computeMeshSize(MSIn, 5, MSOut);

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris,
			MSOut.nBdryQuads, MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms,
			MSOut.nHexes);
	subdividePartMesh(&UM, &UMOut, 5);
	checkExpectedSize(UMOut);
	bool result = UMOut.writeVTKFile("/tmp/test-exa.vtk");
	BOOST_CHECK(result);
}
BOOST_FIXTURE_TEST_SUITE(FaceMatching,MixedMeshFixture)
BOOST_AUTO_TEST_CASE(QuadMatching){
	// run for quad cases of -2 , 1 , -4 
	
	const emInt nDivs=3; 
	const emInt partID=0 ; 
	const emInt remoteID=1; 

	emInt globalForRef_LocalQuad [4] = {1,4,16,13};
	emInt remoteForRef [4] = {10,20,30,40};
	emInt localForRef  [4] = {10,20,30,40}; 
	exa_set<QuadFaceVerts> remoteQuads;
	//Note that passing the same temp local & remote indices 
	QuadFaceVerts refQuad (nDivs,localForRef,globalForRef_LocalQuad,
	remoteForRef,partID,remoteID); 
	refQuad.setCompare(true); 
	SetArtificialIntVertQuad(refQuad,nDivs);
	

	//emInt globalForeRemote [4]= {13,16,4,1}; 
	emInt AllCasesofglobalForeRemote [7][4]={
		{13,16,4,1}, 
		{16,4,1,13},
		{4,1,13,16},
		{1,13,16,4},
		{13,1,4,16},
		{16,13,1,4},
		{4,16,13,1}
	}; 
	for(auto k=0; k<7; k++){
		emInt globalForeRemote[4];
		for (int i = 0; i < 4; i++) {
 			globalForeRemote[i] = AllCasesofglobalForeRemote[k][i];
		}; 
		
		QuadFaceVerts quad (nDivs,localForRef,globalForeRemote
		, remoteForRef,remoteID,partID); 
		quad.setCompare(true); 

		exa_set<QuadFaceVerts> setQuads={quad};
	
		emInt rotation= getQuadRotation(refQuad,setQuads,nDivs); 
		
		SetArtificialIntVertQuad(quad,nDivs);
		remoteQuads.insert(quad); 
		std::unordered_map<emInt,emInt> map; 
		matchQuad(refQuad,rotation,nDivs,remoteQuads,map);
		std::unordered_map<emInt,emInt> expectedMap; 
		setExpectedMapping(rotation,expectedMap); 
		//bool mapsEqual = boost::range::is_permutation(map,expectedMap);
		assert(map==expectedMap);  
		std::cout<<"Passed Test for Rotation case: "<<rotation<<std::endl; 
	}
	// Giving us different rotated Quads
	

		
	// }


}
BOOST_AUTO_TEST_SUITE_END()
BOOST_FIXTURE_TEST_SUITE(MPIFunctions,MixedMeshFixture)
BOOST_AUTO_TEST_CASE(CustomTypeRegisteration){

	boost::mpi::environment env; 
	boost::mpi::communicator world; 
	const size_t containerSize =100; 
	emInt localTri [3]; 
	emInt globalTri [3]; 
	emInt remoteTri [3]; 

	emInt localQuad [4]; 
	emInt globalQuad [4]; 
	emInt remoteQuad [4]; 

	emInt nDivs ; 
	emInt partId ; 
	emInt remoteId; 
	emInt type; 
	emInt elemInd;
	bool globalCompare; 

	double m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax;
	emInt m_first, m_last, m_nParts;

	emInt m_index, m_cellType;
	double m_coords[3];
	if(world.size()>1){
		if(world.rank()==0){

			std::vector<TriFaceVerts> dummytris;
			std::vector<QuadFaceVerts> dummyquads; 

			std::vector<Part> dummyPart; 
			std::vector<CellPartData> dummyCellPartData; 

			for(size_t i =0 ; i<containerSize ; i++){
				setArbitraryTriDataForTesting(i,localTri,globalTri,remoteTri,
				nDivs,partId,remoteId,type,elemInd,globalCompare); 

				setArbitraryQuadDataForTesting(i,localQuad,globalQuad,remoteQuad,
				nDivs,partId,remoteId,type,elemInd,globalCompare); 

				setArbitrary_Part_DataForTesting(i,m_xmin,m_xmax,
				m_ymin,m_ymax,m_zmin,m_zmax,m_first,m_last,m_nParts); 

				setArbitrary_CellPartData_ForTesting(i,m_coords,m_index,m_cellType); 

				Part part(m_first,m_last,m_nParts,m_xmin,
				m_xmax,m_ymin,m_ymax,m_zmin,m_zmax); 

				CellPartData cellPartData (m_index,m_cellType,m_coords[0],m_coords[1],
				m_coords[2]); 



				QuadFaceVerts quad(nDivs,localQuad,globalQuad,remoteQuad,partId,
				remoteId,type,elemInd,globalCompare); 


				TriFaceVerts tri (nDivs,localTri,globalTri,remoteTri,
				partId,remoteId,type,elemInd,globalCompare); 

				dummytris.push_back(tri); 
				dummyquads.push_back(quad); 
				dummyPart.push_back(part); 
				dummyCellPartData.push_back(cellPartData); 
			}
			world.send(1,0,dummytris); 
			world.send(1,0,dummyquads); 
			world.send(1,0,dummyPart); 
			world.send(1,0,dummyCellPartData); 

		}
		if(world.rank()==1){
			std::vector<TriFaceVerts> dummytris(containerSize,TriFaceVerts(1)); 
			std::vector<QuadFaceVerts> dummyquads(containerSize,QuadFaceVerts(1)); 

			std::vector<Part> dummyPart(containerSize,Part()); 
			std::vector<CellPartData> dummyCellPartData(containerSize,CellPartData()); 


			world.recv(0,0,dummytris);
			world.recv(0,0,dummyquads);
			world.recv(0,0,dummyPart); 
			world.recv(0,0,dummyCellPartData); 
			

			for(size_t i=0 ; i<containerSize ; i++){
				setArbitraryTriDataForTesting(i,localTri,globalTri,remoteTri,
				nDivs,partId,remoteId,type,elemInd,globalCompare); 

				setArbitraryQuadDataForTesting(i,localQuad,globalQuad,remoteQuad,
				nDivs,partId,remoteId,type,elemInd,globalCompare);


				setArbitrary_Part_DataForTesting(i,m_xmin,m_xmax,
				m_ymin,m_ymax,m_zmin,m_zmax,m_first,m_last,m_nParts); 

				setArbitrary_CellPartData_ForTesting(i,m_coords,m_index,m_cellType); 


				BOOST_CHECK_EQUAL(dummytris[i].getNumDivs(),nDivs); 
				BOOST_CHECK_EQUAL(dummytris[i].getPartid(),partId); 
				BOOST_CHECK_EQUAL(dummytris[i].getRemoteId(),remoteId); 
				BOOST_CHECK_EQUAL(dummytris[i].getVolElementType(),type); 
				BOOST_CHECK_EQUAL(dummytris[i].getVolElement(),elemInd); 
				BOOST_CHECK_EQUAL(dummytris[i].getGlobalCompare(),globalCompare); 


				BOOST_CHECK_EQUAL(dummyquads[i].getNumDivs(),nDivs); 
				BOOST_CHECK_EQUAL(dummyquads[i].getPartid(),partId); 
				BOOST_CHECK_EQUAL(dummyquads[i].getRemoteId(),remoteId); 
				BOOST_CHECK_EQUAL(dummyquads[i].getVolElementType(),type); 
				BOOST_CHECK_EQUAL(dummyquads[i].getVolElement(),elemInd); 
				BOOST_CHECK_EQUAL(dummyquads[i].getGlobalCompare(),globalCompare); 

				BOOST_CHECK_EQUAL(dummyPart[i].getFirst(),m_first); 
				BOOST_CHECK_EQUAL(dummyPart[i].getLast(), m_last); 
				BOOST_CHECK_EQUAL(dummyPart[i].getXmax(), m_xmax); 
				BOOST_CHECK_EQUAL(dummyPart[i].getXmin(), m_xmin); 
				BOOST_CHECK_EQUAL(dummyPart[i].getYmax(), m_ymax); 
				BOOST_CHECK_EQUAL(dummyPart[i].getYmin(), m_ymin); 
				BOOST_CHECK_EQUAL(dummyPart[i].getZmax(), m_zmax); 
				BOOST_CHECK_EQUAL(dummyPart[i].getZmin(), m_zmin); 

				BOOST_CHECK_EQUAL(dummyCellPartData[i].getCellType(),m_cellType); 
				BOOST_CHECK_EQUAL(dummyCellPartData[i].getCoord(0)  ,m_coords[0]); 
				BOOST_CHECK_EQUAL(dummyCellPartData[i].getCoord(1)  ,m_coords[1]); 
				BOOST_CHECK_EQUAL(dummyCellPartData[i].getCoord(2)  ,m_coords[2]); 


				for (auto k=0 ; k<4 ; k++){
					BOOST_CHECK_EQUAL(dummyquads[i].getGlobalCorner(k), 
					globalQuad[k]); 
					BOOST_CHECK_EQUAL(dummyquads[i].getCorner(k), 
					localQuad[k]); 
					BOOST_CHECK_EQUAL(dummyquads[i].getRemoteIndices(k), 
					remoteQuad[k]); 
				}
			}
		}
	}



}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE_END()


#endif // DO_SUBDIVISION_TESTS
