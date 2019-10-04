/*
 * test-exa.cxx
 *
 *  Created on: Jul. 5, 2019
 *      Author: cfog
 */

#define BOOST_TEST_MODULE test-exa
#include <boost/test/unit_test.hpp>

#include "examesh.h"
#include "UMesh.h"
#include "CubicMesh.h"

#include "TetDivider.h"

#include "Mapping.h"

static void checkExpectedSize(const UMesh& UM) {
	BOOST_CHECK_EQUAL(UM.maxNVerts(), UM.numVerts());
	BOOST_CHECK_EQUAL(UM.maxNBdryTris(), UM.numBdryTris());
	BOOST_CHECK_EQUAL(UM.maxNBdryQuads(), UM.numBdryQuads());
	BOOST_CHECK_EQUAL(UM.maxNTets(), UM.numTets());
	BOOST_CHECK_EQUAL(UM.maxNPyrs(), UM.numPyramids());
	BOOST_CHECK_EQUAL(UM.maxNPrisms(), UM.numPrisms());
	BOOST_CHECK_EQUAL(UM.maxNHexes(), UM.numHexes());
}

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

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris, MSOut.nBdryQuads,
							MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms, MSOut.nHexes);
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

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris, MSOut.nBdryQuads,
							MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms, MSOut.nHexes);
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

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris, MSOut.nBdryQuads,
							MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms, MSOut.nHexes);
	subdividePartMesh(&UM, &UMOut, 3);
	checkExpectedSize(UMOut);
}

BOOST_AUTO_TEST_CASE(SinglePrismN3) {
	UMesh UM(6, 6, 2, 3, 0, 0, 1, 0);
	double coords[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, {
			1, 0, 1 },
													{ 0, 1, 1 } };
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

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris, MSOut.nBdryQuads,
							MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms, MSOut.nHexes);
	subdividePartMesh(&UM, &UMOut, 3);
	checkExpectedSize(UMOut);
}

BOOST_AUTO_TEST_CASE(MixedN2) {
	UMesh UM(11, 11, 6, 6, 1, 1, 1, 1);
	double coords[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, {
			0, 0, 1 },
													{ 0, 0, -1 }, { 1, 0, -1 }, { 1, 1, -1 },
													{ 0, 1, -1 }, { 0, -1, 0 }, { 0, -1, -1 } };
	emInt triVerts[][3] = { { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 }, { 0, 9, 4 }, {
			9, 1, 4 },
													{ 10, 6, 5 } };
	emInt quadVerts[][4] = { { 6, 7, 2, 1 }, { 7, 8, 3, 2 }, { 8, 5, 0, 3 },
														{ 10, 6, 1, 9 }, { 5, 10, 9, 0 }, { 5, 6, 7, 8 } };
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

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris, MSOut.nBdryQuads,
							MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms, MSOut.nHexes);
	subdividePartMesh(&UM, &UMOut, 2);
	checkExpectedSize(UMOut);
}

BOOST_AUTO_TEST_CASE(MixedN3) {
	UMesh UM(11, 11, 6, 6, 1, 1, 1, 1);
	double coords[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, {
			0, 0, 1 },
													{ 0, 0, -1 }, { 1, 0, -1 }, { 1, 1, -1 },
													{ 0, 1, -1 }, { 0, -1, 0 }, { 0, -1, -1 } };
	emInt triVerts[][3] = { { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 }, { 0, 9, 4 }, {
			9, 1, 4 },
													{ 10, 6, 5 } };
	emInt quadVerts[][4] = { { 6, 7, 2, 1 }, { 7, 8, 3, 2 }, { 8, 5, 0, 3 },
														{ 10, 6, 1, 9 }, { 5, 10, 9, 0 }, { 5, 6, 7, 8 } };
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

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris, MSOut.nBdryQuads,
							MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms, MSOut.nHexes);
	subdividePartMesh(&UM, &UMOut, 3);
	checkExpectedSize(UMOut);
	bool result = UMOut.writeVTKFile("/tmp/test-exa.vtk");
	BOOST_CHECK(result);
}

BOOST_AUTO_TEST_CASE(MixedN5) {
	UMesh UM(11, 11, 6, 6, 1, 1, 1, 1);
	double coords[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, {
			0, 0, 1 },
													{ 0, 0, -1 }, { 1, 0, -1 }, { 1, 1, -1 },
													{ 0, 1, -1 }, { 0, -1, 0 }, { 0, -1, -1 } };
	emInt triVerts[][3] = { { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 }, { 0, 9, 4 }, {
			9, 1, 4 },
													{ 10, 6, 5 } };
	emInt quadVerts[][4] = { { 6, 7, 2, 1 }, { 7, 8, 3, 2 }, { 8, 5, 0, 3 },
														{ 10, 6, 1, 9 }, { 5, 10, 9, 0 }, { 5, 6, 7, 8 } };
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

	UMesh UMOut(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris, MSOut.nBdryQuads,
							MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms, MSOut.nHexes);
	subdividePartMesh(&UM, &UMOut, 5);
	checkExpectedSize(UMOut);
	bool result = UMOut.writeVTKFile("/tmp/test-exa.vtk");
	BOOST_CHECK(result);
}

BOOST_AUTO_TEST_SUITE(MappingTests)

	BOOST_AUTO_TEST_CASE(TetMapping) {
		printf("Testing tet mappings.\n");
		UMesh VM(4, 4, 4, 0, 1, 0, 0, 0);
		double vert0[] = { 1, 1, 1 };
		double vert1[] = { 2, 2, 2 };
		double vert2[] = { 1.5, 2, 2 };
		double vert3[] = { 1.5, 1, 3 };

		VM.addVert(vert0);
		VM.addVert(vert1);
		VM.addVert(vert2);
		VM.addVert(vert3);

		emInt verts[] = { 0, 1, 2, 3 };
		VM.addTet(verts);

		TetLengthScaleMapping TLSM(&VM);

		// Test for known cubic functions:
		// x = u^3 + 2 u^2 v + 4 u v^2 - 30 u v w + 5 u^2 + 6 u + 3
		// y = v^3 + 2 v^2 w + 4 v w^2 - 54 u v w + 5 v^2 + 6 v + 7
		// z = w^3 + 2 w^2 u + 4 w u^2 - 72 u v w + 5 w^2 + 6 w + 10
		// so
		// dx/du = 3 u^2 + 4 u v + 4 v^2 - 30 v w + 10 u + 6
		// dx/dv = 2 u^2 + 8 u v - 30 u w
		// dx/dw = -30 u v
		//
		// dy/du = -54 v w
		// dy/dv = 3 v^2 + 4 v w - 54 u w + 4 w^2 + 10 v + 6
		// dy/dw = 2 v^2 + 8 v w - 54 u v
		//
		// dz/du = 2 w^2 + 8 u w -72 v w
		// dz/dv = -54 u w
		// dz/dw = 3 w^2 + 4 u w + 4 u^2 + 10 w + 6 - 54 u v

		//
		// Now just evaluate these at points 0 (0,0,0), 1 (1,0,0), 2 (0,1,0),
		// and 3 (0,0,1) to get:
		{
			double xyz0[] = { 3, 7, 10 };
			double xyz1[] = { 15, 7, 10 };
			double xyz2[] = { 3, 19, 10 };
			double xyz3[] = { 3, 7, 22 };

			double uderiv0[] = { 6, 0, 0 };
			double vderiv0[] = { 0, 6, 0 };
			double wderiv0[] = { 0, 0, 6 };

			double uderiv1[] = { 19, 0, 0 };
			double vderiv1[] = { 2, 6, 0 };
			double wderiv1[] = { 0, 0, 10 };

			double uderiv2[] = { 10, 0, 0 };
			double vderiv2[] = { 0, 19, 0 };
			double wderiv2[] = { 0, 2, 6 };

			double uderiv3[] = { 6, 0, 2 };
			double vderiv3[] = { 0, 10, 0 };
			double wderiv3[] = { 0, 0, 19 };

			TLSM.setPolyCoeffs(xyz0, xyz1, xyz2, xyz3, uderiv0, vderiv0, wderiv0,
												uderiv1, vderiv1, wderiv1, uderiv2, vderiv2, wderiv2,
												uderiv3, vderiv3, wderiv3);

			double uvw[] = { 0.2, 0.3, 0.4 };
			double xyz[3];
			TLSM.computeTransformedCoords(uvw, xyz);
			BOOST_CHECK_CLOSE(xyz[0], 3.784, 1.e-8);
			BOOST_CHECK_CLOSE(xyz[1], 8.245, 1.e-8);
			BOOST_CHECK_CLOSE(xyz[2], 11.664, 1.e-8);
		}

		TLSM.setupCoordMapping(verts);

		double uvw[] = { 0, 0, 0 };
		double xyz[3];

		TLSM.computeTransformedCoords(uvw, xyz);
		BOOST_CHECK_CLOSE(xyz[0], vert0[0], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[1], vert0[1], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[2], vert0[2], 1.e-8);

		uvw[0] = 1;
		TLSM.computeTransformedCoords(uvw, xyz);
		BOOST_CHECK_CLOSE(xyz[0], vert1[0], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[1], vert1[1], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[2], vert1[2], 1.e-8);

		uvw[0] = 0;
		uvw[1] = 1;
		TLSM.computeTransformedCoords(uvw, xyz);
		BOOST_CHECK_CLOSE(xyz[0], vert2[0], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[1], vert2[1], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[2], vert2[2], 1.e-8);

		uvw[1] = 0;
		uvw[2] = 1;
		TLSM.computeTransformedCoords(uvw, xyz);
		BOOST_CHECK_CLOSE(xyz[0], vert3[0], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[1], vert3[1], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[2], vert3[2], 1.e-8);

		// Really a regression test
		uvw[0] = 0.25;
		uvw[1] = 1. / 6;
		uvw[2] = 5. / 12;
		TLSM.computeTransformedCoords(uvw, xyz);
		BOOST_CHECK_CLOSE(xyz[0], 1.284103709571, 1.e-8);
		BOOST_CHECK_CLOSE(xyz[1], 1.108135846478, 1.e-8);
		BOOST_CHECK_CLOSE(xyz[2], 1.902326582363, 1.e-8);

		emInt verts2[] = { 3, 1, 2, 0 };
		TLSM.setupCoordMapping(verts2);

		uvw[0] = uvw[1] = uvw[2] = 0;

		BOOST_CHECK_CLOSE(xyz[0], vert3[0], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[1], vert3[1], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[2], vert3[2], 1.e-8);

		uvw[0] = 1;
		TLSM.computeTransformedCoords(uvw, xyz);
		BOOST_CHECK_CLOSE(xyz[0], vert1[0], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[1], vert1[1], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[2], vert1[2], 1.e-8);

		uvw[0] = 0;
		uvw[1] = 1;
		TLSM.computeTransformedCoords(uvw, xyz);
		BOOST_CHECK_CLOSE(xyz[0], vert2[0], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[1], vert2[1], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[2], vert2[2], 1.e-8);

		uvw[1] = 0;
		uvw[2] = 1;
		TLSM.computeTransformedCoords(uvw, xyz);
		BOOST_CHECK_CLOSE(xyz[0], vert0[0], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[1], vert0[1], 1.e-8);
		BOOST_CHECK_CLOSE(xyz[2], vert0[2], 1.e-8);

		// Really a regression test
		uvw[0] = 0.25;
		uvw[1] = 1. / 6;
		uvw[2] = 1. / 6;
		TLSM.computeTransformedCoords(uvw, xyz);
		BOOST_CHECK_CLOSE(xyz[0], 1.45180347311979, 1.e-8);
		BOOST_CHECK_CLOSE(xyz[1], 1.35250275357685, 1.e-8);
		BOOST_CHECK_CLOSE(xyz[2], 2.11982908459981, 1.e-8);
	}

	BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(LagrangeCubicMapping)
	static void cubicXYZ(const double uvw[3], double xyz[3]) {
		const double& u = uvw[0];
		const double& v = uvw[1];
		const double& w = uvw[2];
		xyz[0] = u * u * u + 2 * u * u * v + 4 * u * v * v - 30 * u * v * w
				+ 5 * u * u + 6 * u + 3;
		xyz[1] = v * v * v + 2 * v * v * w + 4 * v * w * w - 54 * u * v * w
				+ 5 * v * v + 6 * v + 7;
		xyz[2] = w * w * w + 2 * w * w * u + 4 * w * u * u - 72 * u * v * w
				+ 5 * w * w + 6 * w + 10;
	}

	BOOST_AUTO_TEST_CASE(CubicTet) {
		CubicMesh CM(0, 0, 0, 0, 0, 0, 0, 0);
		LagrangeCubicTetMapping LCTM(&CM);
		// Test for known cubic functions:
		// x = u^3 + 2 u^2 v + 4 u v^2 - 30 u v w + 5 u^2 + 6 u + 3
		// y = v^3 + 2 v^2 w + 4 v w^2 - 54 u v w + 5 v^2 + 6 v + 7
		// z = w^3 + 2 w^2 u + 4 w u^2 - 72 u v w + 5 w^2 + 6 w + 10

		// So the nodal values are just these functions evaluated at nodes, where for the
		// canonical tet, we have:
		double uvw[][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 0 }, {
				2 / 3., 1 / 3., 0 },
												{ 1 / 3., 2 / 3., 0 }, { 0, 2 / 3., 1 / 3. },
												{ 0, 1 / 3., 2 / 3. }, { 1 / 3., 0, 2 / 3. },
												{ 2 / 3., 0, 1 / 3. }, { 2 / 3., 0, 0 },
												{ 1 / 3., 0, 0 }, { 0, 2 / 3., 0 },
												{ 0, 1 / 3., 0 }, {
														0, 0, 2 / 3. },
												{ 0, 0, 1 / 3. }, { 1 / 3., 1 / 3., 1 / 3. }, { 1 / 3.,
																																				1 / 3.,
																																				0 },
												{ 0, 1 / 3., 1 / 3. }, { 1 / 3., 0, 1 / 3. } };
		double xyz[20][3] = { 0 };

		// Now the nodal values, by simple functional evaluation.
		for (int ii = 0; ii < 20; ii++) {
			cubicXYZ(uvw[ii], xyz[ii]);
		}
		LCTM.setNodalValues(xyz);

		for (int iFunc = 0; iFunc < 20; iFunc++) {
			for (int iNode = 0; iNode < 20; iNode++) {
				double value = LCTM.computeBasisFunction(iFunc, uvw[iNode]);
				BOOST_CHECK_CLOSE(1 + value, 1 + (iFunc == iNode ? 1 : 0), 1.e-8);
			}
		}

		double testUVW[] = { 1. / M_PI, 1. / M_E, 0.5 * (1 - 1. / M_PI - 1. / M_E) };
		double LCTMxyz[3], funcxyz[3];
		LCTM.computeTransformedCoords(testUVW, LCTMxyz);
		cubicXYZ(testUVW, funcxyz);
		BOOST_CHECK_CLOSE(LCTMxyz[0], funcxyz[0], 1.e-8);
		BOOST_CHECK_CLOSE(LCTMxyz[1], funcxyz[1], 1.e-8);
		BOOST_CHECK_CLOSE(LCTMxyz[2], funcxyz[2], 1.e-8);
	}
	BOOST_AUTO_TEST_SUITE_END()
