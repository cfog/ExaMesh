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

static void checkExpectedSize(const UMesh& UM) {
	BOOST_CHECK_EQUAL(UM.maxNVerts(), UM.numVerts());
	BOOST_CHECK_EQUAL(UM.maxNBdryTris(), UM.numBdryTris());
	BOOST_CHECK_EQUAL(UM.maxNBdryQuads(), UM.numBdryQuads());
	BOOST_CHECK_EQUAL(UM.maxNTets(), UM.numTets());
	BOOST_CHECK_EQUAL(UM.maxNPyrs(), UM.numPyrs());
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
