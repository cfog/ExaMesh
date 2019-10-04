/*
 * ExaMesh.cxx
 *
 *  Created on: Oct. 3, 2019
 *      Author: cfog
 */

#include <assert.h>
#include <vector>
#include <iostream>

using std::cout;
using std::endl;

#include "examesh.h"
#include "GeomUtils.h"

static void triUnitNormal(const double coords0[], const double coords1[],
		const double coords2[], double normal[]) {
	double edge01[] = DIFF(coords1, coords0);
	double edge02[] = DIFF(coords2, coords0);
	CROSS(edge01, edge02, normal);
	NORMALIZE(normal);
}

static double tetVolume(const double coords0[], const double coords1[],
		const double coords2[], const double coords3[]) {
	double edge01[] = DIFF(coords1, coords0);
	double edge02[] = DIFF(coords2, coords0);
	double edge03[] = DIFF(coords3, coords0);
	double normal[3];
	CROSS(edge01, edge02, normal);
	return DOT(normal,edge03) / 6;
}

static void quadUnitNormal(const double coords0[], const double coords1[],
		const double coords2[], const double coords3[], double normal[]) {
	double vecB[3], vecC[3];
	for (int ii = 0; ii < 3; ii++) {
		vecB[ii] = 0.25 * (coords0[ii] + coords3[ii] - coords1[ii] - coords2[ii]);
		vecC[ii] = 0.25 * (coords0[ii] + coords1[ii] - coords3[ii] - coords2[ii]);
	}
	CROSS(vecB, vecC, normal);
	NORMALIZE(normal);
}

static double pyrVolume(const double coords0[], const double coords1[],
		const double coords2[], const double coords3[], double coords4[]) {
	// Point 4 is the apex.
	double vecB[3], vecC[3], vecE[3];
	for (int ii = 0; ii < 3; ii++) {
		vecB[ii] = 0.25 * (coords0[ii] + coords3[ii] - coords1[ii] - coords2[ii]);
		vecC[ii] = 0.25 * (coords0[ii] + coords1[ii] - coords3[ii] - coords2[ii]);
		vecE[ii] = coords4[ii]
				- 0.25 * (coords0[ii] + coords1[ii] + coords2[ii] + coords3[ii]);
	}
	double normal[3];
	CROSS(vecB, vecC, normal);
	return DOT(normal, vecE) / 0.75;
}


// TODO  Transplant into a new ExaMesh.cxx
void ExaMesh::setupLengthScales() {
	if (!m_lenScale) {
		m_lenScale = new double[numVerts()];
	}
	std::vector<double> vertVolume(numVerts(), 0);
	std::vector<double> vertSolidAngle(numVerts(), 0);

	// Iterate over tets
	for (emInt tet = 0; tet < numTets(); tet++) {
		const emInt* const tetVerts = getTetConn(tet);
		double normABC[3], normADB[3], normBDC[3], normCDA[3];
		double coordsA[3], coordsB[3], coordsC[3], coordsD[3];
		getCoords(tetVerts[0], coordsA);
		getCoords(tetVerts[1], coordsB);
		getCoords(tetVerts[2], coordsC);
		getCoords(tetVerts[3], coordsD);
		triUnitNormal(coordsA, coordsB, coordsC, normABC);
		triUnitNormal(coordsA, coordsD, coordsB, normADB);
		triUnitNormal(coordsB, coordsD, coordsC, normBDC);
		triUnitNormal(coordsC, coordsD, coordsA, normCDA);

		// Dihedrals are in the order: 01, 02, 03, 12, 13, 23
		double diheds[6];
		diheds[0] = safe_acos(-DOT(normABC, normADB));
		diheds[1] = safe_acos(-DOT(normABC, normCDA));
		diheds[2] = safe_acos(-DOT(normADB, normCDA));
		diheds[3] = safe_acos(-DOT(normABC, normBDC));
		diheds[4] = safe_acos(-DOT(normADB, normBDC));
		diheds[5] = safe_acos(-DOT(normBDC, normCDA));

		// Solid angles are in the order: 0, 1, 2, 3
		double solids[4];
		solids[0] = diheds[0] + diheds[1] + diheds[2] - M_PI;
		solids[1] = diheds[0] + diheds[3] + diheds[4] - M_PI;
		solids[2] = diheds[1] + diheds[3] + diheds[5] - M_PI;
		solids[3] = diheds[2] + diheds[4] + diheds[5] - M_PI;

		double volume = tetVolume(coordsA, coordsB, coordsC, coordsD);
		assert(volume > 0);
		for (int ii = 0; ii < 4; ii++) {
			vertVolume[tetVerts[ii]] += volume;
			assert(solids[ii] > 0);
			vertSolidAngle[tetVerts[ii]] += solids[ii];
		}
	} // Done looping over tetrahedra

	// Iterate over pyramids
	assert(numPyramids() == 0);

	// Iterate over prisms
	assert(numPrisms() == 0);

	// Iterate over hexahedra
	assert(numHexes() == 0);

	// Now loop over verts computing the length scale
	for (emInt vv = 0; vv < numVerts(); vv++) {
		assert(vertVolume[vv] > 0 && vertSolidAngle[vv] > 0);
		double volume = vertVolume[vv] * (4 * M_PI) / vertSolidAngle[vv];
		double radius = cbrt(volume / (4 * M_PI / 3.));
		m_lenScale[vv] = radius;
	}
}

void ExaMesh::printMeshSizeStats() {
	cout << "Mesh has:" << endl;
	cout.width(16);
	cout << numVerts() << " verts" << endl;
	cout.width(16);
	cout << numBdryTris() << " bdry tris" << endl;
	cout.width(16);
	cout << numBdryQuads() << " bdry quads" << endl;
	cout.width(16);
	cout << numTets() << " tets" << endl;
	cout.width(16);
	cout << numPyramids() << " pyramids" << endl;
	cout.width(16);
	cout << numPrisms() << " prisms" << endl;
	cout.width(16);
	cout << numHexes() << " hexes" << endl;
	cout.width(16);
	cout << numTets() + numPyramids() + numPrisms() + numHexes()
				<< " total cells " << endl;
}



