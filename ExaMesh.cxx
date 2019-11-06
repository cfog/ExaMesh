/*
 * ExaMesh.cxx
 *
 *  Created on: Oct. 3, 2019
 *      Author: cfog
 */

#include <assert.h>
#include <memory>
#include <vector>
#include <iostream>

using std::cout;
using std::endl;

#include "ExaMesh.h"
#include "GeomUtils.h"
#include "Part.h"
#include "UMesh.h"



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
	for (emInt pyr = 0; pyr < numPyramids(); pyr++) {
		const emInt* const pyrVerts = getPyrConn(pyr);
		double norm0123[3], norm014[3], norm124[3], norm234[3], norm304[3];
		double coords0[3], coords1[3], coords2[3], coords3[3], coords4[3];
		getCoords(pyrVerts[0], coords0);
		getCoords(pyrVerts[1], coords1);
		getCoords(pyrVerts[2], coords2);
		getCoords(pyrVerts[3], coords3);
		getCoords(pyrVerts[4], coords4);
		quadUnitNormal(coords0, coords1, coords2, coords3, norm0123);
		triUnitNormal(coords0, coords1, coords4, norm014);
		triUnitNormal(coords1, coords2, coords4, norm124);
		triUnitNormal(coords2, coords3, coords4, norm234);
		triUnitNormal(coords3, coords0, coords4, norm304);

		double diheds[8];
		// Dihedrals are in the order: 01, 04, 12, 14, 23, 24, 30, 34
		diheds[0] = safe_acos(-DOT(norm0123, norm014));
		diheds[1] = safe_acos(-DOT(norm014, norm304));
		diheds[2] = safe_acos(-DOT(norm0123, norm124));
		diheds[3] = safe_acos(-DOT(norm124, norm014));
		diheds[4] = safe_acos(-DOT(norm0123, norm234));
		diheds[5] = safe_acos(-DOT(norm234, norm124));
		diheds[6] = safe_acos(-DOT(norm0123, norm304));
		diheds[7] = safe_acos(-DOT(norm304, norm234));

		// Solid angles are in the order: 0, 1, 2, 3, 4
		double solids[5];
		solids[0] = diheds[0] + diheds[1] + diheds[6] - M_PI;
		solids[1] = diheds[0] + diheds[2] + diheds[3] - M_PI;
		solids[2] = diheds[2] + diheds[4] + diheds[5] - M_PI;
		solids[3] = diheds[4] + diheds[6] + diheds[7] - M_PI;
		solids[4] = diheds[1] + diheds[3] + diheds[5] + diheds[7] - 2 * M_PI;

		double volume = pyrVolume(coords0, coords1, coords2, coords3, coords4);
		assert(volume > 0);
		for (int ii = 0; ii < 5; ii++) {
			vertVolume[pyrVerts[ii]] += volume;
			assert(solids[ii] > 0);
			vertSolidAngle[pyrVerts[ii]] += solids[ii];
		}
	} // Done with pyramids

	// Iterate over prisms
	for (emInt prism = 0; prism < numPyramids(); prism++) {
		const emInt* const prismVerts = getPyrConn(prism);
		double norm1034[3], norm2145[3], norm0253[3], norm012[3], norm543[3];
		double coords0[3], coords1[3], coords2[3], coords3[3], coords4[3],
				coords5[3];
		getCoords(prismVerts[0], coords0);
		getCoords(prismVerts[1], coords1);
		getCoords(prismVerts[2], coords2);
		getCoords(prismVerts[3], coords3);
		getCoords(prismVerts[4], coords4);
		getCoords(prismVerts[5], coords5);
		quadUnitNormal(coords1, coords0, coords3, coords4, norm1034);
		quadUnitNormal(coords2, coords1, coords4, coords5, norm2145);
		quadUnitNormal(coords0, coords2, coords5, coords3, norm0253);
		triUnitNormal(coords0, coords1, coords2, norm012);
		triUnitNormal(coords5, coords4, coords3, norm543);

		double diheds[9];
		// Dihedrals are in the order: 01, 12, 20, 03, 14, 25, 34, 45, 53
		diheds[0] = safe_acos(-DOT(norm1034, norm012));
		diheds[1] = safe_acos(-DOT(norm2145, norm012));
		diheds[2] = safe_acos(-DOT(norm0253, norm012));
		diheds[3] = safe_acos(-DOT(norm0253, norm1034));
		diheds[4] = safe_acos(-DOT(norm1034, norm2145));
		diheds[5] = safe_acos(-DOT(norm2145, norm0253));
		diheds[6] = safe_acos(-DOT(norm1034, norm543));
		diheds[7] = safe_acos(-DOT(norm2145, norm543));
		diheds[8] = safe_acos(-DOT(norm0253, norm543));

		// Solid angles are in the order: 0, 1, 2, 3, 4, 5
		double solids[6];
		solids[0] = diheds[0] + diheds[2] + diheds[3] - M_PI;
		solids[1] = diheds[0] + diheds[1] + diheds[4] - M_PI;
		solids[2] = diheds[1] + diheds[2] + diheds[5] - M_PI;
		solids[3] = diheds[6] + diheds[8] + diheds[3] - M_PI;
		solids[4] = diheds[6] + diheds[7] + diheds[4] - M_PI;
		solids[5] = diheds[7] + diheds[8] + diheds[5] - M_PI;

		double middle[] = { (coords0[0] + coords1[0] + coords2[0] + coords3[0]
													+ coords4[0] + coords5[0])
												/ 6,
												(coords0[1] + coords1[1] + coords2[1] + coords3[1]
													+ coords4[1] + coords5[1])
												/ 6,
												(coords0[2] + coords1[2] + coords2[2] + coords3[2]
													+ coords4[2] + coords5[2])
												/ 6 };
		double volume = tetVolume(coords0, coords1, coords2, middle)
				+ tetVolume(coords5, coords4, coords3, middle)
				+ pyrVolume(coords1, coords0, coords3, coords4, middle)
				+ pyrVolume(coords2, coords1, coords4, coords5, middle)
				+ pyrVolume(coords0, coords2, coords5, coords3, middle);
		assert(volume > 0);
		for (int ii = 0; ii < 6; ii++) {
			vertVolume[prismVerts[ii]] += volume;
			assert(solids[ii] > 0);
			vertSolidAngle[prismVerts[ii]] += solids[ii];
		}
	} // Done with prisms

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

MeshSize ExaMesh::computeFineMeshSize(const int nDivs) const {
	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = numBdryVerts();
	MSIn.nVerts = numVerts();
	MSIn.nBdryTris = numBdryTris();
	MSIn.nBdryQuads = numBdryQuads();
	MSIn.nTets = numTets();
	MSIn.nPyrs = numPyramids();
	MSIn.nPrisms = numPrisms();
	MSIn.nHexes = numHexes();
	bool sizesOK = ::computeMeshSize(MSIn, nDivs, MSOut);
	if (!sizesOK) exit(2);

	return MSOut;
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

void ExaMesh::refineForParallel(const emInt numDivs,
		const emInt maxCellsPerPart) const {
	// Find size of output mesh
	size_t numCells = numTets() + numPyramids() + numHexes() + numPrisms();
	size_t outputCells = numCells * (numDivs * numDivs * numDivs);

	// Calc number of parts.  This funky formula makes it so that, if you need
	// N*maxCells, you'll get N parts.  With N*maxCells + 1, you'll get N+1.
	emInt nParts = (outputCells - 1) / maxCellsPerPart + 1;
	if (nParts > numCells) nParts = numCells;

	// Partition the mesh.
	std::vector<Part> parts;
	std::vector<CellPartData> vecCPD;
	partitionCells(this, nParts, parts, vecCPD);

	// Create new sub-meshes and refine them.
	double totalTime = 0;
	size_t totalCells = 0;
	emInt ii;
//#pragma omp parallel for schedule(dynamic) reduction(+: totalTime, totalCells)
	for (ii = 0; ii < nParts; ii++) {
//		char filename[100];
//		sprintf(filename, "/tmp/submesh%03d.vtk", ii);
//		writeVTKFile(filename);
		printf(
				"Part %3d: cells %5d-%5d.  (%6.3f,%6.3f,%6.3f) (%6.3f,%6.3f,%6.3f)\n",
				ii, parts[ii].getFirst(), parts[ii].getLast(), parts[ii].getXmin(),
				parts[ii].getYmin(), parts[ii].getZmin(), parts[ii].getXmax(),
				parts[ii].getYmax(), parts[ii].getZmax());
		double time;
		size_t cells;
		std::unique_ptr<UMesh> pUM = createFineUMesh(numDivs, parts[ii], vecCPD,
																									time, cells);
		totalTime += time;
		totalCells += cells;
//		char filename[100];
//		sprintf(filename, "/tmp/fine-submesh%03d.vtk", ii);
//		pUM->writeVTKFile(filename);
	}
	fprintf(stderr, "\nDone parallel refinement with %d parts.\n", nParts);
	fprintf(stderr, "CPU time for refinement = %5.2F seconds\n", totalTime);
	fprintf(stderr, "                          %5.2F million cells / minute\n",
					(totalCells / 1000000.) / (totalTime / 60));

}

