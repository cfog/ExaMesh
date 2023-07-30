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

//#include "mpi.h"
#include <boost/mpi.hpp>
//#include <boost/timer/timer.hpp>
#include <chrono>
#include <fstream>


static void triUnitNormal(const double coords0[], const double coords1[],
		const double coords2[], double normal[]) {
	double edge01[] = DIFF(coords1, coords0);
	double edge02[] = DIFF(coords2, coords0);
	CROSS(edge01, edge02, normal);
	NORMALIZE(normal);
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
	for (emInt prism = 0; prism < numPrisms(); prism++) {
		const emInt* const prismVerts = getPrismConn(prism);
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
//		assert(volume > 0);
		for (int ii = 0; ii < 6; ii++) {
			vertVolume[prismVerts[ii]] += volume;
			assert(solids[ii] > 0);
			vertSolidAngle[prismVerts[ii]] += solids[ii];
		}
	} // Done with prisms

	// Iterate over hexahedra
	for (emInt hex = 0; hex < numHexes(); hex++) {
		const emInt* const hexVerts = getHexConn(hex);
		double norm1045[3], norm2156[3], norm3267[3], norm0374[3], norm0123[3],
				norm7654[3];
		double coords0[3], coords1[3], coords2[3], coords3[3], coords4[3],
				coords5[3], coords6[3], coords7[3];
		getCoords(hexVerts[0], coords0);
		getCoords(hexVerts[1], coords1);
		getCoords(hexVerts[2], coords2);
		getCoords(hexVerts[3], coords3);
		getCoords(hexVerts[4], coords4);
		getCoords(hexVerts[5], coords5);
		getCoords(hexVerts[6], coords6);
		getCoords(hexVerts[7], coords7);
		quadUnitNormal(coords1, coords0, coords4, coords5, norm1045);
		quadUnitNormal(coords2, coords1, coords5, coords6, norm2156);
		quadUnitNormal(coords3, coords2, coords6, coords7, norm3267);
		quadUnitNormal(coords0, coords3, coords7, coords4, norm0374);
		quadUnitNormal(coords0, coords1, coords2, coords3, norm0123);
		quadUnitNormal(coords7, coords6, coords5, coords4, norm7654);

		double diheds[12];
		// Dihedrals are in the order: 01, 12, 23, 30, 04, 15, 26, 37, 45, 56, 67, 74
		diheds[0] = safe_acos(-DOT(norm1045, norm0123));
		diheds[1] = safe_acos(-DOT(norm2156, norm0123));
		diheds[2] = safe_acos(-DOT(norm3267, norm0123));
		diheds[3] = safe_acos(-DOT(norm0374, norm0123));
		diheds[4] = safe_acos(-DOT(norm1045, norm0374));
		diheds[5] = safe_acos(-DOT(norm2156, norm1045));
		diheds[6] = safe_acos(-DOT(norm3267, norm2156));
		diheds[7] = safe_acos(-DOT(norm0374, norm3267));
		diheds[8] = safe_acos(-DOT(norm1045, norm7654));
		diheds[9] = safe_acos(-DOT(norm2156, norm7654));
		diheds[10] = safe_acos(-DOT(norm3267, norm7654));
		diheds[11] = safe_acos(-DOT(norm0374, norm7654));


		// Solid angles are in the order: 0, 1, 2, 3, 4, 5, 6, 7
		double solids[8];
		solids[0] = diheds[3] + diheds[0] + diheds[4] - M_PI;
		solids[1] = diheds[0] + diheds[1] + diheds[5] - M_PI;
		solids[2] = diheds[1] + diheds[2] + diheds[6] - M_PI;
		solids[3] = diheds[2] + diheds[3] + diheds[7] - M_PI;
		solids[4] = diheds[11] + diheds[8] + diheds[4] - M_PI;
		solids[5] = diheds[8] + diheds[9] + diheds[5] - M_PI;
		solids[6] = diheds[9] + diheds[10] + diheds[6] - M_PI;
		solids[7] = diheds[10] + diheds[11] + diheds[7] - M_PI;

		double middle[] = { (coords0[0] + coords1[0] + coords2[0] + coords3[0]
													+ coords4[0] + coords5[0] + coords6[0] + coords7[0])
												/ 8,
												(coords0[1] + coords1[1] + coords2[1] + coords3[1]
													+ coords4[1] + coords5[1] + coords6[1] + coords7[1])
												/ 8,
												(coords0[2] + coords1[2] + coords2[2] + coords3[2]
													+ coords4[2] + coords5[2] + coords6[2] + coords7[2])
												/ 8 };
		double volume = pyrVolume(coords1, coords0, coords4, coords5, middle)
				+ pyrVolume(coords2, coords1, coords5, coords6, middle)
				+ pyrVolume(coords3, coords2, coords6, coords7, middle)
				+ pyrVolume(coords0, coords3, coords7, coords4, middle)
				+ pyrVolume(coords0, coords1, coords2, coords3, middle)
				+ pyrVolume(coords7, coords6, coords5, coords4, middle);
//		assert(volume > 0);
		for (int ii = 0; ii < 8; ii++) {
			vertVolume[hexVerts[ii]] += volume;
//			assert(solids[ii] > 0);
			vertSolidAngle[hexVerts[ii]] += solids[ii];
		}
	} // Done with hexahedra

	// Now loop over verts computing the length scale
	for (emInt vv = 0; vv < numVerts(); vv++) {
//		assert(vertVolume[vv] > 0 && vertSolidAngle[vv] > 0);
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

void ExaMesh::prettyPrintCellCount(size_t cells, const char* prefix) const {
	if (cells == 0) return;
	printf("%s = ", prefix);
	if (cells >> 30) {
		printf("%.2f B\n", cells / 1.e9);
	}
	else if (cells >> 20) {
		printf("%.2f M\n", cells / 1.e6);
	}
	else if (cells >> 10) {
		printf("%.2f K\n", cells / 1.e3);
	}
	else {
		printf("%lu \n", cells);
	}
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
	double start = exaTime();
	partitionCells(this, nParts, parts, vecCPD);
	double partitionTime = exaTime() - start;

	// Create new sub-meshes and refine them.
	double totalRefineTime = 0;
	double totalExtractTime = 0;
	size_t totalCells = 0;
	size_t totalTets = 0, totalPyrs = 0, totalPrisms = 0, totalHexes = 0;
	size_t totalFileSize = 0;
	struct RefineStats RS;
	double totalTime = partitionTime;
	emInt ii;
//#pragma omp parallel for schedule(dynamic) reduction(+: totalRefineTime, totalExtractTime, totalTets, totalPyrs, totalPrisms, totalHexes, totalCells) num_threads(8)
	for (ii = 0; ii < nParts; ii++) {
		start = exaTime();
//		char filename[100];
//		sprintf(filename, "/tmp/submesh%03d.vtk", ii);
//		writeVTKFile(filename);
		printf("Part %3d: cells %5d-%5d.\n", ii, parts[ii].getFirst(),
						parts[ii].getLast());
		std::unique_ptr<UMesh> pUM = createFineUMesh(numDivs, parts[ii], vecCPD,
																									RS);
		totalRefineTime += RS.refineTime;
		totalExtractTime += RS.extractTime;
		totalCells += RS.cells;
		totalTets += pUM->numTets();
		totalPyrs += pUM->numPyramids();
		totalPrisms += pUM->numPrisms();
		totalHexes += pUM->numHexes();
		totalFileSize += pUM->getFileImageSize();
		totalTime += exaTime() - start;
		printf("\nCPU time for refinement = %5.2F seconds\n",
						RS.refineTime);
		printf("                          %5.2F million cells / minute\n",
						(RS.cells / 1000000.) / (RS.refineTime / 60));

		// char filename[100];
		// sprintf(filename, "TestCases/Testsubmesh%03d.vtk", ii);
		// pUM->writeVTKFile(filename);
	}
	printf("\nDone parallel refinement with %d parts.\n", nParts);
	printf("Time for partitioning:           %10.3F seconds\n",
					partitionTime);
	printf("Time for coarse mesh extraction: %10.3F seconds\n",
					totalExtractTime);
	printf("Time for refinement:             %10.3F seconds\n",
					totalRefineTime);
	printf("Rate (refinement only):  %5.2F million cells / minute\n",
					(totalCells / 1000000.) / (totalRefineTime / 60));
	printf("Rate (overall):          %5.2F million cells / minute\n",
					(totalCells / 1000000.) / (totalTime / 60));

	if (totalFileSize >> 37) {
		printf("Total ugrid file size = %lu GB\n", totalFileSize >> 30);
	}
	else if (totalFileSize >> 30) {
		printf("Total ugrid file size = %.2f GB\n",
						(totalFileSize >> 20) / 1024.);
	}
	else {
		printf("Total ugrid file size = %lu MB\n", totalFileSize >> 20);
	}

	prettyPrintCellCount(totalCells, "Total cells");
	prettyPrintCellCount(totalTets, "Total tets");
	prettyPrintCellCount(totalPyrs, "Total pyrs");
	prettyPrintCellCount(totalPrisms, "Total prisms");
	prettyPrintCellCount(totalHexes, "Total hexes");
}

//void ExaMesh::buildFaceCellConnectivity() {
//	fprintf(stderr, "Starting to build face cell connectivity\n");
//	// Create a multimap that will hold all of the face data, in duplicate.
//	// The sort key will be a TriFaceVerts / QuadFaceVerts object; the payload
//	// will be a CellInfo object.
//	exa_multimap<TriFaceVerts, CellInfo> triFaceData;
//	exa_multimap<QuadFaceVerts, CellInfo> quadFaceData;
//
//	// Use the linear tag for element type, as this does no harm.
//	for (emInt ii = 0; ii < numTets(); ii++) {
//		const emInt *conn = getTetConn(ii);
//		triFaceData.insert(
//				std::make_pair(TriFaceVerts(conn[0], conn[1], conn[2]),
//												CellInfo(ii, TETRAHEDRON, 0)));
//		triFaceData.insert(
//				std::make_pair(TriFaceVerts(conn[0], conn[1], conn[3]),
//												CellInfo(ii, TETRAHEDRON, 1)));
//		triFaceData.insert(
//				std::make_pair(TriFaceVerts(conn[1], conn[2], conn[3]),
//												CellInfo(ii, TETRAHEDRON, 2)));
//		triFaceData.insert(
//				std::make_pair(TriFaceVerts(conn[2], conn[0], conn[3]),
//												CellInfo(ii, TETRAHEDRON, 3)));
//	}
//
//	for (emInt ii = 0; ii < numPyramids(); ii++) {
//		const emInt *conn = getPyrConn(ii);
//		quadFaceData.insert(
//				std::make_pair(QuadFaceVerts(conn[0], conn[1], conn[2], conn[3]),
//												CellInfo(ii, PYRAMID, 0)));
//		triFaceData.insert(
//				std::make_pair(TriFaceVerts(conn[0], conn[4], conn[1]),
//												CellInfo(ii, PYRAMID, 1)));
//		triFaceData.insert(
//				std::make_pair(TriFaceVerts(conn[1], conn[4], conn[2]),
//												CellInfo(ii, PYRAMID, 2)));
//		triFaceData.insert(
//				std::make_pair(TriFaceVerts(conn[2], conn[4], conn[3]),
//												CellInfo(ii, PYRAMID, 3)));
//		triFaceData.insert(
//				std::make_pair(TriFaceVerts(conn[3], conn[4], conn[0]),
//												CellInfo(ii, PYRAMID, 4)));
//	}
//
//	for (emInt ii = 0; ii < numPrisms(); ii++) {
//		const emInt *conn = getPrismConn(ii);
//		quadFaceData.insert(
//				std::make_pair(QuadFaceVerts(conn[2], conn[1], conn[4], conn[5]),
//												CellInfo(ii, PRISM, 0)));
//		quadFaceData.insert(
//				std::make_pair(QuadFaceVerts(conn[1], conn[0], conn[3], conn[4]),
//												CellInfo(ii, PRISM, 1)));
//		quadFaceData.insert(
//				std::make_pair(QuadFaceVerts(conn[0], conn[2], conn[5], conn[3]),
//												CellInfo(ii, PRISM, 2)));
//		triFaceData.insert(
//				std::make_pair(TriFaceVerts(conn[0], conn[1], conn[2]),
//												CellInfo(ii, PRISM, 3)));
//		triFaceData.insert(
//				std::make_pair(TriFaceVerts(conn[5], conn[4], conn[3]),
//												CellInfo(ii, PRISM, 4)));
//	}
//
//	for (emInt ii = 0; ii < numHexes(); ii++) {
//		const emInt *conn = getHexConn(ii);
//		quadFaceData.insert(
//				std::make_pair(QuadFaceVerts(conn[0], conn[1], conn[2], conn[3]),
//												CellInfo(ii, HEXAHEDRON, 0)));
//		quadFaceData.insert(
//				std::make_pair(QuadFaceVerts(conn[7], conn[6], conn[5], conn[4]),
//												CellInfo(ii, HEXAHEDRON, 1)));
//		quadFaceData.insert(
//				std::make_pair(QuadFaceVerts(conn[0], conn[4], conn[5], conn[1]),
//												CellInfo(ii, HEXAHEDRON, 2)));
//		quadFaceData.insert(
//				std::make_pair(QuadFaceVerts(conn[1], conn[5], conn[6], conn[2]),
//												CellInfo(ii, HEXAHEDRON, 3)));
//		quadFaceData.insert(
//				std::make_pair(QuadFaceVerts(conn[2], conn[6], conn[7], conn[3]),
//												CellInfo(ii, HEXAHEDRON, 4)));
//		quadFaceData.insert(
//				std::make_pair(QuadFaceVerts(conn[3], conn[7], conn[4], conn[0]),
//												CellInfo(ii, HEXAHEDRON, 5)));
//	}
//
//	for (emInt ii = 0; ii < numBdryTris(); ii++) {
//		const emInt *conn = getBdryTriConn(ii);
//		triFaceData.insert(
//				std::make_pair(TriFaceVerts(conn[0], conn[1], conn[2]),
//												CellInfo(ii, TRIANGLE, 0)));
//	}
//
//	for (emInt ii = 0; ii < numBdryQuads(); ii++) {
//		const emInt *conn = getBdryQuadConn(ii);
//		quadFaceData.insert(
//				std::make_pair(QuadFaceVerts(conn[0], conn[1], conn[2], conn[3]),
//												CellInfo(ii, QUADRILATERAL, 0)));
//	}
//	fprintf(stderr, "Done filling up face maps\n");
//	fprintf(stderr, "Done building face cell connectivity\n");
//}
void
ExaMesh:: TestMPI(const emInt &nDivs, const emInt &nParts, ParallelTester* tester, 
const char MeshType)
{
	vecPart parts;
	vecCellPartData vecCPD;

	std::set<int> triRotations;
	std::set<int> quadRotations;

	vecHashTri tris;
	vecHashQuad quads;

	std::vector<std::shared_ptr<ExaMesh>> submeshes;
	vecSharePtrUmesh refinedUMeshes;

	std::vector<TableTri2TableIndex2Index>  matchedTrisAllParts  (nParts); 
	std::vector<TableQuad2TableIndex2Index> matchedQuadsAllParts (nParts);



	//std::vector< TriFaceVerts , std::unordered_map<emInt,emInt>>  matchedTrisAllParts (nParts); 
	//std::vector< QuadFaceVerts, std::unordered_map<emInt,emInt>>  matchedQuadAllParts (nParts); 


	partitionCells(this, nParts, parts, vecCPD);


	
	this->partFaceMatching(parts, vecCPD, tris, quads);
	

	

	for (auto i = 0; i < nParts; i++)
	{
		auto coarse = this->extractCoarseMesh(parts[i], vecCPD, nDivs, tris[i], quads[i], i);
		submeshes.emplace_back(coarse.release()); 
	}

	assert(submeshes.size() == static_cast<std::size_t>(nParts));
	// Refine the mesh
	for (auto i = 0; i < nParts; i++)
	{
		

		if(MeshType=='C')
		{
			auto inputMesh = dynamic_cast<CubicMesh*> ((submeshes[i].get()));
					auto refineUmesh = std::make_shared<UMesh>(
			*(inputMesh), nDivs, i);
			refinedUMeshes.push_back(refineUmesh);
		}
		if(MeshType=='U')
		{
			auto inputMesh = dynamic_cast<UMesh*> ((submeshes[i].get()));
					auto refineUmesh = std::make_shared<UMesh>(
			*(inputMesh), nDivs, i);
			refinedUMeshes.push_back(refineUmesh);
		}


	}
	for (auto iPart = 0; iPart < nParts; iPart++)
	{
		TableTri2TableIndex2Index  matchedTris; 
		TableQuad2TableIndex2Index matchedQuads; 
		exa_set<TriFaceVerts> tri =
			refinedUMeshes[iPart]->getRefinedPartTris();
		exa_set<QuadFaceVerts> quads = refinedUMeshes[iPart]->getRefinedPartQuads();
		for (auto it = tri.begin(); it != tri.end(); it++)
		{
			std::unordered_map<emInt, emInt> localRemote;
			emInt whereToLook = it->getRemoteId();
			exa_set<TriFaceVerts> remoteTriSet =
				refinedUMeshes[whereToLook]->getRefinedPartTris();
			emInt rotation = getTriRotation(*it, remoteTriSet, nDivs);
			triRotations.insert(rotation);

			matchTri(*it, rotation, nDivs, remoteTriSet, localRemote);
			matchedTris.emplace(*it,localRemote); 
			//printMatchedTris(localRemote,iPart); 
			for (auto itmap = localRemote.begin();
				 itmap != localRemote.end(); itmap++)
			{
				assert(abs(refinedUMeshes[iPart]->getX(itmap->first) -
						   refinedUMeshes[whereToLook]->getX(itmap->second)) < TOLTEST);
				assert(abs(refinedUMeshes[iPart]->getY(itmap->first) -
						   refinedUMeshes[whereToLook]->getY(itmap->second)) < TOLTEST);
				assert(abs(refinedUMeshes[iPart]->getZ(itmap->first) -
						   refinedUMeshes[whereToLook]->getZ(itmap->second)) < TOLTEST);
			}
		}
		for (auto it = quads.begin(); it != quads.end(); it++)
		{
			std::unordered_map<emInt, emInt> localRemote;
			emInt whereToLook = it->getRemoteId();
			exa_set<QuadFaceVerts> remoteQuadSet =
				refinedUMeshes[whereToLook]->getRefinedPartQuads();
			emInt rotation = getQuadRotation(*it, remoteQuadSet, nDivs);
			quadRotations.insert(rotation);

			matchQuad(*it, rotation, nDivs, remoteQuadSet, localRemote);
			matchedQuads.emplace(*it,localRemote); 
		
			//printMatchedQuads(localRemote,iPart); 
			for (auto itmap = localRemote.begin();
				 itmap != localRemote.end(); itmap++)
			{
				assert(abs(refinedUMeshes[iPart]->getX(itmap->first) -
						   refinedUMeshes[whereToLook]->getX(itmap->second)) < TOLTEST);
				assert(abs(refinedUMeshes[iPart]->getY(itmap->first) -
						   refinedUMeshes[whereToLook]->getY(itmap->second)) < TOLTEST);
				assert(abs(refinedUMeshes[iPart]->getZ(itmap->first) -
						   refinedUMeshes[whereToLook]->getZ(itmap->second)) < TOLTEST);
			}
		}
		matchedTrisAllParts[iPart] = matchedTris; 
		matchedQuadsAllParts[iPart]= matchedQuads; 
	}
	// Which rotation cases are covered?
	std::cout << "Covered Tri Rotations: " << std::endl;
	for (auto it = triRotations.begin();
		 it != triRotations.end(); it++)
	{
		std::cout << *it << " ";
	};
	std::cout << std::endl;
	std::cout << "Covered Quad Rotations: " << std::endl;
	for (auto it = quadRotations.begin();
		 it != quadRotations.end(); it++)
	{
		std::cout << *it << " ";
	};
	std::cout << std::endl;

	tester->setMatchedTris (matchedTrisAllParts); 
	tester->setMatchedQuads(matchedQuadsAllParts); 
}

void 
ExaMesh::refineForMPI( const int numDivs ,ParallelTester* tester,const char MeshType) 
const
{
	boost::mpi::environment   env; 
	boost::mpi::communicator  world;
	
	auto start = exaTime();
	
	std::vector<boost::mpi::request>      triReqs;
	std::vector<boost::mpi::request>      quadReqs;
	std::vector<std::unique_ptr<UMesh>>   refinedMeshVec;  

	vecPart         parts; 
	vecCellPartData vecCPD; 
	
	std::size_t    vecCPDSize; 
	std::size_t    trisSize; 
	std::size_t    quadsSize; 
	std::size_t    nParts = world.size();

	hashTri  trisS; 
	hashQuad quadsS;

	vecTri   triV;
	vecQuad  quadV;  

	int MASTER = 0; 
	int tag    = 0; 

	intToVecTri  remoteTovecTris;
	intToVecQuad remoteTovecQuads; 

	std::set<int> triNeighbrs;  
	std::set<int> quadNeighbrs; 

	
	std::vector<vecTri>  trisTobeRcvd; 
	std::vector<vecQuad> quadsTobeRcvd;

	hashTri  recvdTris;
	hashQuad recvdQuads; 

	TableTri2TableIndex2Index   matchedTris; 
	TableQuad2TableIndex2Index  matchedQuads; 

	if(world.rank()==MASTER)
	{
		partitionCells(this, nParts, parts,vecCPD); 

		vecCPDSize = vecCPD.size(); 

		assert(vecCPDSize>0); 

		vecHashTri  VectrisHash; 
		vecHashQuad VecquadsHash;

		vecVecTri   VecTriVec; 
		vecVecQuad  vecQuadVec; 

		this->partFaceMatching(parts,vecCPD,VectrisHash,VecquadsHash); 
		for(auto  itri=0 ; itri<VectrisHash.size(); itri++)
		{
			vecTri TriVec; 
			SetToVector(VectrisHash[itri],TriVec); 
			VecTriVec.emplace_back(TriVec); 
		}
		for(auto iquad=0 ; iquad<VecquadsHash.size(); iquad++)
		{
			vecQuad QuadVec; 
			SetToVector(VecquadsHash[iquad],QuadVec); 
			vecQuadVec.emplace_back(QuadVec); 
		}

		trisS = VectrisHash[0]; // For MASTER 
		quadsS= VecquadsHash[0];
		
		for(auto irank=1 ; irank<world.size();irank++)
		{
			world.send(irank,tag,parts[irank]); 
			world.send(irank,tag,VecTriVec[irank]); 
			world.send(irank,tag,vecQuadVec[irank]); 
		}
	}
	else
	{
		parts.resize(world.size()); 
		world.recv(MASTER,tag,parts[world.rank()]);
		world.recv(MASTER,tag,triV); 
		world.recv(MASTER,tag,quadV);
		vectorToSet(triV,trisS);
		vectorToSet(quadV,quadsS); 
	}

	boost::mpi::broadcast(world,vecCPDSize,MASTER);

	if(world.rank()==MASTER)
	{
		for(auto irank=1 ; irank<world.size();irank++)
		{
			world.send(irank,tag,vecCPD); 
		}
	}
	if(world.rank()!= MASTER)
	{
		vecCPD.resize(vecCPDSize); 
		world.recv(MASTER,tag,vecCPD); 
	}

	auto coarse= this->extractCoarseMesh(parts[world.rank()],vecCPD,numDivs, 
	 trisS,quadsS,world.rank()); 
	
	if(MeshType=='C')
	{
		auto refinedMesh = std::make_unique<UMesh>(
	 		*(dynamic_cast<CubicMesh*> (coarse.release())), 
			numDivs, world.rank());
			refinedMeshVec.emplace_back(refinedMesh.release()); 
	}
	if(MeshType=='U')
	{
		auto refinedMesh = std::make_unique<UMesh>(
	 		*(dynamic_cast<UMesh*> (coarse.release())), 
			numDivs, world.rank());
			refinedMeshVec.emplace_back(refinedMesh.release()); 
	}


	auto tris  = refinedMeshVec[0]->getRefinedPartTris();
	auto quads = refinedMeshVec[0]->getRefinedPartQuads();

	
	buildTrisMap(tris,remoteTovecTris,triNeighbrs);
	buildQuadsMap(quads,remoteTovecQuads,quadNeighbrs); 
	

	trisTobeRcvd.resize (triNeighbrs.size()); 
	quadsTobeRcvd.resize(quadNeighbrs.size());

	int Trijj=0 ; 
	for(auto isource:triNeighbrs)
	{
		auto findSource = remoteTovecTris.find(isource); 
		// same amount of message will be sent and received from other prcoeszsor
		// Be cuatious if this assumption later on will be break
		trisTobeRcvd[Trijj].resize(findSource->second.size()); 
		Trijj++; 
	}

	int Quadjj=0; 
	for(auto isource:quadNeighbrs)
	{
		auto findSource = remoteTovecQuads.find(isource); 
		// same amount of message will be sent and received from other prcoeszsor
		// Be cuatious if this assumption later on will be break
		quadsTobeRcvd[Quadjj].resize(findSource->second.size()); 
		Quadjj++; 
	}

	for(auto tri:remoteTovecTris)
	{
		boost::mpi::request req= world.isend(tri.first,tag,tri.second);
		triReqs.push_back(req);
	}

	for(auto quad:remoteTovecQuads)
	{
		boost::mpi::request req = world.isend(quad.first,tag,quad.second);
		quadReqs.push_back(req);  
	}

	int Quadkk=0 ; 
	for(auto source: quadNeighbrs)
	{
		boost::mpi::request req = world.irecv(source,tag,quadsTobeRcvd[Quadkk]);
		quadReqs.push_back(req); 
		Quadkk++; 
	}


	int Trikk=0 ; 
	for(auto source: triNeighbrs)
	{
		
		boost::mpi::request req= world.irecv(source,tag,trisTobeRcvd[Trikk]);
		triReqs.push_back(req);
		Trikk++; 
	}

	boost::mpi::wait_all(triReqs.begin(),triReqs.end()); 

	for(const auto& tri: trisTobeRcvd)
	{

		recvdTris.insert(tri.begin(),tri.end()); 
	}

	for(auto it=tris.begin(); it!=tris.end(); it++)
	{
		std::unordered_map<emInt, emInt> localRemote;
		//int rotation = getTriRotation(*it,recvdTris,numDivs);
		//matchTri(*it,rotation,numDivs,recvdTris,localRemote); 
		findRotationAndMatchTris(*it,numDivs,recvdTris,localRemote); 
		
#ifndef NDEBUG
		//matchedTris.emplace(*it,localRemote);
#endif //NDEBUG
	}

	boost::mpi::wait_all(quadReqs.begin(),quadReqs.end()); 

	for(const auto& quad:quadsTobeRcvd)
	{
		recvdQuads.insert(quad.begin(),quad.end()); 
	}

	for(auto it=quads.begin(); it!=quads.end();it++)
	{
		std::unordered_map<emInt,emInt> localRemote; 
		//int rotation = getQuadRotation(*it,recvdQuads,numDivs); 
		//matchQuad(*it,rotation,numDivs,recvdQuads,localRemote);
		findRotationAndMatchQuads(*it,numDivs,recvdQuads,localRemote); 
#ifndef NDEBUG			
		//matchedQuads.emplace(*it,localRemote); 
#endif
	}
	
	double time = exaTime() - start;
	double maximumTime; 
	boost::mpi::reduce(world, time, maximumTime,boost:: mpi::maximum<double>(), MASTER);
	
	if (world.rank() == MASTER) 
	{	
    	fprintf(stderr, "MAXIMUM CPU time for refinement = %5.2F seconds, for numDivs of %d\n", maximumTime,numDivs);
 	} 
	world.barrier();
	fprintf(stderr, "Rank of: %d has %lu tris and has %lu quads.\n", world.rank(), tris.size(), quads.size());
#ifndef NDEBUG				
	tester->testMatchedTris(matchedTris,world.rank()); 
	tester->testMatchedQuads(matchedQuads,world.rank()); 				
#endif

}