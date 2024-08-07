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

#include "exa-defs.h"
#include "ExaMesh.h"
#include "GeomUtils.h"
#include "Part.h"
#include "UMesh.h"

//#include "mpi.h"
#include <boost/mpi.hpp>
//#include <boost/timer/timer.hpp>
#include <chrono>
#include <fstream>
#include "resultGenerator.cxx"

std::unique_ptr<ExaMesh> ExaMesh::readMeshFromFile(
		const std::string fileName,
		const std::string fileSuffix,
		const std::string fileInfix) {
	std::unique_ptr<ExaMesh> retVal;
	if (fileSuffix == "cgns") {
		retVal = std::make_unique<CubicMesh>(fileName);
	}
	else {
		// This calls the UnstructuredMeshAnalyzer FileWrapper,
		// which can read a lot more than ugrid files.
		retVal = std::make_unique<UMesh>(fileName, fileSuffix, fileInfix);
	}
	return retVal;
}
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
		vecB[ii] = 0.25
				* (coords0[ii] + coords3[ii] - coords1[ii] - coords2[ii]);
		vecC[ii] = 0.25
				* (coords0[ii] + coords1[ii] - coords3[ii] - coords2[ii]);
	}
	CROSS(vecB, vecC, normal);
	NORMALIZE(normal);
}

void ExaMesh::setupLengthScales() {
	if (!m_lenScale) {
		m_lenScale = new double[numVerts()];
	}
	for (emInt ii = 0; ii < numVerts(); ii++) {
		m_lenScale[ii] = DBL_MAX;
	}
	std::vector<double> vertVolume(numVerts(), DBL_MAX / 100);
	std::vector<double> vertSolidAngle(numVerts(), 1);

	// Iterate over tets
	for (emInt tet = 0; tet < numTets(); tet++) {
		const emInt *const tetVerts = getTetConn(tet);
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
		const emInt *const pyrVerts = getPyrConn(pyr);
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
		const emInt *const prismVerts = getPrismConn(prism);
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
				+ coords4[0] + coords5[0]) / 6, (coords0[1] + coords1[1]
				+ coords2[1] + coords3[1] + coords4[1] + coords5[1]) / 6,
				(coords0[2] + coords1[2] + coords2[2] + coords3[2] + coords4[2]
						+ coords5[2]) / 6 };
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
		const emInt *const hexVerts = getHexConn(hex);
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
				+ coords4[0] + coords5[0] + coords6[0] + coords7[0]) / 8,
				(coords0[1] + coords1[1] + coords2[1] + coords3[1] + coords4[1]
						+ coords5[1] + coords6[1] + coords7[1]) / 8, (coords0[2]
						+ coords1[2] + coords2[2] + coords3[2] + coords4[2]
						+ coords5[2] + coords6[2] + coords7[2]) / 8 };
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
	printf("Done computing length scale\n");
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
	if (!sizesOK)
		exit(2);

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

void ExaMesh::prettyPrintCellCount(size_t cells, const char *prefix) const {
	if (cells == 0)
		return;
	printf("%s = ", prefix);
	if (cells >> 30) {
		printf("%.2f B\n", cells / 1.e9);
	} else if (cells >> 20) {
		printf("%.2f M\n", cells / 1.e6);
	} else if (cells >> 10) {
		printf("%.2f K\n", cells / 1.e3);
	} else {
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
	size_t nParts = (outputCells - 1) / maxCellsPerPart + 1;
	if (nParts > numCells)
		nParts = numCells;

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
	size_t ii;
//#pragma omp parallel for schedule(dynamic) reduction(+: totalRefineTime, totalExtractTime, totalTets, totalPyrs, totalPrisms, totalHexes, totalCells) num_threads(8)
	for (ii = 0; ii < nParts; ii++) {
		start = exaTime();
//		char filename[100];
//		sprintf(filename, "/tmp/submesh%03d.vtk", ii);
//		writeVTKFile(filename);
		printf("Part %zu: cells %5d-%5d.\n", ii, parts[ii].getFirst(),
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
		printf("\nCPU time for refinement = %5.2F seconds\n", RS.refineTime);
		printf("                          %5.2F million cells / minute\n",
				(RS.cells / 1000000.) / (RS.refineTime / 60));

		// char filename[100];
		// sprintf(filename, "TestCases/Testsubmesh%03d.vtk", ii);
		// pUM->writeVTKFile(filename);
	}
	printf("\nDone parallel refinement with %zu parts.\n", nParts);
	printf("Time for partitioning:           %10.3F seconds\n", partitionTime);
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
	} else if (totalFileSize >> 30) {
		printf("Total ugrid file size = %.2f GB\n",
				(totalFileSize >> 20) / 1024.);
	} else {
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

// TODO These tests should be in testExa instead.
//void
//ExaMesh:: TestMPI(const emInt &nDivs, const emInt &nParts, ParallelTester* tester,
//const char MeshType)
//{
//	vecPart parts;
//	vecCellPartData vecCPD;
//
//	std::set<int> triRotations;
//	std::set<int> quadRotations;
//
//	vecHashTri tris;
//	vecHashQuad quads;
//
//	std::vector<std::shared_ptr<ExaMesh>> submeshes;
//	vecSharePtrUmesh refinedUMeshes;
//
//	std::vector<TableTri2TableIndex2Index>  matchedTrisAllParts  (nParts);
//	std::vector<TableQuad2TableIndex2Index> matchedQuadsAllParts (nParts);
//
//
//
//	//std::vector< TriFaceVerts , std::unordered_map<emInt,emInt>>  matchedTrisAllParts (nParts);
//	//std::vector< QuadFaceVerts, std::unordered_map<emInt,emInt>>  matchedQuadAllParts (nParts);
//
//
//	partitionCells(this, nParts, parts, vecCPD);
//
//
//	size_t totalTriSize;
//	size_t totalQuadSize;
//
//	this->partFaceMatching(parts, vecCPD, tris, quads,totalTriSize,totalQuadSize);
//
//
//
//
//	for (emInt i = 0; i < nParts; i++)
//	{
//		auto coarse = this->extractCoarseMeshPseudoParallel(parts[i], vecCPD, nDivs, tris[i], quads[i], i);
//		submeshes.emplace_back(coarse.release());
//	}
//
//	assert(submeshes.size() == static_cast<std::size_t>(nParts));
//	// Refine the mesh
//	for (emInt i = 0; i < nParts; i++)
//	{
//
//
//		if(MeshType=='C')
//		{
//			auto inputMesh = dynamic_cast<CubicMesh*> ((submeshes[i].get()));
//					auto refineUmesh = std::make_shared<UMesh>(
//			*(inputMesh), nDivs, i);
//			refinedUMeshes.push_back(refineUmesh);
//		}
//		if(MeshType=='U')
//		{
//			auto inputMesh = dynamic_cast<UMesh*> ((submeshes[i].get()));
//					auto refineUmesh = std::make_shared<UMesh>(
//			*(inputMesh), nDivs, i);
//			refinedUMeshes.push_back(refineUmesh);
//		}
//
//
//	}
//	for (emInt iPart = 0; iPart < nParts; iPart++)
//	{
//		TableTri2TableIndex2Index  matchedTris;
//		TableQuad2TableIndex2Index matchedQuads;
//		exa_set<TriFaceVerts> tri =
//			refinedUMeshes[iPart]->getRefinedPartTris();
//		exa_set<QuadFaceVerts> quads = refinedUMeshes[iPart]->getRefinedPartQuads();
//		for (auto it = tri.begin(); it != tri.end(); it++)
//		{
//			std::unordered_map<emInt, emInt> localRemote;
//			emInt whereToLook = it->getRemoteId();
//			exa_set<TriFaceVerts> remoteTriSet =
//				refinedUMeshes[whereToLook]->getRefinedPartTris();
//			emInt rotation = getTriRotation(*it, remoteTriSet, nDivs);
//			triRotations.insert(rotation);
//
//			matchTri(*it, rotation, nDivs, remoteTriSet, localRemote);
//			matchedTris.emplace(*it,localRemote);
//			//printMatchedTris(localRemote,iPart);
//			for (auto itmap = localRemote.begin();
//				 itmap != localRemote.end(); itmap++)
//			{
//				assert(abs(refinedUMeshes[iPart]->getX(itmap->first) -
//						   refinedUMeshes[whereToLook]->getX(itmap->second)) < TOLTEST);
//				assert(abs(refinedUMeshes[iPart]->getY(itmap->first) -
//						   refinedUMeshes[whereToLook]->getY(itmap->second)) < TOLTEST);
//				assert(abs(refinedUMeshes[iPart]->getZ(itmap->first) -
//						   refinedUMeshes[whereToLook]->getZ(itmap->second)) < TOLTEST);
//			}
//		}
//		for (auto it = quads.begin(); it != quads.end(); it++)
//		{
//			std::unordered_map<emInt, emInt> localRemote;
//			emInt whereToLook = it->getRemoteId();
//			exa_set<QuadFaceVerts> remoteQuadSet =
//				refinedUMeshes[whereToLook]->getRefinedPartQuads();
//			emInt rotation = getQuadRotation(*it, remoteQuadSet, nDivs);
//			quadRotations.insert(rotation);
//
//			matchQuad(*it, rotation, nDivs, remoteQuadSet, localRemote);
//			matchedQuads.emplace(*it,localRemote);
//
//			//printMatchedQuads(localRemote,iPart);
//			for (auto itmap = localRemote.begin();
//				 itmap != localRemote.end(); itmap++)
//			{
//				assert(abs(refinedUMeshes[iPart]->getX(itmap->first) -
//						   refinedUMeshes[whereToLook]->getX(itmap->second)) < TOLTEST);
//				assert(abs(refinedUMeshes[iPart]->getY(itmap->first) -
//						   refinedUMeshes[whereToLook]->getY(itmap->second)) < TOLTEST);
//				assert(abs(refinedUMeshes[iPart]->getZ(itmap->first) -
//						   refinedUMeshes[whereToLook]->getZ(itmap->second)) < TOLTEST);
//			}
//		}
//		matchedTrisAllParts[iPart] = matchedTris;
//		matchedQuadsAllParts[iPart]= matchedQuads;
//	}
//	// Which rotation cases are covered?
//	std::cout << "Covered Tri Rotations: " << std::endl;
//	for (auto it = triRotations.begin();
//		 it != triRotations.end(); it++)
//	{
//		std::cout << *it << " ";
//	};
//	std::cout << std::endl;
//	std::cout << "Covered Quad Rotations: " << std::endl;
//	for (auto it = quadRotations.begin();
//		 it != quadRotations.end(); it++)
//	{
//		std::cout << *it << " ";
//	};
//	std::cout << std::endl;
//
//	tester->setMatchedTris (matchedTrisAllParts);
//	tester->setMatchedQuads(matchedQuadsAllParts);
//}

//void
//ExaMesh::refineForMPI( const int numDivs ,
//ParallelTester* tester,
//const char MeshType,
//std::string mshName,
//FILE* fileAllTimes)
//const
//{
//	boost::mpi::environment   env;
//	boost::mpi::communicator  world;
//
//	double  MAXTotalTime(0);
//	double  MAXSyncTime(0);
//
//	double  MAXExtractionTime(0);
//	double  MAXRefineTime(0);
//	double  MAXTriTime(0);
//	double  MAXQuadTime(0);
//	double  MAXSerialTime(0);
//	double  MAXFaceExchangeTime(0);
//
//	double  StartPartitionTime(0);
//	double  StartPartFaceMatching(0);
//	double  StartExtractionTime(0);
//	double  StartRefineTime(0);
//	double  StartTriTime(0);
//	double  StartQuadTime(0);
//	double  StartSyncTimeForTri(0);
//	double  StartSyncTimeForQuad(0);
//	double  StartFaceExchangeTime(0);
//	//double  StartSerialTime;
//
//	double  PartitionTime(0);
//	double  ExtractionTime(0);
//	double  RefineTime(0);
//	double  TriTime(0);
//	double  QuadTime(0);
//	double  PartFaceMatchingTime(0);
//	double  SyncTimeForTri(0);
//	double  SyncTimeForQuad(0);
//	double  TotalSyncTime(0);
//	double  faceExchangeTime(0);
//
//	double  serialTime(0);
//
//	std::vector<boost::mpi::request>      triReqs;
//	std::vector<boost::mpi::request>      quadReqs;
//	std::vector<std::unique_ptr<UMesh>>   refinedMeshVec;
//
//	vecPart         parts;
//	vecCellPartData vecCPD;
//
//	std::size_t    vecCPDSize;
//	std::size_t    trisSize;
//	std::size_t    quadsSize;
//	std::size_t    nParts = world.size();
//	size_t         nCells(0);
//	size_t         totalCells(0);
//	size_t         totalPartTriSize(0);
//	size_t         totalPartQuadSize(0);
//
//	hashTri  trisS;
//	hashQuad quadsS;
//
//	vecTri   triV;
//	vecQuad  quadV;
//
//	//int MASTER = 0;
//	int tag    = 0;
//
//	intToVecTri  remoteTovecTris;
//	intToVecQuad remoteTovecQuads;
//
//	std::set<int> triNeighbrs;
//	std::set<int> quadNeighbrs;
//
//
//	std::vector<vecTri>  trisTobeRcvd;
//	std::vector<vecQuad> quadsTobeRcvd;
//
//	hashTri  recvdTris;
//	hashQuad recvdQuads;
//
//	TableTri2TableIndex2Index   matchedTris;
//	TableQuad2TableIndex2Index  matchedQuads;
//
//	auto    StartTotalTime = exaTime();
//
//	double StartSerialTime = exaTime();
//
//	if(world.rank()==MASTER)
//	{
//
//		StartPartitionTime= exaTime();
//		partitionCells(this, nParts, parts,vecCPD);
//		PartitionTime= exaTime()- StartPartitionTime;
//
//		vecCPDSize = vecCPD.size();
//
//		assert(vecCPDSize>0);
//
//		vecHashTri  VectrisHash;
//		vecHashQuad VecquadsHash;
//
//		vecVecTri   VecTriVec;
//		vecVecQuad  vecQuadVec;
//
//		StartPartFaceMatching= exaTime();
//		this->partFaceMatching(parts,vecCPD,VectrisHash,VecquadsHash,totalPartTriSize,totalPartQuadSize);
//		PartFaceMatchingTime= exaTime() -StartPartFaceMatching;
//
//
//		for(size_t  itri=0 ; itri<VectrisHash.size(); itri++)
//		{
//			vecTri TriVec;
//			SetToVector(VectrisHash[itri],TriVec);
//			VecTriVec.emplace_back(TriVec);
//		}
//		for(size_t iquad=0 ; iquad<VecquadsHash.size(); iquad++)
//		{
//			vecQuad QuadVec;
//			SetToVector(VecquadsHash[iquad],QuadVec);
//			vecQuadVec.emplace_back(QuadVec);
//		}
//
//		trisS = VectrisHash[0]; // For MASTER
//		quadsS= VecquadsHash[0];
//
//		for(auto irank=1 ; irank<world.size();irank++)
//		{
//			world.send(irank,tag,parts[irank]);
//			world.send(irank,tag,VecTriVec[irank]);
//			world.send(irank,tag,vecQuadVec[irank]);
//		}
//	}
//	else
//	{
//		parts.resize(world.size());
//		world.recv(MASTER,tag,parts[world.rank()]);
//		world.recv(MASTER,tag,triV);
//		world.recv(MASTER,tag,quadV);
//		vectorToSet(triV,trisS);
//		vectorToSet(quadV,quadsS);
//	}
//
//	boost::mpi::broadcast(world,vecCPDSize,MASTER);
//
//	if(world.rank()==MASTER)
//	{
//		for(auto irank=1 ; irank<world.size();irank++)
//		{
//			world.send(irank,tag,vecCPD);
//		}
//	}
//	if(world.rank()!= MASTER)
//	{
//		vecCPD.resize(vecCPDSize);
//		world.recv(MASTER,tag,vecCPD);
//	}
//	serialTime = exaTime()- StartSerialTime;
//	StartExtractionTime = exaTime();
//	auto coarse= this->extractCoarseMeshPseudoParallel(parts[world.rank()],vecCPD,numDivs,
//	 trisS,quadsS,world.rank());
//	ExtractionTime = exaTime()- StartExtractionTime;
//
//
//	StartRefineTime = exaTime();
//	if(MeshType=='C')
//	{
//		auto refinedMesh = coarse->subdivideMesh(numDivs, world.rank());
//		refinedMeshVec.emplace_back(refinedMesh.release());
//	}
//	if(MeshType=='U')
//	{
//		coarse->subdivideMesh(numDivs, world.rank());
//		refinedMeshVec.emplace_back(refinedMesh.release());
//	}
//	RefineTime = exaTime()- StartRefineTime;
//
//	//writeEachRankMeshStatics(world.rank(),refinedMeshVec[0]->numCells(),numDivs,world.size(),
//	//mshName);
//
//	StartFaceExchangeTime = exaTime();
//	auto tris  = refinedMeshVec[0]->getRefinedPartTris();
//	auto quads = refinedMeshVec[0]->getRefinedPartQuads();
//	nCells     = refinedMeshVec[0]->numCells();
//
//
//	buildTrisMap(tris,remoteTovecTris,triNeighbrs);
//	buildQuadsMap(quads,remoteTovecQuads,quadNeighbrs);
//
//
//	trisTobeRcvd.resize (triNeighbrs.size());
//	quadsTobeRcvd.resize(quadNeighbrs.size());
//
//	int Trijj=0 ;
//	for(auto isource:triNeighbrs)
//	{
//		auto findSource = remoteTovecTris.find(isource);
//		// same amount of message will be sent and received from other prcoeszsor
//		// Be cuatious if this assumption later on will be break
//		trisTobeRcvd[Trijj].resize(findSource->second.size());
//		Trijj++;
//	}
//
//	int Quadjj=0;
//	for(auto isource:quadNeighbrs)
//	{
//		auto findSource = remoteTovecQuads.find(isource);
//		// same amount of message will be sent and received from other prcoeszsor
//		// Be cuatious if this assumption later on will be break
//		quadsTobeRcvd[Quadjj].resize(findSource->second.size());
//		Quadjj++;
//	}
//
//	for(auto tri:remoteTovecTris)
//	{
//		boost::mpi::request req= world.isend(tri.first,tag,tri.second);
//		triReqs.push_back(req);
//	}
//
//	for(auto quad:remoteTovecQuads)
//	{
//		boost::mpi::request req = world.isend(quad.first,tag,quad.second);
//		quadReqs.push_back(req);
//	}
//
//	int Quadkk=0 ;
//	for(auto source: quadNeighbrs)
//	{
//		boost::mpi::request req = world.irecv(source,tag,quadsTobeRcvd[Quadkk]);
//		quadReqs.push_back(req);
//		Quadkk++;
//	}
//
//
//	int Trikk=0 ;
//	for(auto source: triNeighbrs)
//	{
//
//		boost::mpi::request req= world.irecv(source,tag,trisTobeRcvd[Trikk]);
//		triReqs.push_back(req);
//		Trikk++;
//	}
//
//	faceExchangeTime = exaTime()-StartFaceExchangeTime;
//
//	StartSyncTimeForTri= exaTime();
//
//	boost::mpi::wait_all(triReqs.begin(),triReqs.end());
//
//	SyncTimeForTri = exaTime()- StartSyncTimeForTri;
//
//	StartTriTime = exaTime();
//	for(const auto& tri: trisTobeRcvd)
//	{
//		// I'm collecting the whole data into a set
//		// I should have received the whole data from it
//
//		recvdTris.insert(tri.begin(),tri.end());
//	}
//
//	for(auto it=tris.begin(); it!=tris.end(); it++)
//	{
//		std::unordered_map<emInt, emInt> localRemote;
//		//int rotation = getTriRotation(*it,recvdTris,numDivs);
//		//matchTri(*it,rotation,numDivs,recvdTris,localRemote);
//		findRotationAndMatchTris(*it,numDivs,recvdTris,localRemote);
//
//#ifndef NDEBUG
//		//matchedTris.emplace(*it,localRemote);
//#endif //NDEBUG
//	}
//	TriTime = exaTime()-StartTriTime;
//
//	StartSyncTimeForQuad= exaTime();
//
//	boost::mpi::wait_all(quadReqs.begin(),quadReqs.end());
//
//	SyncTimeForQuad = exaTime()- StartSyncTimeForQuad;
//
//	StartQuadTime = exaTime();
//	for(const auto& quad:quadsTobeRcvd)
//	{
//		recvdQuads.insert(quad.begin(),quad.end());
//	}
//
//	for(auto it=quads.begin(); it!=quads.end();it++)
//	{
//		std::unordered_map<emInt,emInt> localRemote;
//		//int rotation = getQuadRotation(*it,recvdQuads,numDivs);
//		//matchQuad(*it,rotation,numDivs,recvdQuads,localRemote);
//		findRotationAndMatchQuads(*it,numDivs,recvdQuads,localRemote);
//#ifndef NDEBUG
//		//matchedQuads.emplace(*it,localRemote);
//#endif
//	}
//	QuadTime         = exaTime()- StartQuadTime;
//
//	double Totaltime = exaTime() - StartTotalTime;
//	TotalSyncTime    = SyncTimeForQuad+SyncTimeForTri;
//
//	FILE *TimeEachRank    = fopen((mshName+"TimeEachRank.txt").c_str(), "a");
//
//	setlocale(LC_NUMERIC, "");
//	FILE *meshStaticsRank = fopen((mshName+"staticsEachRank.txt").c_str(), "a");
//
//	writeEachRankMeshStatics(meshStaticsRank,world.rank(),refinedMeshVec[0]->numTets(),
//	refinedMeshVec[0]->numPyramids(),refinedMeshVec[0]->numPrisms(),
//	refinedMeshVec[0]->numHexes(), tris.size(), quads.size(), refinedMeshVec[0]->numCells());
//
//	printTimeEachRank(TimeEachRank,world.rank(),PartitionTime,PartFaceMatchingTime,serialTime,
//	ExtractionTime,RefineTime,faceExchangeTime,TotalSyncTime, TriTime,QuadTime,Totaltime);
//
//	//if(world.rank()==0)
//	//{
//		boost::mpi::reduce(world, ExtractionTime   , MAXExtractionTime ,boost:: mpi::maximum<double>(), MASTER);
//
//		boost::mpi::reduce(world, TotalSyncTime    , MAXSyncTime         ,boost:: mpi::maximum<double>(), MASTER);
//
//		boost::mpi::reduce(world, Totaltime        , MAXTotalTime        ,boost:: mpi::maximum<double>(), MASTER);
//
//		boost::mpi::reduce(world, RefineTime       , MAXRefineTime       ,boost:: mpi::maximum<double>(), MASTER);
//
//		boost::mpi::reduce(world, TriTime          , MAXTriTime          ,boost:: mpi::maximum<double>(), MASTER);
//
//		boost::mpi::reduce(world, QuadTime         , MAXQuadTime         ,boost:: mpi::maximum<double>(), MASTER);
//
//		boost::mpi::reduce(world, serialTime       , MAXSerialTime       ,boost:: mpi::maximum<double>(), MASTER);
//
//		boost::mpi::reduce(world, faceExchangeTime , MAXFaceExchangeTime ,boost:: mpi::maximum<double>(), MASTER);
//	//}
//
//
//
//	boost::mpi::reduce(world, size_t(nCells) , totalCells        , std::plus<size_t>()          , MASTER);
//	mshName = mshName + "-weakScaliblity.txt" ;
//	if (world.rank() == MASTER)
//	{
//		writeAllTimeResults(fileAllTimes,world.size(),PartitionTime,
//		PartFaceMatchingTime,MAXSerialTime,MAXExtractionTime,MAXRefineTime,
//		MAXFaceExchangeTime,MAXSyncTime,MAXTriTime,
//		MAXQuadTime,MAXTotalTime,totalPartTriSize,totalPartQuadSize,totalCells);
//		FILE *fileWeakScaliblity = fopen(mshName.c_str(), "a");
//    	if (fileWeakScaliblity == NULL)
//		{
//        	fprintf(stderr, "Error opening the file!\n");
//   		}
//		writeAllTimeResults(fileWeakScaliblity,world.size(),PartitionTime,
//		PartFaceMatchingTime, MAXSerialTime ,MAXExtractionTime,MAXRefineTime,MAXFaceExchangeTime,
//		MAXSyncTime, MAXTriTime,
//		MAXQuadTime,MAXTotalTime,totalPartTriSize,totalPartQuadSize,totalCells);
//
//
//
//	}
//	//world.barrier();
//	//fprintf(stderr, "Rank of: %d has %lu tris and has %lu quads.\n", world.rank(), tris.size(), quads.size());
//#ifndef NDEBUG
//	//tester->testMatchedTris(matchedTris,world.rank());
//	//tester->testMatchedQuads(matchedQuads,world.rank());
//#endif
//
//}

emInt ExaMesh::FastpartFaceMatching(const emInt nParts,
		const std::vector<std::vector<emInt>> &part2cells,
		const std::vector<emInt> &cell2part, vecVecTri &tris,
		vecVecQuad &quads) const {
	tris.resize(nParts);
	quads.resize(nParts);
	emInt numDivs = 1;

	// Mark the cells that are cell parts in fact
	std::unordered_set<std::pair<emInt, emInt>, pairHash> cellParts;
	for (emInt icell = 0; icell < cell2part.size(); icell++) {
		emInt ipart = cell2part[icell];
		emInt ineighsize = getCellConnSize(icell);
		std::vector<emInt> vneighs;
		for (emInt ineigh = 0; ineigh < ineighsize; ineigh++) {
			vneighs.push_back(getCellConn(icell, ineigh));
		}
		// icell belongs to ipart
		// Check whether all of my neibours are in this part
		for (emInt ineigh = 0; ineigh < ineighsize; ineigh++) {
			emInt neighCellID = vneighs[ineigh];
			emInt neighPartID = cell2part[neighCellID];
			if (neighPartID != ipart) {
				if (neighCellID < icell) {
					cellParts.emplace(neighCellID, icell);
				} else {
					cellParts.emplace(icell, neighCellID);
				}

			}
		}

	}

	for (const auto &icellPart : cellParts) {
		emInt cellID1 = icellPart.first;
		emInt cellID2 = icellPart.second;

		emInt partID1 = cell2part[cellID1];
		emInt partID2 = cell2part[cellID2];

//		printf("Cell 1: (%5u %5u)  Cell 2: (%5u %5u)\n",
//				cellID1, partID1, cellID2, partID2);
//		printf("  Locals: (%5u %5u)  (%5u %5u)\n",
//				cellID2cellTypeLocalID[cellID1].first,
//				cellID2cellTypeLocalID[cellID1].second,
//				cellID2cellTypeLocalID[cellID2].first,
//				cellID2cellTypeLocalID[cellID2].second);

		emInt ind1 = cellID2cellTypeLocalID[cellID1].second;
		emInt ind2 = cellID2cellTypeLocalID[cellID2].second;

		emInt type1 = cellID2cellTypeLocalID[cellID1].first;
		emInt type2 = cellID2cellTypeLocalID[cellID2].first;

		std::vector<TriFaceVerts> tris1, tris2;
		std::vector<QuadFaceVerts> quads1, quads2;

		getFaceLists(ind1, type1, partID1, 1, tris1, quads1);
		getFaceLists(ind2, type2, partID2, 1, tris2, quads2);

		setTri partBdryTris;
		setQuad partBdryQuads;

		for (emInt itri = 0; itri < tris1.size(); itri++) {
			partBdryTris.insert(tris1[itri]);
		}
		for (emInt itri = 0; itri < tris2.size(); itri++) {
			partBdryTris.insert(tris2[itri]);
		}

		for (emInt iquad = 0; iquad < quads1.size(); iquad++) {
			partBdryQuads.insert(quads1[iquad]);
		}

		for (emInt iquad = 0; iquad < quads2.size(); iquad++) {
			partBdryQuads.insert(quads2[iquad]);
		}
		preMatchingPartBdryTris(numDivs, partBdryTris, tris);
		preMatchingPartBdryQuads(numDivs, partBdryQuads, quads);
	}
	return cellParts.size();

}
;

void ExaMesh::getFaceLists(const emInt ind, const emInt type,
		const emInt partID, const emInt numDivs,
		std::vector<TriFaceVerts> &tris,
		std::vector<QuadFaceVerts> &quads) const {
	const emInt *conn;
	switch (type) {
	default:
		// Panic! Should never get here.
		printf("Type: %u  Ind: %u \n", type, ind);
		assert(0);
		break;
	case CGNS_ENUMV(TRI_3):
	case CGNS_ENUMV(TRI_10):
		break;
	case CGNS_ENUMV(QUAD_4):
	case CGNS_ENUMV(QUAD_16):
		break;
	case CGNS_ENUMV(TETRA_4):
	case CGNS_ENUMV(TETRA_20): {

		conn = getTetConn(ind);

		emInt global012[3] = { conn[0], conn[1], conn[2] };
		emInt global013[3] = { conn[0], conn[1], conn[3] };
		emInt global123[3] = { conn[1], conn[2], conn[3] };
		emInt global203[3] = { conn[2], conn[0], conn[3] };
		TriFaceVerts T012(numDivs, global012, partID);
		TriFaceVerts T013(numDivs, global013, partID);
		TriFaceVerts T123(numDivs, global123, partID);
		TriFaceVerts T203(numDivs, global203, partID);
		tris.emplace_back(T012);
		tris.emplace_back(T013);
		tris.emplace_back(T123);
		tris.emplace_back(T203);

		break;
	}
	case CGNS_ENUMV(PYRA_5):
	case CGNS_ENUMV(PYRA_30): {

		conn = getPyrConn(ind);

		emInt global0123[4] = { conn[0], conn[1], conn[2], conn[3] };
		emInt global014[3] = { conn[0], conn[1], conn[4] };
		emInt global124[3] = { conn[1], conn[2], conn[4] };
		emInt global234[3] = { conn[2], conn[3], conn[4] };
		emInt global304[3] = { conn[3], conn[0], conn[4] };
		TriFaceVerts T014(numDivs, global014, partID);
		TriFaceVerts T124(numDivs, global124, partID);
		TriFaceVerts T234(numDivs, global234, partID);
		TriFaceVerts T304(numDivs, global304, partID);
		QuadFaceVerts Q0123(numDivs, global0123, partID);
		tris.emplace_back(T014);
		tris.emplace_back(T124);
		tris.emplace_back(T234);
		tris.emplace_back(T304);
		quads.emplace_back(Q0123);

		break;
	}
	case CGNS_ENUMV(PENTA_6):
	case CGNS_ENUMV(PENTA_40): {

		conn = getPrismConn(ind);

		emInt global0143[4] = { conn[0], conn[1], conn[4], conn[3] };
		emInt global1254[4] = { conn[1], conn[2], conn[5], conn[4] };
		emInt global2035[4] = { conn[2], conn[0], conn[3], conn[5] };

		emInt global012[3] = { conn[0], conn[1], conn[2] };
		emInt global345[3] = { conn[3], conn[4], conn[5] };

		TriFaceVerts T012(numDivs, global012, partID);
		TriFaceVerts T345(numDivs, global345, partID);
		QuadFaceVerts Q0143(numDivs, global0143, partID);
		QuadFaceVerts Q1254(numDivs, global1254, partID);
		QuadFaceVerts Q2035(numDivs, global2035, partID);
		tris.emplace_back(T012);
		tris.emplace_back(T345);
		quads.emplace_back(Q0143);
		quads.emplace_back(Q1254);
		quads.emplace_back(Q2035);
		break;
	}
	case CGNS_ENUMV(HEXA_8):
	case CGNS_ENUMV(HEXA_64): {

		conn = getHexConn(ind);

		emInt global0154[4] = { conn[0], conn[1], conn[5], conn[4] };
		emInt global1265[4] = { conn[1], conn[2], conn[6], conn[5] };
		emInt global2376[4] = { conn[2], conn[3], conn[7], conn[6] };
		emInt global3047[4] = { conn[3], conn[0], conn[4], conn[7] };
		emInt global0123[4] = { conn[0], conn[1], conn[2], conn[3] };
		emInt global4567[4] = { conn[4], conn[5], conn[6], conn[7] };

		QuadFaceVerts Q0154(numDivs, global0154, partID);
		QuadFaceVerts Q1265(numDivs, global1265, partID);
		QuadFaceVerts Q2376(numDivs, global2376, partID);
		QuadFaceVerts Q3047(numDivs, global3047, partID);
		QuadFaceVerts Q0123(numDivs, global0123, partID);
		QuadFaceVerts Q4567(numDivs, global4567, partID);
		quads.emplace_back(Q0154);
		quads.emplace_back(Q1265);
		quads.emplace_back(Q2376);
		quads.emplace_back(Q3047);
		quads.emplace_back(Q0123);
		quads.emplace_back(Q4567);
		break;
	}
	}
}

void ExaMesh::buildCellToCellConnectivity() {

	// Now go through the cell-vert connectivity to set up
	// face-cell connectivity
	{
		multimpFace2Cell face2cell;
		std::vector<emInt> faceVerts;
		std::vector<emInt> sortedTriVerts(3);
		std::vector<emInt> sortedQuadVerts(4);

		cellID2cellTypeLocalID.resize(numCells());

		emInt offset = 0;

		emInt type = getTriType();
		for (emInt ii = 0; ii < numBdryTris(); ii++) {
			const emInt *triConn = getBdryTriConn(ii);
			faceVerts = { triConn[0], triConn[1], triConn[2] };
			sortVerts3(faceVerts.data(), sortedTriVerts.data());
			face2cell.emplace(sortedTriVerts,
					std::make_pair(ii + offset, type));
			cellID2cellTypeLocalID[ii + offset] = std::make_pair(
					type, ii);
		}
		offset += numBdryTris();

		type = getQuadType();
		for (emInt ii = 0; ii < numBdryQuads(); ii++) {
			const emInt *quadConn = getBdryQuadConn(ii);
			faceVerts = { quadConn[0], quadConn[1], quadConn[2], quadConn[3] };
			sortVerts4(faceVerts.data(), sortedQuadVerts.data());
			face2cell.emplace(sortedQuadVerts,
					std::make_pair(ii + offset, type));
			cellID2cellTypeLocalID[ii + offset] = std::make_pair(
					type, ii);
		}
		offset += numBdryQuads();

		type = getTetType();
		for (emInt ii = 0; ii < numTets(); ii++) {
			const emInt *tetConn = getTetConn(ii);
			faceVerts = { tetConn[0], tetConn[1], tetConn[2] };
			sortVerts3(faceVerts.data(), sortedTriVerts.data());
			face2cell.emplace(sortedTriVerts,
					std::make_pair(ii + offset, type));

			faceVerts = { tetConn[0], tetConn[1], tetConn[3] };
			sortVerts3(faceVerts.data(), sortedTriVerts.data());
			face2cell.emplace(sortedTriVerts,
					std::make_pair(ii + offset, type));

			faceVerts = { tetConn[1], tetConn[2], tetConn[3] };
			sortVerts3(faceVerts.data(), sortedTriVerts.data());
			face2cell.emplace(sortedTriVerts,
					std::make_pair(ii + offset, type));

			faceVerts = { tetConn[2], tetConn[0], tetConn[3] };
			sortVerts3(faceVerts.data(), sortedTriVerts.data());
			face2cell.emplace(sortedTriVerts,
					std::make_pair(ii + offset, type));

			cellID2cellTypeLocalID[ii + offset] = std::make_pair(
					type, ii);
		}
		offset += numTets();

		type = getPyrType();
		for (emInt ii = 0; ii < numPyramids(); ii++) {
			const emInt *pyrConn = getPyrConn(ii);
			faceVerts = { pyrConn[0], pyrConn[1], pyrConn[4] };
			sortVerts3(faceVerts.data(), sortedTriVerts.data());
			face2cell.emplace(sortedTriVerts,
					std::make_pair(ii + offset, type));

			faceVerts = { pyrConn[1], pyrConn[2], pyrConn[4] };
			sortVerts3(faceVerts.data(), sortedTriVerts.data());
			face2cell.emplace(sortedTriVerts,
					std::make_pair(ii + offset, type));

			faceVerts = { pyrConn[2], pyrConn[3], pyrConn[4] };
			sortVerts3(faceVerts.data(), sortedTriVerts.data());
			face2cell.emplace(sortedTriVerts,
					std::make_pair(ii + offset, type));

			faceVerts = { pyrConn[3], pyrConn[0], pyrConn[4] };
			sortVerts3(faceVerts.data(), sortedTriVerts.data());
			face2cell.emplace(sortedTriVerts,
					std::make_pair(ii + offset, type));

			faceVerts = { pyrConn[0], pyrConn[1], pyrConn[2], pyrConn[3] };
			sortVerts4(faceVerts.data(), sortedQuadVerts.data());
			face2cell.emplace(sortedQuadVerts,
					std::make_pair(ii + offset, type));

			cellID2cellTypeLocalID[ii + offset] = std::make_pair(
					type, ii);
		}
		offset += numPyramids();

		type = getPrismType();
		for (emInt ii = 0; ii < numPrisms(); ii++) {
			const emInt *prismConn = getPrismConn(ii);
			faceVerts = { prismConn[0], prismConn[1], prismConn[2] };
			sortVerts3(faceVerts.data(), sortedTriVerts.data());
			face2cell.emplace(sortedTriVerts,
					std::make_pair(ii + offset, type));

			faceVerts = { prismConn[3], prismConn[4], prismConn[5] };
			sortVerts3(faceVerts.data(), sortedTriVerts.data());
			face2cell.emplace(sortedTriVerts,
					std::make_pair(ii + offset, type));

			faceVerts =
					{ prismConn[0], prismConn[1], prismConn[4], prismConn[3] };
			sortVerts4(faceVerts.data(), sortedQuadVerts.data());
			face2cell.emplace(sortedQuadVerts,
					std::make_pair(ii + offset, type));

			faceVerts =
					{ prismConn[1], prismConn[2], prismConn[5], prismConn[4] };
			sortVerts4(faceVerts.data(), sortedQuadVerts.data());
			face2cell.emplace(sortedQuadVerts,
					std::make_pair(ii + offset, type));

			faceVerts =
					{ prismConn[2], prismConn[0], prismConn[3], prismConn[5] };
			sortVerts4(faceVerts.data(), sortedQuadVerts.data());
			face2cell.emplace(sortedQuadVerts,
					std::make_pair(ii + offset, type));

			cellID2cellTypeLocalID[ii + offset] = std::make_pair(
					type, ii);
		}
		offset += numPrisms();

		type = getHexType();
		for (emInt ii = 0; ii < numHexes(); ii++) {
			const emInt *hexConn = getHexConn(ii);
			faceVerts = { hexConn[0], hexConn[1], hexConn[2], hexConn[3] };
			sortVerts4(faceVerts.data(), sortedQuadVerts.data());
			face2cell.emplace(sortedQuadVerts,
					std::make_pair(ii + offset, type));

			faceVerts = { hexConn[4], hexConn[5], hexConn[6], hexConn[7] };
			sortVerts4(faceVerts.data(), sortedQuadVerts.data());
			face2cell.emplace(sortedQuadVerts,
					std::make_pair(ii + offset, type));

			faceVerts = { hexConn[0], hexConn[1], hexConn[5], hexConn[4] };
			sortVerts4(faceVerts.data(), sortedQuadVerts.data());
			face2cell.emplace(sortedQuadVerts,
					std::make_pair(ii + offset, type));

			faceVerts = { hexConn[1], hexConn[2], hexConn[6], hexConn[5] };
			sortVerts4(faceVerts.data(), sortedQuadVerts.data());
			face2cell.emplace(sortedQuadVerts,
					std::make_pair(ii + offset, type));

			faceVerts = { hexConn[2], hexConn[3], hexConn[7], hexConn[6] };
			sortVerts4(faceVerts.data(), sortedQuadVerts.data());
			face2cell.emplace(sortedQuadVerts,
					std::make_pair(ii + offset, type));

			faceVerts = { hexConn[3], hexConn[0], hexConn[4], hexConn[7] };
			sortVerts4(faceVerts.data(), sortedQuadVerts.data());
			face2cell.emplace(sortedQuadVerts,
					std::make_pair(ii + offset, type));

			cellID2cellTypeLocalID[ii + offset] = std::make_pair(
					type, ii);
		}

		// Debugging diagnostics
//		for (auto cellPair : cellID2cellTypeLocalID) {
//			printf("ID pairs: %5u %5u\n",
//					cellPair.first, cellPair.second);
//		}
//
//		for (auto faceEntry : face2cell) {
//			printf("Num verts: %lu (%2u %2u %2u %2u)  Cell: %3u  Type: %2u\n",
//					faceEntry.first.size(), faceEntry.first[0],
//					faceEntry.first[1], faceEntry.first[2], faceEntry.first[3],
//					faceEntry.second.first, faceEntry.second.second);
//		}
		printf("Number of face-cell entries: %lu\n", face2cell.size());
		buildCell2CellConnFromFace2Cell(face2cell, numCells());
	}
}

void ExaMesh::buildCell2CellConnFromFace2Cell(multimpFace2Cell &face2cell,
		const emInt nCells) {

	auto it = face2cell.begin();
	vcell2cell.resize(nCells);
	while (it != face2cell.end()) {
		auto current = it;
		auto next = std::next(it);
		if (next != face2cell.end() && current->first == next->first) {
			emInt currentCellID = current->second.first;
			//emInt currentCellType   = current->second.second;

			emInt nextCellID = next->second.first;
			//emInt nextCellType      = next->second.second;

			auto &currentCellVector = vcell2cell[currentCellID];
			auto &nextCellVector = vcell2cell[nextCellID];

			currentCellVector.emplace_back(nextCellID);
			nextCellVector.emplace_back(currentCellID);
		}

		it++;
	}
}

void ExaMesh::testCell2CellConn(emInt nCells) {

	// for(auto icell= cell2cell.begin(); icell!=cell2cell.end();icell++)
	// {
	// 	const int nNeigbrs = icell->second.size();
	// 	assert(nNeigbrs!=0); // number of neighbours can't be zero
	// 	const auto cellType = icell->first.second;
	// 	auto Itrcell =
	// 	cell2faces.find(std::make_pair(icell->first.first,icell->first.second));
	// 	auto const& faces = Itrcell->second;

	// }
	// for (const auto& entry : cell2cell) {
	//     const auto& key = entry.first;
	//     const auto& values = entry.second;

	//     std::cout << "Cell ID: {" << key.first << ", " << key.second << "} Neighbors: {";
	//     for (const auto& value : values) {
	//         std::cout << value << " ";
	//     }
	//     std::cout << "}\n";
	// }

}
void ExaMesh::buidCell2FacesConn(std::pair<emInt, emInt> cellInfo, emInt v0,
		emInt v1, emInt v2) {
	std::set<emInt> faceVerts = { v0, v1, v2 };
	emInt cellID = cellInfo.first;
	emInt cellType = cellInfo.second;
	cell2faces[std::make_pair(cellID, cellType)].insert(faceVerts);
}
void ExaMesh::buidCell2FacesConn(std::pair<emInt, emInt> cellInfo, emInt v0,
		emInt v1, emInt v2, emInt v3) {
	std::set<emInt> faceVerts = { v0, v1, v2, v3 };
	emInt cellID = cellInfo.first;
	emInt cellType = cellInfo.second;
	cell2faces[std::make_pair(cellID, cellType)].insert(faceVerts);
}

void ExaMesh::testCell2FaceConn(emInt nCells) {

	assert(cell2faces.size() == nCells);
	// Checking whethee info of all cells have captured
	// Then, checking whether each cell type has the correct number of faces
	for (auto icell = cell2faces.begin(); icell != cell2faces.end(); icell++) {
		const auto &cellType = (icell->first).second;
		switch (cellType) {
		case CGNS_ENUMV(TETRA_4):
			assert(icell->second.size() == 4);
			break;
		case CGNS_ENUMV(PYRA_5):
			assert(icell->second.size() == 5);
			break;
		case CGNS_ENUMV(PENTA_6):
			assert(icell->second.size() == 5);
			break;
		case CGNS_ENUMV(HEXA_8):
			assert(icell->second.size() == 6);
			break;
			// case CGNS_ENUMV(TRI_3):
			// 	assert(icell->second.size()==3);
			// 	break;
			// case CGNS_ENUMV(QUAD_4):
			// 	assert(icell->second.size()==5);
			// 	break;
			// default:
			// 	assert(0);
		}

	}

}

