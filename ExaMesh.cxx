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
		m_lenScale[ii] = std::numeric_limits<double>::max();
	}
	std::vector<double> vertVolume(numVerts(), std::numeric_limits<double>::max() / 100);
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
