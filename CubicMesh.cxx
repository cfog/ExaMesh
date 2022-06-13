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
 * CubicMesh.cxx
 *
 *  Created on: Oct. 3, 2019
 *      Author: cfog
 */

#include <algorithm>
#include <cstdio>
#include <memory>
#include <stdlib.h>

#include "exa_config.h"

#if (HAVE_CGNS == 1)
#include <cgnslib.h>

#define CHECK_STATUS if (status != CG_OK) cg_error_exit()
#endif

#include "CubicMesh.h"
#include "UMesh.h"

CubicMesh::CubicMesh(const emInt nVerts, const emInt nBdryVerts,
		const emInt nBdryTris, const emInt nBdryQuads, const emInt nTets,
		const emInt nPyramids, const emInt nPrisms, const emInt nHexes) :
		m_vert(0), m_tri(0), m_quad(0), m_tet(0), m_pyr(0), m_prism(0), m_hex(0),
				m_nVerts(nVerts), m_nBdryVerts(nBdryVerts), m_nTri10(nBdryTris),
				m_nQuad16(nBdryQuads), m_nTet20(nTets), m_nPyr30(nPyramids),
				m_nPrism40(nPrisms), m_nHex64(nHexes), m_nVertNodes(0) {
	m_xcoords = new double[m_nVerts];
	m_ycoords = new double[m_nVerts];
	m_zcoords = new double[m_nVerts];

	m_Tri10Conn = new emInt[m_nTri10][10];
	m_Quad16Conn = new emInt[m_nQuad16][16];
	m_Tet20Conn = new emInt[m_nTet20][20];
	m_Pyr30Conn = new emInt[m_nPyr30][30];
	m_Prism40Conn = new emInt[m_nPrism40][40];
	m_Hex64Conn = new emInt[m_nHex64][64];
}

void CubicMesh::decrementVertIndices(emInt connSize, emInt* const connect) {
	for (emInt ii = 0; ii < connSize; ii++) {
		connect[ii]--;
	}
}

#if (HAVE_CGNS == 1)
void CubicMesh::readCGNSfile(const char CGNSfilename[]) {
	int status;
	int index_file;
	status = cg_open(CGNSfilename, CG_MODE_READ, &index_file);
	CHECK_STATUS;
	fprintf(stderr, "Opened CGNS file %s\n", CGNSfilename);
	int nBases = -1;
	status = cg_nbases(index_file, &nBases);
	CHECK_STATUS;
	if (nBases != 1) {
		fprintf(stderr, "Can only handle one base node\n");
		exit(1);
	}
	int topoDim = -1, geomDim = -1;
	char baseName[33];
	status = cg_base_read(index_file, 1, baseName, &topoDim, &geomDim);
	CHECK_STATUS;
	fprintf(stderr, "Got base node %s, dims %d/%d\n", baseName, topoDim, geomDim);
	int nZones = -1;
	status = cg_nzones(index_file, 1, &nZones);
	CHECK_STATUS;
	if (nZones != 1) {
		fprintf(stderr, "Can only handle one zone\n");
		exit(1);
	}
	ZoneType_t zoneType;
	status = cg_zone_type(index_file, 1, 1, &zoneType);
	CHECK_STATUS;
	if (zoneType != CGNS_ENUMV(Unstructured)) {
		fprintf(stderr, "Bad zone type %d\n", zoneType);
		exit(1);
	}
	int zoneSize[3];
	char zoneName[33];
	status = cg_zone_read(index_file, 1, 1, zoneName, zoneSize);
	CHECK_STATUS;
	fprintf(stderr, "Got zone %s, size: %u verts, %u cells\n", zoneName,
					zoneSize[0], zoneSize[1]);
	m_nVerts = zoneSize[0];
	m_nBdryVerts = zoneSize[2];
	int nSections = -1;
	cg_nsections(index_file, 1, 1, &nSections);
	CHECK_STATUS;
	fprintf(stderr, "Zone has %d sections\n", nSections);
	// Now parse the sections to see how many of what kind of elements there are.
	size_t elementCounts[CGNS_ENUMV(HEXA_125) + 1] = { 0 };
	for (int iSec = 1; iSec <= nSections; iSec++) {
		int start, end, count;
		int nBdry, parentFlag;
		ElementType_t eType;
		char sectionName[33];
		status = cg_section_read(index_file, 1, 1, iSec, sectionName, &eType,
															&start, &end, &nBdry, &parentFlag);
		CHECK_STATUS;
		count = end - start + 1;
		fprintf(stderr, "Scanned section %3d (%20s).  %10u elements of type %d.\n",
						iSec, sectionName, count, eType);
		elementCounts[eType] += count;
	}
	for (int type = 0; type <= CGNS_ENUMV(HEXA_125); type++) {
		if (elementCounts[type] != 0) {
			fprintf(stderr, "%10lu elements of type %d.\n", elementCounts[type],
							type);
		}
	}
	m_nTri10 = elementCounts[CGNS_ENUMV(TRI_10)];
	m_nQuad16 = elementCounts[CGNS_ENUMV(QUAD_16)];
	m_nTet20 = elementCounts[CGNS_ENUMV(TETRA_20)];
	m_nPyr30 = elementCounts[CGNS_ENUMV(PYRA_30)];
	m_nPrism40 = elementCounts[CGNS_ENUMV(PENTA_40)];
	m_nHex64 = elementCounts[CGNS_ENUMV(HEXA_64)];
	fprintf(stderr, "Allocating space for connectivity.\n");
	// We'd better hope that an emInt and a cgsize_t are compatible.
	m_Tri10Conn = new emInt[m_nTri10][10];
	m_Quad16Conn = new emInt[m_nQuad16][16];
	m_Tet20Conn = new emInt[m_nTet20][20];
	m_Pyr30Conn = new emInt[m_nPyr30][30];
	m_Prism40Conn = new emInt[m_nPrism40][40];
	m_Hex64Conn = new emInt[m_nHex64][64];
	fprintf(stderr, "Reading connectivity.\n");
	emInt tri10count = 0, quad16count = 0, tet20count = 0, pyr30count = 0,
			prism40count = 0, hex64count = 0;
	for (int iSec = 1; iSec <= nSections; iSec++) {
		int start, end, count, totalCount;
		int nBdry, parentFlag;
		ElementType_t eType;
		char sectionName[33];
		status = cg_section_read(index_file, 1, 1, iSec, sectionName, &eType,
															&start, &end, &nBdry, &parentFlag);
		CHECK_STATUS;
		count = end - start + 1;
		switch (eType) {
			case CGNS_ENUMV(TRI_10):
				status = cg_elements_partial_read(
						index_file, 1, 1, iSec, start, end,
						reinterpret_cast<int*>(m_Tri10Conn + tri10count), nullptr);
				totalCount = tri10count += count;
				break;
			case CGNS_ENUMV(QUAD_16):
				status = cg_elements_partial_read(
						index_file, 1, 1, iSec, start, end,
						reinterpret_cast<int*>(m_Quad16Conn + quad16count), nullptr);
				totalCount = quad16count += count;
				break;
			case CGNS_ENUMV(TETRA_20):
				status = cg_elements_partial_read(
						index_file, 1, 1, iSec, start, end,
						reinterpret_cast<int*>(m_Tet20Conn + tet20count), nullptr);
				totalCount = tet20count += count;
				break;
			case CGNS_ENUMV(PYRA_30):
				status = cg_elements_partial_read(
						index_file, 1, 1, iSec, start, end,
						reinterpret_cast<int*>(m_Pyr30Conn + pyr30count), nullptr);
				totalCount = pyr30count += count;
				break;
			case CGNS_ENUMV(PENTA_40):
				status = cg_elements_partial_read(
						index_file, 1, 1, iSec, start, end,
						reinterpret_cast<int*>(m_Prism40Conn + prism40count), nullptr);
				totalCount = prism40count += count;
				break;
			case CGNS_ENUMV(HEXA_64):
				status = cg_elements_partial_read(
						index_file, 1, 1, iSec, start, end,
						reinterpret_cast<int*>(m_Hex64Conn + hex64count), nullptr);
				totalCount = hex64count += count;
				break;
			default:
				assert(0);
				break;
		}
		fprintf(stderr,
						"Read section %3d (%20s).  %10u elements of type %d, %10u total.\n",
						iSec, sectionName, count, eType, totalCount);
		fflush(stderr);
	}
	fprintf(stderr, "Allocating space for coordinates.\n");
	m_xcoords = new double[zoneSize[0]];
	m_ycoords = new double[zoneSize[0]];
	m_zcoords = new double[zoneSize[0]];
	DataType_t dataType = CGNS_ENUMT(RealDouble);
	int min = 1, max = zoneSize[0];
	status = cg_coord_read(index_file, 1, 1, "CoordinateX", dataType, &min, &max,
													m_xcoords);
	CHECK_STATUS;
	fprintf(stderr, "Read x coords, min %u, max %u\n", min, max);
	status = cg_coord_read(index_file, 1, 1, "CoordinateY", dataType, &min, &max,
													m_ycoords);
	CHECK_STATUS;
	fprintf(stderr, "Read y coords, min %u, max %u\n", min, max);
	status = cg_coord_read(index_file, 1, 1, "CoordinateZ", dataType, &min, &max,
													m_zcoords);
	CHECK_STATUS;
	fprintf(stderr, "Read z coords, min %u, max %u\n", min, max);
	status = cg_close(index_file);
	CHECK_STATUS;
//	fprintf(stderr, "First few points are:\n");
//	for (int ii = 0; ii < 10; ii++) {
//		fprintf(stderr, "%2d (%.8f %.8f %.8f)\n", ii, m_xcoords[ii], m_ycoords[ii],
//						m_zcoords[ii]);
//	}
//	fprintf(stderr, "First few cells are:\n");
//	for (int ii = 0; ii < 10; ii++) {
//		fprintf(stderr, "%2d %8d %8d %8d %8d\n", ii, m_Tet20Conn[ii][0],
//						m_Tet20Conn[ii][1], m_Tet20Conn[ii][2], m_Tet20Conn[ii][3]);
//	}

// Need to decrement all the connectivity, because CGNS is 1-based, and
// everything else in the non-Fortran parts of the universe is 0-based.
// With a little casting, this can be factored and easy.
	decrementVertIndices(10 * m_nTri10, reinterpret_cast<emInt*>(m_Tri10Conn));
	decrementVertIndices(16 * m_nQuad16, reinterpret_cast<emInt*>(m_Quad16Conn));
	decrementVertIndices(20 * m_nTet20, reinterpret_cast<emInt*>(m_Tet20Conn));
	decrementVertIndices(30 * m_nPyr30, reinterpret_cast<emInt*>(m_Pyr30Conn));
	decrementVertIndices(40 * m_nPrism40,
												reinterpret_cast<emInt*>(m_Prism40Conn));
	decrementVertIndices(64 * m_nHex64, reinterpret_cast<emInt*>(m_Hex64Conn));

	fprintf(stderr, "Counting bdry verts: ");
	// Now tag all bdry verts
	bool *isBdryVert = new bool[m_nVerts];
	for (emInt ii = 0; ii < m_nVerts; ii++) {
		isBdryVert[ii] = false;
	}
	for (emInt iTri = 0; iTri < m_nTri10; iTri++) {
		isBdryVert[m_Tri10Conn[iTri][0]] = true;
		isBdryVert[m_Tri10Conn[iTri][1]] = true;
		isBdryVert[m_Tri10Conn[iTri][2]] = true;
	}
	for (emInt iQuad = 0; iQuad < m_nQuad16; iQuad++) {
		isBdryVert[m_Quad16Conn[iQuad][0]] = true;
		isBdryVert[m_Quad16Conn[iQuad][1]] = true;
		isBdryVert[m_Quad16Conn[iQuad][2]] = true;
		isBdryVert[m_Quad16Conn[iQuad][3]] = true;
	}
	m_nBdryVerts = 0;
	for (emInt ii = 0; ii < m_nVerts; ii++) {
		if (isBdryVert[ii]) {
			m_nBdryVerts++;
		}
	}
	delete[] isBdryVert;
	fprintf(stderr, "%'u\n", m_nBdryVerts);
}

CubicMesh::CubicMesh(const char CGNSfilename[]) {
	readCGNSfile(CGNSfilename);
	m_vert = m_nVerts;
	m_tri = m_nTri10;
	m_quad = m_nQuad16;
	m_tet = m_nTet20;
	m_pyr = m_nPyr30;
	m_prism = m_nPrism40;
	m_hex = m_nHex64;
	reorderCubicMesh();
//	buildFaceCellConnectivity();
}
#endif

void CubicMesh::renumberNodes(emInt thisSize, emInt* aliasConn,
		emInt* newNodeInd) {
	emInt* cloneConn = new emInt[thisSize];
	std::copy(aliasConn, aliasConn + thisSize, cloneConn);
	for (emInt ii = 0; ii < thisSize; ii++) {
		aliasConn[ii] = newNodeInd[cloneConn[ii]];
	}
	delete[] cloneConn;
}

void CubicMesh::reorderCubicMesh() {
	emInt *newNodeInd = new emInt[m_nVerts];
	bool *isVertexNode = new bool[m_nVerts];

	for (emInt ii = 0; ii < m_nVerts; ii++) {
		newNodeInd[ii] = EMINT_MAX;
		isVertexNode[ii] = false;
	}

	fprintf(stderr, "Tagging vertex nodes.\n");
	for (emInt ii = 0; ii < m_nTet20; ii++) {
		for (emInt jj = 0; jj < 4; jj++) {
			emInt node = m_Tet20Conn[ii][jj];
			isVertexNode[node] = true;
		}
	}
	for (emInt ii = 0; ii < m_nPyr30; ii++) {
		for (emInt jj = 0; jj < 5; jj++) {
			emInt node = m_Pyr30Conn[ii][jj];
			isVertexNode[node] = true;
		}
	}
	for (emInt ii = 0; ii < m_nPrism40; ii++) {
		for (emInt jj = 0; jj < 6; jj++) {
			emInt node = m_Prism40Conn[ii][jj];
			isVertexNode[node] = true;
		}
	}
	for (emInt ii = 0; ii < m_nHex64; ii++) {
		for (emInt jj = 0; jj < 8; jj++) {
			emInt node = m_Hex64Conn[ii][jj];
			isVertexNode[node] = true;
		}
	}

	fprintf(stderr, "Renumbering vertex nodes\n");
	emInt node = 0;
	for (emInt ii = 0; ii < m_nVerts; ii++) {
		if (isVertexNode[ii]) {
			assert(newNodeInd[ii] == EMINT_MAX);
			newNodeInd[ii] = node;
			node++;
		}
	}
	m_nVertNodes = node;
	fprintf(stderr, "%'u vertex nodes.\nRenumbering other nodes.\n", node);
	for (emInt ii = 0; ii < m_nVerts; ii++) {
		if (!isVertexNode[ii]) {
			assert(newNodeInd[ii] == EMINT_MAX);
			newNodeInd[ii] = node;
			node++;
		}
	}
	assert(node == m_nVerts);
	delete[] isVertexNode;

	// Clone and re-order the coordinates.  This is an out-of-place
	// re-ordering, but I'm about to go nuts on memory anyway, so who cares
	// about the overhead?
	fprintf(stderr, "Permuting coordinate storage order: x");
	double *cloneCoords = new double[m_nVerts];
	std::copy(m_xcoords, m_xcoords + m_nVerts, cloneCoords);
	for (emInt ii = 0; ii < m_nVerts; ii++) {
		m_xcoords[newNodeInd[ii]] = cloneCoords[ii];
	}
	fprintf(stderr, " y");
	std::copy(m_ycoords, m_ycoords + m_nVerts, cloneCoords);
	for (emInt ii = 0; ii < m_nVerts; ii++) {
		m_ycoords[newNodeInd[ii]] = cloneCoords[ii];
	}
	fprintf(stderr, " z");
	std::copy(m_zcoords, m_zcoords + m_nVerts, cloneCoords);
	for (emInt ii = 0; ii < m_nVerts; ii++) {
		m_zcoords[newNodeInd[ii]] = cloneCoords[ii];
	}
	fprintf(stderr, "\n");
	delete[] cloneCoords;

	// Now update the connectivity.  This is done with a "mapping" to a 1D
	// array to reduce loop overhead and, at least as importantly, make it so
	// that I can factor the whole thing.
	renumberNodes(10 * m_nTri10, reinterpret_cast<emInt*>(m_Tri10Conn),
								newNodeInd);
	renumberNodes(16 * m_nQuad16, reinterpret_cast<emInt*>(m_Quad16Conn),
								newNodeInd);
	renumberNodes(20 * m_nTet20, reinterpret_cast<emInt*>(m_Tet20Conn),
								newNodeInd);
	renumberNodes(30 * m_nPyr30, reinterpret_cast<emInt*>(m_Pyr30Conn),
								newNodeInd);
	renumberNodes(40 * m_nPrism40, reinterpret_cast<emInt*>(m_Prism40Conn),
								newNodeInd);
	renumberNodes(64 * m_nHex64, reinterpret_cast<emInt*>(m_Hex64Conn),
								newNodeInd);

	delete[] newNodeInd;
}

CubicMesh::~CubicMesh() {
	delete[] m_xcoords;
	delete[] m_ycoords;
	delete[] m_zcoords;

	delete[] m_Tri10Conn;
	delete[] m_Quad16Conn;
	delete[] m_Tet20Conn;
	delete[] m_Pyr30Conn;
	delete[] m_Prism40Conn;
	delete[] m_Hex64Conn;
}

static void remapIndices(const emInt nPts, const emInt newIndices[],
		const emInt* conn, emInt* newConn) {
	for (emInt jj = 0; jj < nPts; jj++) {
		newConn[jj] = newIndices[conn[jj]];
	}
}

std::unique_ptr<CubicMesh> CubicMesh::extractCoarseMesh(Part& P,
		std::vector<CellPartData>& vecCPD, const int numDivs) const {
	CALLGRIND_TOGGLE_COLLECT
	;

	// Count the number of tris, quads, tets, pyrs, prisms and hexes.
	const emInt first = P.getFirst();
	const emInt last = P.getLast();

	exa_set<TriFaceVerts> partBdryTris;
	exa_set<QuadFaceVerts> partBdryQuads;

	emInt nTris(0), nQuads(0), nTets(0), nPyrs(0), nPrisms(0), nHexes(0);
	const emInt *conn;

//	std::vector<bool> isVertUsed(numVerts(), false);
	bool *isVertUsed = reinterpret_cast<bool*>(calloc(numVerts(), sizeof(bool)));
	exa_set<emInt> bdryVerts;
	exa_set<emInt> cornerNodes;

	for (emInt ii = first; ii < last; ii++) {
		emInt type = vecCPD[ii].getCellType();
		emInt ind = vecCPD[ii].getIndex();
		switch (type) {
			default:
				// Panic! Should never get here.
				assert(0);
				break;
			case TETRA_20: {
				nTets++;
				conn = getTetConn(ind);
				TriFaceVerts TFV012(numDivs, conn[0], conn[1], conn[2], TETRA_20, ind);
				TriFaceVerts TFV013(numDivs, conn[0], conn[1], conn[3], TETRA_20, ind);
				TriFaceVerts TFV123(numDivs, conn[1], conn[2], conn[3], TETRA_20, ind);
				TriFaceVerts TFV203(numDivs, conn[2], conn[0], conn[3], TETRA_20, ind);
				addUniquely(partBdryTris, TFV012);
				addUniquely(partBdryTris, TFV013);
				addUniquely(partBdryTris, TFV123);
				addUniquely(partBdryTris, TFV203);
//				vertsUsed.insert(conn, conn + 20);
				for (int jj = 0; jj < 20; jj++) {
					isVertUsed[conn[jj]] = true;
				}
				cornerNodes.insert(conn, conn + 4);
				break;
			}
			case PYRA_30: {
				nPyrs++;
				conn = getPyrConn(ind);
				QuadFaceVerts QFV0123(numDivs, conn[0], conn[1], conn[2], conn[3], PYRA_30, ind);
				TriFaceVerts TFV014(numDivs, conn[0], conn[1], conn[4], PYRA_30, ind);
				TriFaceVerts TFV124(numDivs, conn[1], conn[2], conn[4], PYRA_30, ind);
				TriFaceVerts TFV234(numDivs, conn[2], conn[3], conn[4], PYRA_30, ind);
				TriFaceVerts TFV304(numDivs, conn[3], conn[0], conn[4], PYRA_30, ind);
				addUniquely(partBdryQuads, QFV0123);
				addUniquely(partBdryTris, TFV014);
				addUniquely(partBdryTris, TFV124);
				addUniquely(partBdryTris, TFV234);
				addUniquely(partBdryTris, TFV304);
//				vertsUsed.insert(conn, conn + 30);
				for (int jj = 0; jj < 30; jj++) {
					isVertUsed[conn[jj]] = true;
				}
				cornerNodes.insert(conn, conn + 5);
				break;
			}
			case PENTA_40: {
				nPrisms++;
				conn = getPrismConn(ind);
				QuadFaceVerts QFV0143(numDivs, conn[0], conn[1], conn[4], conn[3], PENTA_40,
															ind);
				QuadFaceVerts QFV1254(numDivs, conn[1], conn[2], conn[5], conn[4], PENTA_40,
															ind);
				QuadFaceVerts QFV2035(numDivs, conn[2], conn[0], conn[3], conn[5], PENTA_40,
															ind);
				TriFaceVerts TFV012(numDivs, conn[0], conn[1], conn[2], PENTA_40, ind);
				TriFaceVerts TFV345(numDivs, conn[3], conn[4], conn[5], PENTA_40, ind);
				addUniquely(partBdryQuads, QFV0143);
				addUniquely(partBdryQuads, QFV1254);
				addUniquely(partBdryQuads, QFV2035);
				addUniquely(partBdryTris, TFV012);
				addUniquely(partBdryTris, TFV345);
//				vertsUsed.insert(conn, conn + 40);

				for (int jj = 0; jj < 40; jj++) {
					isVertUsed[conn[jj]] = true;
				}
				cornerNodes.insert(conn, conn + 6);
				break;
			}
			case HEXA_64: {
				nHexes++;
				conn = getHexConn(ind);
				QuadFaceVerts QFV0154(numDivs, conn[0], conn[1], conn[5], conn[4], HEXA_64, ind);
				QuadFaceVerts QFV1265(numDivs, conn[1], conn[2], conn[6], conn[5], HEXA_64, ind);
				QuadFaceVerts QFV2376(numDivs, conn[2], conn[3], conn[7], conn[6], HEXA_64, ind);
				QuadFaceVerts QFV3047(numDivs, conn[3], conn[0], conn[6], conn[7], HEXA_64, ind);
				QuadFaceVerts QFV0123(numDivs, conn[0], conn[1], conn[2], conn[3], HEXA_64, ind);
				QuadFaceVerts QFV4567(numDivs, conn[4], conn[5], conn[6], conn[7], HEXA_64, ind);
				addUniquely(partBdryQuads, QFV0154);
				addUniquely(partBdryQuads, QFV1265);
				addUniquely(partBdryQuads, QFV2376);
				addUniquely(partBdryQuads, QFV3047);
				addUniquely(partBdryQuads, QFV0123);
				addUniquely(partBdryQuads, QFV4567);
//				vertsUsed.insert(conn, conn + 64);
				for (int jj = 0; jj < 64; jj++) {
					isVertUsed[conn[jj]] = true;
				}
				cornerNodes.insert(conn, conn + 8);
				break;
			}
		} // end switch
	} // end loop to gather information

	// Now check to see which bdry entities are in this part.  That'll be the
	// ones whose verts are all marked as used.  Unfortunately, this requires
	// searching through -all- the bdry entities for each part.
	std::vector<emInt> realBdryTris;
	std::vector<emInt> realBdryQuads;
	for (emInt ii = 0; ii < numBdryTris(); ii++) {
		conn = getBdryTriConn(ii);
		if (isVertUsed[conn[0]] && isVertUsed[conn[1]] && isVertUsed[conn[2]]) {
			TriFaceVerts TFV(numDivs, conn[0], conn[1], conn[2]);
			auto iter = partBdryTris.find(TFV);
			// If this bdry tri is an unmatched tri from this part, match it, and
			// add the bdry tri to the list of things to copy to the part coarse
			// mesh.  Otherwise, do nothing.  This will keep the occasional wrong
			// bdry face from slipping through.
			if (iter != partBdryTris.end()) {
				partBdryTris.erase(iter);
				bdryVerts.insert(conn[0]);
				bdryVerts.insert(conn[1]);
				bdryVerts.insert(conn[2]);
				realBdryTris.push_back(ii);
				nTris++;
			}
		}
	}
	for (emInt ii = 0; ii < numBdryQuads(); ii++) {
		conn = getBdryQuadConn(ii);
		if (isVertUsed[conn[0]] && isVertUsed[conn[1]] && isVertUsed[conn[2]]
				&& isVertUsed[conn[3]]) {
			QuadFaceVerts QFV(numDivs, conn[0], conn[1], conn[2], conn[3]);
			auto iter = partBdryQuads.find(QFV);
			// If this bdry tri is an unmatched tri from this part, match it, and
			// add the bdry tri to the list of things to copy to the part coarse
			// mesh.  Otherwise, do nothing.  This will keep the occasional wrong
			// bdry face from slipping through.
			if (iter != partBdryQuads.end()) {
				partBdryQuads.erase(iter);
				bdryVerts.insert(conn[0]);
				bdryVerts.insert(conn[1]);
				bdryVerts.insert(conn[2]);
				bdryVerts.insert(conn[3]);
				realBdryQuads.push_back(ii);
				nQuads++;
			}
		}
	}

	emInt nPartBdryTris = partBdryTris.size();
	emInt nPartBdryQuads = partBdryQuads.size();

	for (auto tri : partBdryTris) {
		bdryVerts.insert(tri.getCorner(0));
		bdryVerts.insert(tri.getCorner(1));
		bdryVerts.insert(tri.getCorner(2));
	}
	for (auto quad : partBdryQuads) {
		bdryVerts.insert(quad.getCorner(0));
		bdryVerts.insert(quad.getCorner(1));
		bdryVerts.insert(quad.getCorner(2));
		bdryVerts.insert(quad.getCorner(3));
	}
	emInt nBdryVerts = 0, nNodes = 0;
	emInt nVertNodes = 0;
	emInt nVerts = numVerts();
	for (emInt ii = 0; ii < nVerts; ii++) {
		if (isVertUsed[ii]) nNodes++;
	}
	nBdryVerts = bdryVerts.size();
	nVertNodes = cornerNodes.size();

	// Now set up the data structures for the new coarse UMesh
	auto UCM = std::make_unique<CubicMesh>(nNodes, nBdryVerts,
																					nTris + nPartBdryTris,
																					nQuads + nPartBdryQuads, nTets, nPyrs,
																					nPrisms, nHexes);
	UCM->setNVertNodes(nVertNodes);

	// Store the vertices, while keeping a mapping from the full list of verts
	// to the restricted list so the connectivity can be copied properly.
	emInt *newIndices = new emInt[nVerts];
	for (emInt ii = 0; ii < nVerts; ii++) {
		if (isVertUsed[ii]) {
			double coords[3];
			getCoords(ii, coords);
			newIndices[ii] = UCM->addVert(coords);
			// Copy length scale for vertices from the parent; otherwise, there will be
			// mismatches in the refined meshes.
			UCM->setLengthScale(newIndices[ii], getLengthScale(ii));
		}
	}
	free(isVertUsed);

	// Now copy connectivity.
	emInt newConn[64];
	for (emInt ii = first; ii < last; ii++) {
		emInt type = vecCPD[ii].getCellType();
		emInt ind = vecCPD[ii].getIndex();
		switch (type) {
			default:
				// Panic! Should never get here.
				assert(0);
				break;
			case TETRA_20: {
				conn = getTetConn(ind);
				remapIndices(20, newIndices, conn, newConn);
				UCM->addTet(newConn);
				break;
			}
			case PYRA_30: {
				conn = getPyrConn(ind);
				remapIndices(30, newIndices, conn, newConn);
				UCM->addPyramid(newConn);
				break;
			}
			case PENTA_40: {
				conn = getPrismConn(ind);
				remapIndices(40, newIndices, conn, newConn);
				UCM->addPrism(newConn);
				break;
			}
			case HEXA_64: {
				conn = getHexConn(ind);
				remapIndices(64, newIndices, conn, newConn);
				UCM->addHex(newConn);
				break;
			}
		} // end switch
	} // end loop to copy most connectivity

	for (emInt ii = 0; ii < realBdryTris.size(); ii++) {
		conn = getBdryTriConn(realBdryTris[ii]);
		remapIndices(10, newIndices, conn, newConn);
		UCM->addBdryTri(newConn);
	}
	for (emInt ii = 0; ii < realBdryQuads.size(); ii++) {
		conn = getBdryQuadConn(realBdryQuads[ii]);
		remapIndices(16, newIndices, conn, newConn);
		UCM->addBdryQuad(newConn);
	}

	// Now, finally, the part bdry connectivity.
	// TODO: Currently, there's nothing in the data structure that marks which
	// are part bdry faces.
	for (auto tri : partBdryTris) {
		emInt cellInd = tri.getVolElement();
		emInt conn[10];
		switch (tri.getVolElementType()) {
			case TETRA_20: {
				emInt *elemConn = m_Tet20Conn[cellInd];
				// Identify which face this is.  Has to be 012, 013, 123, or 203.
				if (tri.getCorner(2) == elemConn[2]) {
					// Has to be 012
					assert(tri.getCorner(0) == elemConn[0]);
					assert(tri.getCorner(1) == elemConn[1]);
					conn[0] = elemConn[0];
					conn[1] = elemConn[1];
					conn[2] = elemConn[2];
					conn[3] = elemConn[4];
					conn[4] = elemConn[5];
					conn[5] = elemConn[6];
					conn[6] = elemConn[7];
					conn[7] = elemConn[8];
					conn[8] = elemConn[9];
					conn[9] = elemConn[16];
				}
				else if (tri.getCorner(0) == elemConn[0]) {
					// Has to be 013
					assert(tri.getCorner(1) == elemConn[1]);
					assert(tri.getCorner(2) == elemConn[3]);
					conn[0] = elemConn[0];
					conn[1] = elemConn[1];
					conn[2] = elemConn[3];
					// Between 0 and 1
					conn[3] = elemConn[4];
					conn[4] = elemConn[5];
					// Between 1 and 3
					conn[5] = elemConn[12];
					conn[6] = elemConn[13];
					// Between 3 and 0
					conn[7] = elemConn[11];
					conn[8] = elemConn[10];
					// On face
					conn[9] = elemConn[17];
				}
				else if (tri.getCorner(0) == elemConn[1]) {
					// Has to be 123
					assert(tri.getCorner(1) == elemConn[2]);
					assert(tri.getCorner(2) == elemConn[3]);
					conn[0] = elemConn[1];
					conn[1] = elemConn[2];
					conn[2] = elemConn[3];
					// Between 1 and 2
					conn[3] = elemConn[6];
					conn[4] = elemConn[7];
					// Between 2 and 3
					conn[5] = elemConn[14];
					conn[6] = elemConn[15];
					// Between 3 and 1
					conn[7] = elemConn[13];
					conn[8] = elemConn[12];
					// On face
					conn[9] = elemConn[18];
				}
				else if (tri.getCorner(0) == elemConn[2]) {
					// Has to be 203
					assert(tri.getCorner(1) == elemConn[0]);
					assert(tri.getCorner(2) == elemConn[3]);
					conn[0] = elemConn[2];
					conn[1] = elemConn[0];
					conn[2] = elemConn[3];
					// Between 2 and 0
					conn[3] = elemConn[8];
					conn[4] = elemConn[9];
					// Between 0 and 3
					conn[5] = elemConn[10];
					conn[6] = elemConn[11];
					// Between 3 and 2
					conn[7] = elemConn[15];
					conn[8] = elemConn[14];
					// On face
					conn[9] = elemConn[19];
				}
				else {
					// Should never get here
					assert(0);
				}
				break;
			}
			case PYRA_30: {
				emInt *elemConn = m_Pyr30Conn[cellInd];
				if (tri.getCorner(0) == elemConn[0]) {
					assert(tri.getCorner(1) == elemConn[1]);
					assert(tri.getCorner(2) == elemConn[4]);
					conn[0] = elemConn[0];
					conn[1] = elemConn[1];
					conn[2] = elemConn[4];
					// Between 0 and 1
					conn[3] = elemConn[5];
					conn[4] = elemConn[6];
					// Between 1 and 4
					conn[5] = elemConn[15];
					conn[6] = elemConn[16];
					// Between 4 and 0
					conn[7] = elemConn[14];
					conn[8] = elemConn[13];
					// On face
					conn[9] = elemConn[25];
				}
				else if (tri.getCorner(0) == elemConn[1]) {
					assert(tri.getCorner(1) == elemConn[2]);
					assert(tri.getCorner(2) == elemConn[4]);
					conn[0] = elemConn[1];
					conn[1] = elemConn[2];
					conn[2] = elemConn[4];
					// Between 1 and 2
					conn[3] = elemConn[7];
					conn[4] = elemConn[8];
					// Between 2 and 4
					conn[5] = elemConn[17];
					conn[6] = elemConn[18];
					// Between 4 and 1
					conn[7] = elemConn[16];
					conn[8] = elemConn[15];
					// On face
					conn[9] = elemConn[26];
				}
				else if (tri.getCorner(0) == elemConn[2]) {
					assert(tri.getCorner(1) == elemConn[3]);
					assert(tri.getCorner(2) == elemConn[4]);
					conn[0] = elemConn[2];
					conn[1] = elemConn[3];
					conn[2] = elemConn[4];
					// Between 2 and 3
					conn[3] = elemConn[9];
					conn[4] = elemConn[10];
					// Between 3 and 4
					conn[5] = elemConn[19];
					conn[6] = elemConn[20];
					// Between 4 and 1
					conn[7] = elemConn[18];
					conn[8] = elemConn[17];
					// On face
					conn[9] = elemConn[27];
				}
				else if (tri.getCorner(0) == elemConn[3]) {
					assert(tri.getCorner(1) == elemConn[0]);
					assert(tri.getCorner(2) == elemConn[4]);
					conn[0] = elemConn[3];
					conn[1] = elemConn[0];
					conn[2] = elemConn[4];
					// Between 3 and 0
					conn[3] = elemConn[13];
					conn[4] = elemConn[14];
					// Between 0 and 4
					conn[5] = elemConn[17];
					conn[6] = elemConn[18];
					// Between 4 and 3
					conn[7] = elemConn[20];
					conn[8] = elemConn[19];
					// On face
					conn[9] = elemConn[28];
				}
				else {
					// Should never get here
					assert(0);
				}
				break;
			}
			case PENTA_40: {
				emInt *elemConn = m_Prism40Conn[cellInd];
				if (tri.getCorner(0)  == elemConn[0]) {
					assert(tri.getCorner(1) == elemConn[1]);
					assert(tri.getCorner(2) == elemConn[2]);
					conn[0] = elemConn[0];
					conn[1] = elemConn[1];
					conn[2] = elemConn[2];
					// Between 0 and 1
					conn[3] = elemConn[6];
					conn[4] = elemConn[7];
					// Between 1 and 2
					conn[5] = elemConn[8];
					conn[6] = elemConn[9];
					// Between 2 and 0
					conn[7] = elemConn[10];
					conn[8] = elemConn[11];
					// On face
					conn[9] = elemConn[24];
				}
				else if (tri.getCorner(0) == elemConn[3]) {
					assert(tri.getCorner(1) == elemConn[4]);
					assert(tri.getCorner(2) == elemConn[5]);
					conn[0] = elemConn[3];
					conn[1] = elemConn[4];
					conn[2] = elemConn[5];
					// Between 3 and 4
					conn[3] = elemConn[18];
					conn[4] = elemConn[19];
					// Between 4 and 5
					conn[5] = elemConn[20];
					conn[6] = elemConn[21];
					// Between 5 and 3
					conn[7] = elemConn[22];
					conn[8] = elemConn[23];
					// On face
					conn[9] = elemConn[37];
				}
				else {
					// Should never get here
					assert(0);
				}
				break;
			}
			default: {
				// Should never get here.
				assert(0);
			}
		}
		emInt newConn[10];
		remapIndices(10, newIndices, conn, newConn);
		UCM->addBdryTri(newConn);
	}

	for (auto quad : partBdryQuads) {
		emInt cellInd = quad.getVolElement();
		emInt conn[16];
		switch (quad.getVolElementType()) {
			case PYRA_30: {
				// Only one quad here, so it had better be the right one.
				emInt *elemConn = m_Pyr30Conn[cellInd];
				assert(quad.getCorner(0) == elemConn[0]);
				assert(quad.getCorner(1) == elemConn[1]);
				assert(quad.getCorner(2) == elemConn[2]);
				assert(quad.getCorner(3) == elemConn[3]);

				conn[0] = elemConn[0];
				conn[1] = elemConn[1];
				conn[2] = elemConn[2];
				conn[3] = elemConn[3];
				// Between 0 and 1
				conn[4] = elemConn[5];
				conn[5] = elemConn[6];
				// Between 1 and 2
				conn[6] = elemConn[7];
				conn[7] = elemConn[8];
				// Between 2 and 3
				conn[8] = elemConn[9];
				conn[9] = elemConn[10];
				// Between 3 and 0
				conn[10] = elemConn[11];
				conn[11] = elemConn[12];
				// On face
				conn[12] = elemConn[21];
				conn[13] = elemConn[22];
				conn[14] = elemConn[23];
				conn[15] = elemConn[24];
				break;
			}
			case PENTA_40: {
				emInt *elemConn = m_Prism40Conn[cellInd];

				// Three possible quads: 0143 1254 2035
				if (quad.getCorner(0) == elemConn[0]) {
					// 0143
					assert(quad.getCorner(1) == elemConn[1]);
					assert(quad.getCorner(2) == elemConn[4]);
					assert(quad.getCorner(3) == elemConn[3]);

					conn[0] = elemConn[0];
					conn[1] = elemConn[1];
					conn[2] = elemConn[4];
					conn[3] = elemConn[3];
					// Between 0 and 1
					conn[4] = elemConn[6];
					conn[5] = elemConn[7];
					// Between 1 and 4
					conn[6] = elemConn[14];
					conn[7] = elemConn[15];
					// Between 4 and 3
					conn[8] = elemConn[19];
					conn[9] = elemConn[18];
					// Between 3 and 0
					conn[10] = elemConn[13];
					conn[11] = elemConn[12];
					// On face
					conn[12] = elemConn[25];
					conn[13] = elemConn[26];
					conn[14] = elemConn[27];
					conn[15] = elemConn[28];
				}
				else if (quad.getCorner(0) == elemConn[1]) {
					// 1254
					assert(quad.getCorner(1) == elemConn[2]);
					assert(quad.getCorner(2) == elemConn[5]);
					assert(quad.getCorner(3) == elemConn[4]);

					conn[0] = elemConn[1];
					conn[1] = elemConn[2];
					conn[2] = elemConn[5];
					conn[3] = elemConn[4];
					// Between 1 and 2
					conn[4] = elemConn[8];
					conn[5] = elemConn[9];
					// Between 2 and 5
					conn[6] = elemConn[16];
					conn[7] = elemConn[17];
					// Between 5 and 4
					conn[8] = elemConn[21];
					conn[9] = elemConn[20];
					// Between 4 and 1
					conn[10] = elemConn[15];
					conn[11] = elemConn[14];
					// On face
					conn[12] = elemConn[29];
					conn[13] = elemConn[30];
					conn[14] = elemConn[31];
					conn[15] = elemConn[32];
				}
				else if (quad.getCorner(0) == elemConn[2]) {
					// 2035
					assert(quad.getCorner(1) == elemConn[0]);
					assert(quad.getCorner(2) == elemConn[3]);
					assert(quad.getCorner(3) == elemConn[5]);

					conn[0] = elemConn[2];
					conn[1] = elemConn[0];
					conn[2] = elemConn[3];
					conn[3] = elemConn[5];
					// Between 2 and 0
					conn[4] = elemConn[10];
					conn[5] = elemConn[11];
					// Between 0 and 3
					conn[6] = elemConn[12];
					conn[7] = elemConn[13];
					// Between 3 and 5
					conn[8] = elemConn[23];
					conn[9] = elemConn[22];
					// Between 5 and 0
					conn[10] = elemConn[17];
					conn[11] = elemConn[16];
					// On face
					conn[12] = elemConn[33];
					conn[13] = elemConn[34];
					conn[14] = elemConn[35];
					conn[15] = elemConn[36];
				}
				else {
					// Should never get here
					assert(0);
				}
				break;
			}
			case HEXA_64: {
				emInt *elemConn = m_Hex64Conn[cellInd];

				// Six quads: 0154 1265 2376 3047 0123 4567
				if (quad.getCorner(2) == elemConn[2]) {
					// Bottom: 0123
					assert(quad.getCorner(0) == elemConn[0]);
					assert(quad.getCorner(1) == elemConn[1]);
					assert(quad.getCorner(3) == elemConn[3]);

					conn[0] = elemConn[0];
					conn[1] = elemConn[1];
					conn[2] = elemConn[2];
					conn[3] = elemConn[3];
					// Between 0 and 1
					conn[4] = elemConn[8];
					conn[5] = elemConn[9];
					// Between 1 and 2
					conn[6] = elemConn[10];
					conn[7] = elemConn[11];
					// Between 2 and 3
					conn[8] = elemConn[12];
					conn[9] = elemConn[13];
					// Between 3 and 0
					conn[10] = elemConn[14];
					conn[11] = elemConn[15];
					// On face
					conn[12] = elemConn[32];
					conn[13] = elemConn[33];
					conn[14] = elemConn[34];
					conn[15] = elemConn[35];
				}
				else if (quad.getCorner(0) == elemConn[4]) {
					// Top: 4567
					assert(quad.getCorner(1) == elemConn[5]);
					assert(quad.getCorner(2) == elemConn[6]);
					assert(quad.getCorner(3) == elemConn[7]);

					conn[0] = elemConn[4];
					conn[1] = elemConn[5];
					conn[2] = elemConn[6];
					conn[3] = elemConn[7];
					// Between 4 and 5
					conn[4] = elemConn[24];
					conn[5] = elemConn[25];
					// Between 5 and 6
					conn[6] = elemConn[26];
					conn[7] = elemConn[27];
					// Between 6 and 7
					conn[8] = elemConn[28];
					conn[9] = elemConn[29];
					// Between 7 and 4
					conn[10] = elemConn[30];
					conn[11] = elemConn[31];
					// On face
					conn[12] = elemConn[52];
					conn[13] = elemConn[53];
					conn[14] = elemConn[54];
					conn[15] = elemConn[55];
				}
				else if (quad.getCorner(0) == elemConn[0]) {
					// Side: 0154
					assert(quad.getCorner(1) == elemConn[1]);
					assert(quad.getCorner(2) == elemConn[5]);
					assert(quad.getCorner(3) == elemConn[4]);

					conn[0] = elemConn[0];
					conn[1] = elemConn[1];
					conn[2] = elemConn[5];
					conn[3] = elemConn[4];
					// Between 0 and 1
					conn[4] = elemConn[8];
					conn[5] = elemConn[9];
					// Between 1 and 5
					conn[6] = elemConn[18];
					conn[7] = elemConn[19];
					// Between 5 and 4
					conn[8] = elemConn[25];
					conn[9] = elemConn[24];
					// Between 4 and 1
					conn[10] = elemConn[17];
					conn[11] = elemConn[16];
					// On face
					conn[12] = elemConn[36];
					conn[13] = elemConn[37];
					conn[14] = elemConn[38];
					conn[15] = elemConn[39];
				}
				else if (quad.getCorner(0) == elemConn[1]) {
					// Side: 1254
					assert(quad.getCorner(1) == elemConn[2]);
					assert(quad.getCorner(2) == elemConn[6]);
					assert(quad.getCorner(3) == elemConn[5]);

					conn[0] = elemConn[1];
					conn[1] = elemConn[2];
					conn[2] = elemConn[6];
					conn[3] = elemConn[5];
					// Between 1 and 2
					conn[4] = elemConn[10];
					conn[5] = elemConn[11];
					// Between 2 and 6
					conn[6] = elemConn[20];
					conn[7] = elemConn[21];
					// Between 6 and 5
					conn[8] = elemConn[27];
					conn[9] = elemConn[26];
					// Between 5 and 1
					conn[10] = elemConn[19];
					conn[11] = elemConn[18];
					// On face
					conn[12] = elemConn[40];
					conn[13] = elemConn[41];
					conn[14] = elemConn[42];
					conn[15] = elemConn[43];
				}
				else if (quad.getCorner(0) == elemConn[2]) {
					// Side: 2376
					assert(quad.getCorner(1) == elemConn[3]);
					assert(quad.getCorner(2) == elemConn[7]);
					assert(quad.getCorner(3) == elemConn[6]);

					conn[0] = elemConn[2];
					conn[1] = elemConn[3];
					conn[2] = elemConn[7];
					conn[3] = elemConn[6];
					// Between 2 and 3
					conn[4] = elemConn[12];
					conn[5] = elemConn[13];
					// Between 3 and 7
					conn[6] = elemConn[22];
					conn[7] = elemConn[23];
					// Between 7 and 6
					conn[8] = elemConn[29];
					conn[9] = elemConn[28];
					// Between 6 and 2
					conn[10] = elemConn[21];
					conn[11] = elemConn[20];
					// On face
					conn[12] = elemConn[44];
					conn[13] = elemConn[45];
					conn[14] = elemConn[46];
					conn[15] = elemConn[47];
				}
				else if (quad.getCorner(0) == elemConn[3]) {
					// Side: 3047
					assert(quad.getCorner(1) == elemConn[0]);
					assert(quad.getCorner(2) == elemConn[4]);
					assert(quad.getCorner(3) == elemConn[7]);

					conn[0] = elemConn[3];
					conn[1] = elemConn[0];
					conn[2] = elemConn[4];
					conn[3] = elemConn[7];
					// Between 3 and 0
					conn[4] = elemConn[14];
					conn[5] = elemConn[15];
					// Between 0 and 4
					conn[6] = elemConn[16];
					conn[7] = elemConn[17];
					// Between 4 and 7
					conn[8] = elemConn[31];
					conn[9] = elemConn[30];
					// Between 7 and 3
					conn[10] = elemConn[23];
					conn[11] = elemConn[22];
					// On face
					conn[12] = elemConn[48];
					conn[13] = elemConn[49];
					conn[14] = elemConn[50];
					conn[15] = elemConn[51];
				}
				else {
					// Should never get here
					assert(0);
				}

				break;
			}
			default:
				// Should never get here.
				assert(0);
		}
		emInt newConn[10];
		remapIndices(10, newIndices, conn, newConn);
		UCM->addBdryQuad(newConn);
	}
	delete[] newIndices;
	CALLGRIND_TOGGLE_COLLECT
	;
	return UCM;
}

emInt CubicMesh::addVert(const double newCoords[3]) {
	assert(m_vert < m_nVerts);
	m_xcoords[m_vert] = newCoords[0];
	m_ycoords[m_vert] = newCoords[1];
	m_zcoords[m_vert] = newCoords[2];
	return (m_vert++);
}

emInt CubicMesh::addBdryTri(const emInt verts[]) {
	assert(m_tri < m_nTri10);
	for (int ii = 0; ii < 10; ii++) {
		assert(verts[ii] < m_nVerts);
		m_Tri10Conn[m_tri][ii] = verts[ii];
	}
	return (m_tri++);
}

emInt CubicMesh::addBdryQuad(const emInt verts[]) {
	assert(m_quad < m_nQuad16);
	for (int ii = 0; ii < 16; ii++) {
		assert(verts[ii] < m_nVerts);
		m_Quad16Conn[m_quad][ii] = verts[ii];
	}
	return (m_quad++);
}

emInt CubicMesh::addTet(const emInt verts[]) {
	assert(m_tet < m_nTet20);
	for (int ii = 0; ii < 20; ii++) {
		assert(verts[ii] < m_nVerts);
		m_Tet20Conn[m_tet][ii] = verts[ii];
	}
	return (m_tet++);
}

emInt CubicMesh::addPyramid(const emInt verts[]) {
	assert(m_pyr < m_nPyr30);
	for (int ii = 0; ii < 30; ii++) {
		assert(verts[ii] < m_nVerts);
		m_Pyr30Conn[m_pyr][ii] = verts[ii];
	}
	return (m_pyr++);
}

emInt CubicMesh::addPrism(const emInt verts[]) {
	assert(m_prism < m_nPrism40);
	for (int ii = 0; ii < 40; ii++) {
		assert(verts[ii] < m_nVerts);
		m_Prism40Conn[m_prism][ii] = verts[ii];
	}
	return (m_prism++);
}

emInt CubicMesh::addHex(const emInt verts[]) {
	assert(m_hex < m_nHex64);
	for (int ii = 0; ii < 64; ii++) {
		assert(verts[ii] < m_nVerts);
		m_Hex64Conn[m_hex][ii] = verts[ii];
	}
	return (m_hex++);
}

std::unique_ptr<UMesh> CubicMesh::createFineUMesh(const emInt numDivs, Part& P,
		std::vector<CellPartData>& vecCPD, struct RefineStats& RS) const {
	// Create a coarse
	double start = exaTime();
	auto coarse = extractCoarseMesh(P, vecCPD, numDivs);
	double middle = exaTime();
	RS.extractTime = middle - start;

	// For some reason, I needed the helper variable to keep the compiler happy here.
	auto UUM = std::make_unique<UMesh>(*coarse, numDivs);
	RS.cells = UUM->numCells();
	RS.refineTime = exaTime() - middle;
	return UUM;
}

void CubicMesh::setupCellDataForPartitioning(std::vector<CellPartData>& vecCPD,
		double &xmin, double& ymin, double& zmin, double& xmax, double& ymax,
		double& zmax) const {
	// Partitioning only cells, not bdry faces.  Also, currently no
	// cost differential for different cell types.
	for (emInt ii = 0; ii < numTets(); ii++) {
		const emInt* verts = getTetConn(ii);
		addCellToPartitionData(verts, 20, ii, TETRA_20, vecCPD, xmin, ymin, zmin,
														xmax, ymax, zmax);
	}
	for (emInt ii = 0; ii < numPyramids(); ii++) {
		const emInt* verts = getPyrConn(ii);
		addCellToPartitionData(verts, 30, ii, PYRA_30, vecCPD, xmin, ymin, zmin,
														xmax, ymax, zmax);
	}
	for (emInt ii = 0; ii < numPrisms(); ii++) {
		const emInt* verts = getPrismConn(ii);
		addCellToPartitionData(verts, 40, ii, PENTA_40, vecCPD, xmin, ymin, zmin,
														xmax, ymax, zmax);
	}
	for (emInt ii = 0; ii < numHexes(); ii++) {
		const emInt* verts = getHexConn(ii);
		addCellToPartitionData(verts, 64, ii, HEXA_64, vecCPD, xmin, ymin, zmin,
														xmax, ymax, zmax);
	}
}

