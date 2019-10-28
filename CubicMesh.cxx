/*
 * CubicMesh.cxx
 *
 *  Created on: Oct. 3, 2019
 *      Author: cfog
 */

#include <algorithm>
#include <cstdio>
#include <memory>

#include <cgnslib.h>

#include "CubicMesh.h"

#define CHECK_STATUS if (status != CG_OK) cg_error_exit()

CubicMesh::CubicMesh(const emInt nVerts, const emInt nBdryVerts,
		const emInt nBdryTris, const emInt nBdryQuads, const emInt nTets,
		const emInt nPyramids, const emInt nPrisms, const emInt nHexes) :
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
	reorderCubicMesh();
}

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
	fprintf(stderr, "Permuting coordinates: x");
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
	delete cloneCoords;

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

std::unique_ptr<CubicMesh> CubicMesh::extractCoarseMesh(Part& P,
		std::vector<CellPartData>& vecCPD) const {
	auto CUM = std::make_unique<CubicMesh>(0, 0, 0, 0, 0, 0, 0, 0);
	return CUM;
}
