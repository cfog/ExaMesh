/*
 * CubicMesh.cxx
 *
 *  Created on: Oct. 3, 2019
 *      Author: cfog
 */

#include <cstdio>

#include <cgnslib.h>

#include "CubicMesh.h"

#define CHECK_STATUS if (status != CG_OK) cg_error_exit()

CubicMesh::CubicMesh(const emInt nVerts, const emInt nBdryVerts,
		const emInt nBdryTris, const emInt nBdryQuads, const emInt nTets,
		const emInt nPyramids, const emInt nPrisms, const emInt nHexes) :
		m_nVerts(nVerts), m_nBdryVerts(nBdryVerts), m_nTri10(nBdryTris),
				m_nQuad16(nBdryQuads), m_nTet20(nTets), m_nPyr29(nPyramids),
				m_nPrism38(nPrisms), m_nHex56(nHexes) {
	m_xcoords = new double[m_nVerts];
	m_ycoords = new double[m_nVerts];
	m_zcoords = new double[m_nVerts];

	m_Tri10Conn = new emInt[m_nTri10][10];
	m_Quad16Conn = new emInt[m_nQuad16][16];
	m_Tet20Conn = new emInt[m_nTet20][20];
	m_Pyr29Conn = new emInt[m_nPyr29][29];
	m_Prism38Conn = new emInt[m_nPrism38][38];
	m_Hex56Conn = new emInt[m_nHex56][56];
}

CubicMesh::CubicMesh(const char CGNSfilename[]) {
	int status;
	int index_file;

	status = cg_open(CGNSfilename,
	CG_MODE_READ,
										&index_file);
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

	cgsize_t zoneSize[3];
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
		cgsize_t start, end, count;
		int nBdry, parentFlag;
		ElementType_t eType;
		char sectionName[33];
		status = cg_section_read(index_file, 1, 1, iSec, sectionName, &eType,
															&start, &end, &nBdry, &parentFlag);
		CHECK_STATUS;
		count = end - start + 1;

		fprintf(stderr, "Read section %3d (%20s).  %10u elements of type %d.\n",
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
	m_nPyr29 = elementCounts[CGNS_ENUMV(PYRA_29)];
	m_nPrism38 = elementCounts[CGNS_ENUMV(PENTA_38)];
	m_nHex56 = elementCounts[CGNS_ENUMV(HEXA_56)];

	fprintf(stderr, "Allocating space for connectivity.\n");
	// For now, cheat and just allocate for TETRA_20's.
	m_Tri10Conn = nullptr;
	m_Quad16Conn = nullptr;
	m_Tet20Conn = nullptr;
	m_Pyr29Conn = nullptr;
	m_Prism38Conn = nullptr;
	m_Hex56Conn = nullptr;
	// We'd better hope that an emInt and a cgsize_t are compatible.
	m_Tet20Conn = new emInt[elementCounts[CGNS_ENUMV(TETRA_20)]][20];

	fprintf(stderr, "Reading connectivity.\n");
	// Also cheating on the section number.
	// status = cg_elements_read(index_file, 1, 1, 1,
	// 			    reinterpret_cast<cgsize_t*>(tetra4Connect),
	// 			    nullptr);
	status = cg_elements_read(index_file, 1, 1, 1,
														reinterpret_cast<cgsize_t*>(m_Tet20Conn),
														nullptr);

	fprintf(stderr, "Allocating space for coordinates.\n");
	m_xcoords = new double[zoneSize[0]];
	m_ycoords = new double[zoneSize[0]];
	m_zcoords = new double[zoneSize[0]];

	DataType_t dataType = CGNS_ENUMT(RealDouble);
	cgsize_t min = 1, max = zoneSize[0];
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

	fprintf(stderr, "First few points are:\n");
	for (int ii = 0; ii < 10; ii++) {
		fprintf(stderr, "%2d (%.8f %.8f %.8f)\n", ii, m_xcoords[ii], m_ycoords[ii],
						m_zcoords[ii]);
	}

	fprintf(stderr, "First few cells are:\n");
	for (int ii = 0; ii < 10; ii++) {
		fprintf(stderr, "%2d %8d %8d %8d %8d\n", ii, m_Tet20Conn[ii][0],
						m_Tet20Conn[ii][1], m_Tet20Conn[ii][2],
						m_Tet20Conn[ii][3]);
	}

	for (int ii = 0; ii < zoneSize[1]; ii++) {
		for (int jj = 0; jj < 20; jj++) {
			m_Tet20Conn[ii][jj]--;
		}
	}

	fprintf(stderr, "Finding max and min node indices: ");

	cgsize_t maxInd = -1, minInd = 0x7fffffff;
	for (int ii = 0; ii < zoneSize[1]; ii++) {
		for (int jj = 0; jj < 20; jj++) {
			cgsize_t thisInd = m_Tet20Conn[ii][jj];
			if (thisInd < minInd) minInd = thisInd;
			if (thisInd > maxInd) maxInd = thisInd;
		}
	}
	fprintf(stderr, "%d %d\n", minInd, maxInd);
}

CubicMesh::~CubicMesh() {
	delete[] m_xcoords;
	delete[] m_ycoords;
	delete[] m_zcoords;

	delete[] m_Tri10Conn;
	delete[] m_Quad16Conn;
	delete[] m_Tet20Conn;
	delete[] m_Pyr29Conn;
	delete[] m_Prism38Conn;
	delete[] m_Hex56Conn;
}


