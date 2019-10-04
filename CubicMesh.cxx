/*
 * CubicMesh.cxx
 *
 *  Created on: Oct. 3, 2019
 *      Author: cfog
 */

#include "CubicMesh.h"

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

CubicMesh::CubicMesh(const char CGNSFileName[]) {
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


