/*
 * CubicMesh.h
 *
 *  Created on: Oct. 1, 2019
 *      Author: cfog
 */

#ifndef SRC_CUBICMESH_H_
#define SRC_CUBICMESH_H_

#include <assert.h>

#include "examesh.h"

// This data structure is organized to read and write easily to/from CGNS files.
class CubicMesh: public ExaMesh {
	emInt m_nVerts, m_nBdryVerts, m_nTri10, m_nQuad16, m_nTet20, m_nPyr29,
			m_nPrism38, m_nHex56;
	double *m_xcoords, *m_ycoords, *m_zcoords;
	emInt (*m_Tri10Conn)[10];
	emInt (*m_Quad16Conn)[16];
	emInt (*m_Tet20Conn)[20];
	emInt (*m_Pyr29Conn)[29];
	emInt (*m_Prism38Conn)[38];
	emInt (*m_Hex56Conn)[56];

	CubicMesh(const CubicMesh&);
	CubicMesh& operator=(const CubicMesh&);

	// Length scales
public:
	CubicMesh(const emInt nVerts, const emInt nBdryVerts, const emInt nBdryTris,
			const emInt nBdryQuads, const emInt nTets, const emInt nPyramids,
			const emInt nPrisms, const emInt nHexes);
	CubicMesh(const char CGNSFileName[]);
	virtual ~CubicMesh();

	// Will eventually want to create a fine CubicMesh from a coarse CubicMesh
	// (for the mapping) and a fine UMesh (for the new cell definitions).
	// Not urgent.

	virtual emInt numVerts() const {
		return m_nVerts;
	}
	virtual emInt numBdryVerts() const {
		return m_nBdryVerts;
	}
	virtual emInt numBdryTris() const {
		return m_nTri10;
	}
	virtual emInt numBdryQuads() const {
		return m_nQuad16;
	}
	virtual emInt numTets() const {
		return m_nTet20;
	}
	virtual emInt numPyramids() const {
		return m_nPyr29;
	}
	virtual emInt numPrisms() const {
		return m_nPrism38;
	}
	virtual emInt numHexes() const {
		return m_nHex56;
	}

	double getX(const emInt vert) const {
		assert(vert < m_nVerts);
		return m_xcoords[vert];
	}

	double getY(const emInt vert) const {
		assert(vert < m_nVerts);
		return m_ycoords[vert];
	}

	double getZ(const emInt vert) const {
		assert(vert < m_nVerts);
		return m_zcoords[vert];
	}

	void getCoords(const emInt vert, double coords[3]) const {
		assert(vert < m_nVerts);
		coords[0] = m_xcoords[vert];
		coords[1] = m_ycoords[vert];
		coords[2] = m_zcoords[vert];
	}

	const emInt* getBdryTriConn(const emInt bdryTri) const {
		assert(bdryTri < m_nTri10);
		return m_Tri10Conn[bdryTri];
	}

	const emInt* getBdryQuadConn(const emInt bdryQuad) const {
		assert(bdryQuad < m_nQuad16);
		return m_Quad16Conn[bdryQuad];
	}

	const emInt* getTetConn(const emInt tet) const {
		assert(tet < m_nTet20);
		return m_Tet20Conn[tet];
	}

	const emInt* getPyrConn(const emInt pyr) const {
		assert(pyr < m_nPyr29);
		return m_Pyr29Conn[pyr];
	}

	const emInt* getPrismConn(const emInt prism) const {
		assert(prism < m_nPrism38);
		return m_Prism38Conn[prism];
	}

	const emInt* getHexConn(const emInt hex) const {
		assert(hex < m_nHex56);
		return m_Hex56Conn[hex];
	}

	Mapping::MappingType getDefaultMappingType() const {
		return Mapping::Lagrange;
	}
};



#endif /* SRC_CUBICMESH_H_ */
