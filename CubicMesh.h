/*
 * CubicMesh.h
 *
 *  Created on: Oct. 1, 2019
 *      Author: cfog
 */

#ifndef SRC_CUBICMESH_H_
#define SRC_CUBICMESH_H_

#include <assert.h>

#include "ExaMesh.h"

// This data structure is organized to read and write easily to/from CGNS files.
class CubicMesh: public ExaMesh {
	emInt m_vert, m_tri, m_quad, m_tet, m_pyr, m_prism, m_hex;
	emInt m_nVerts, m_nBdryVerts, m_nTri10, m_nQuad16, m_nTet20, m_nPyr30,
			m_nPrism40, m_nHex64, m_nVertNodes;
	double *m_xcoords, *m_ycoords, *m_zcoords;
	emInt (*m_Tri10Conn)[10];
	emInt (*m_Quad16Conn)[16];
	emInt (*m_Tet20Conn)[20];
	emInt (*m_Pyr30Conn)[30];
	emInt (*m_Prism40Conn)[40];
	emInt (*m_Hex64Conn)[64];

	CubicMesh(const CubicMesh&);
	CubicMesh& operator=(const CubicMesh&);
	void readCGNSfile(const char CGNSfilename[]);
	void reorderCubicMesh();
	void renumberNodes(emInt thisSize, emInt* aliasConn, emInt* newNodeInd);
	void decrementVertIndices(emInt connSize, emInt* const connect);

	// Length scales
public:
	CubicMesh(const emInt nVerts, const emInt nBdryVerts, const emInt nBdryTris,
			const emInt nBdryQuads, const emInt nTets, const emInt nPyramids,
			const emInt nPrisms, const emInt nHexes);
	CubicMesh(const char CGNSFileName[]);
	virtual ~CubicMesh();

	// Will eventually want to create a fine CubicMesh from a coarse CubicMesh
	// (for the mapping). Not urgent.

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
		return m_nPyr30;
	}
	virtual emInt numPrisms() const {
		return m_nPrism40;
	}
	virtual emInt numHexes() const {
		return m_nHex64;
	}
	virtual emInt numVertsToCopy() const {
		return m_nVertNodes;
	}

	emInt addVert(const double newCoords[3], const emInt coarseGlobalIndex =
	EMINT_MAX);
	emInt addBdryTri(const emInt verts[]);
	emInt addBdryQuad(const emInt verts[]);
	emInt addTet(const emInt verts[]);
	emInt addPyramid(const emInt verts[]);
	emInt addPrism(const emInt verts[]);
	emInt addHex(const emInt verts[]);

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
		assert(pyr < m_nPyr30);
		return m_Pyr30Conn[pyr];
	}

	const emInt* getPrismConn(const emInt prism) const {
		assert(prism < m_nPrism40);
		return m_Prism40Conn[prism];
	}

	const emInt* getHexConn(const emInt hex) const {
		assert(hex < m_nHex64);
		return m_Hex64Conn[hex];
	}

	Mapping::MappingType getDefaultMappingType() const {
		return Mapping::Lagrange;
	}

	std::unique_ptr<CubicMesh> extractCoarseMesh(Part& P,
			std::vector<CellPartData>& vecCPD) const;

	virtual std::unique_ptr<UMesh> createFineUMesh(const emInt numDivs, Part& P,
			std::vector<CellPartData>& vecCPD) const;

	void setupCellDataForPartitioning(std::vector<CellPartData>& vecCPD,
			double &xmin, double& ymin, double& zmin, double& xmax, double& ymax,
			double& zmax) const;

	void setNVertNodes(emInt nVertNodes) {
		m_nVertNodes = nVertNodes;
	}
};



#endif /* SRC_CUBICMESH_H_ */
