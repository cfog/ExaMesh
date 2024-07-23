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
 * CubicMesh.h
 *
 *  Created on: Oct. 1, 2019
 *      Author: cfog
 */

#ifndef SRC_CUBICMESH_H_
#define SRC_CUBICMESH_H_

#include <assert.h>

#include "exa-defs.h"
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
	exa_set <TriFaceVerts> TemppartTris;;
	exa_set <QuadFaceVerts> TemppartQuads; 
	exa_set<TriFaceVerts>  partTris; 
	exa_set<QuadFaceVerts> partQuads;
	exa_set<TriFaceVerts>  refinedPartTris;

	CubicMesh(const CubicMesh&);
	CubicMesh& operator=(const CubicMesh&);
	void readCGNSfile(const char CGNSfilename[]);
	void reorderCubicMesh();
	void renumberNodes(emInt thisSize, emInt* aliasConn, emInt* newNodeInd);
	void decrementVertIndices(emInt connSize, emInt* const connect);

	// Confirm positive volume for all subelements
	bool verifyTetValidity() const;
	bool verifyPyramidValidity() const;
	bool verifyPrismValidity() const;
	bool verifyHexValidity() const;

	// Length scales
public:
	CubicMesh(const emInt nVerts, const emInt nBdryVerts, const emInt nBdryTris,
			const emInt nBdryQuads, const emInt nTets, const emInt nPyramids,
			const emInt nPrisms, const emInt nHexes);
#if (HAVE_CGNS == 1)
	CubicMesh(const char CGNSFileName[]);
#endif
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
	virtual emInt numCells() const {
		// This needs to be changed and having number of boundary quads and tris as well 
		return numTets() + numPyramids() + numPrisms() + numHexes();
	}
	virtual emInt numVertsToCopy() const {
		return m_nVertNodes;
	}
	emInt numBdryTrisFromReader()  const
	{
		// return m_nTrisFromReader;
		// it needs to be implemented 
	} 
	emInt numBdryQuadsFromReader() const
	{
		// return m_nQuadsFromReader;
		// it needs to be implemented 
	}

	emInt addVert(const double newCoords[3]);
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
	// TODO ; NOT SET FOR CUBIC MESH 
	// emInt getSizePartTris()const{
	// 	return TemppartTris.size();
	// }
	// emInt getSizePartQuads()const{
	// 	return TemppartQuads.size();
	// }
	// exa_set <QuadFaceVerts> getTempQuadPart() const{
	// 	return TemppartQuads; 
	// }
	// exa_set <TriFaceVerts> getTempTriPart() const {
	// 	return TemppartTris; 
	// }
	// exa_set<QuadFaceVerts> getQuadPart() const{
	// 	return partQuads; 
	// }
	// exa_set<TriFaceVerts> getTriPart() const {
	// 	return partTris; 
	// }
	// exa_set<TriFaceVerts> getRefinedPartTris() const{
	// 	return refinedPartTris; 
	// }
	void partFaceMatching(
		 std::vector<Part>& parts, const std::vector<CellPartData>& vecCPD,	
		 std::vector<std::unordered_set<TriFaceVerts>>  &tris,
		 std::vector<std::unordered_set<QuadFaceVerts>> &quads,size_t &totalTriSize, size_t &totalQuadSize )const;

	std::unique_ptr<ExaMesh> extractCoarseMesh(Part& P,
			std::vector<CellPartData>& vecCPD, const int numDivs,			
			const std::unordered_set<TriFaceVerts> &tris= std::unordered_set<TriFaceVerts>(), 
			const std::unordered_set<QuadFaceVerts> &quads= std::unordered_set<QuadFaceVerts>(), 
			const emInt partID=-1) const;

	virtual std::unique_ptr<UMesh> createFineUMesh(const emInt numDivs, Part& P,
			std::vector<CellPartData>& vecCPD, struct RefineStats& RS) const;
	void setupCellDataForPartitioning(std::vector<CellPartData>& vecCPD,
			double &xmin, double& ymin, double& zmin, double& xmax, double& ymax,
			double& zmax) const;

	void setNVertNodes(emInt nVertNodes) {
		m_nVertNodes = nVertNodes;
	}
};

#endif /* SRC_CUBICMESH_H_ */
