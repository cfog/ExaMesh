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
 * UMesh.h
 *
 *  Created on: Jul. 9, 2019
 *      Author: cfog
 */

#ifndef SRC_UMESH_H_
#define SRC_UMESH_H_

#include <assert.h>

#include "CubicMesh.h"
#include "ExaMesh.h"

class UMesh: public ExaMesh {
	emInt m_nVerts, m_nBdryVerts, m_nTris, m_nQuads, m_nTets, m_nPyrs, m_nPrisms,
			m_nHexes;
	enum {
		eVert = 0, eTri, eQuad, eTet, ePyr, ePrism, eHex
	};
	size_t m_fileImageSize;
	emInt *m_header;
	double (*m_coords)[3];
	emInt (*m_TriConn)[3];
	emInt (*m_QuadConn)[4];
	emInt *m_TriBC;
	emInt *m_QuadBC;
	emInt (*m_TetConn)[4];
	emInt (*m_PyrConn)[5];
	emInt (*m_PrismConn)[6];
	emInt (*m_HexConn)[8];
	char *m_buffer, *m_fileImage;
	// exa_set<QuadFaceVerts> TemppartQuads; 
	// exa_set<TriFaceVerts>  TemppartTris; 
	// exa_set<TriFaceVerts>  partTris; 
	// exa_set<TriFaceVerts>  refinedPartTris; 
	// exa_set<QuadFaceVerts> refinedPartQuads; 
	// exa_set<QuadFaceVerts> partQuads; 
	UMesh(const UMesh&);
	UMesh& operator=(const UMesh&);

public:
	void partFaceMatching(
		 std::vector<Part>& parts, const std::vector<CellPartData>& vecCPD,	
		 std::vector<std::unordered_set<TriFaceVerts>>  &tris,
		 std::vector<std::unordered_set<QuadFaceVerts>> &quads )const;		 
	UMesh(const emInt nVerts, const emInt nBdryVerts, const emInt nBdryTris,
			const emInt nBdryQuads, const emInt nTets, const emInt nPyramids,
			const emInt nPrisms, const emInt nHexes);
	UMesh(const char baseFileName[], const char type[], const char ugridInfix[]);
	UMesh(const UMesh& UM_in, const int nDivs, const emInt partID=-1);
	UMesh(const CubicMesh& CM, const int nDivs, const emInt partID=-1);
	~UMesh();
	emInt maxNVerts() const {
		return m_nVerts;
	}
	emInt maxNBdryTris() const {
		return m_nTris;
	}
	emInt maxNBdryQuads() const {
		return m_nQuads;
	}
	emInt maxNTets() const {
		return m_nTets;
	}
	emInt maxNPyrs() const {
		return m_nPyrs;
	}
	emInt maxNPrisms() const {
		return m_nPrisms;
	}
	emInt maxNHexes() const {
		return m_nHexes;
	}
	emInt numVerts() const {
		return m_header[eVert];
	}
	emInt numBdryVerts() const {
		return m_nBdryVerts;
	}
	emInt numBdryTris() const {
		return m_header[eTri];
	}
	emInt numBdryQuads() const {
		return m_header[eQuad];
	}
	emInt numTets() const {
		return m_header[eTet];
	}
	emInt numPyramids() const {
		return m_header[ePyr];
	}
	emInt numPrisms() const {
		return m_header[ePrism];
	}
	emInt numHexes() const {
		return m_header[eHex];
	}
	emInt numCells() const {
		return numTets() + numPyramids() + numPrisms() + numHexes();
	}

	emInt addVert(const double newCoords[3]);
	emInt addBdryTri(const emInt verts[]);
	emInt addBdryQuad(const emInt verts[]);
	emInt addTet(const emInt verts[]);
	emInt addPyramid(const emInt verts[]);
	emInt addPrism(const emInt verts[]);
	emInt addHex(const emInt verts[]);

	virtual void getCoords(const emInt vert, double coords[3]) const {
		assert(vert < m_nVerts && vert < m_header[eVert]);
		const double* const tmp = m_coords[vert];
		coords[0] = tmp[0];
		coords[1] = tmp[1];
		coords[2] = tmp[2];
	}

	double getX(const emInt vert) const {
		assert(vert < m_nVerts && vert < m_header[eVert]);
		return m_coords[vert][0];
	}

	double getY(const emInt vert) const {
		assert(vert < m_nVerts && vert < m_header[eVert]);
		return m_coords[vert][1];
	}

	double getZ(const emInt vert) const {
		assert(vert < m_nVerts && vert < m_header[eVert]);
		return m_coords[vert][2];
	}

	const emInt* getBdryTriConn(const emInt bdryTri) const {
		assert(bdryTri < m_nTris && bdryTri < m_header[eTri]);
		return m_TriConn[bdryTri];
	}

	const emInt* getBdryQuadConn(const emInt bdryQuad) const {
		assert(bdryQuad < m_nQuads && bdryQuad < m_header[eQuad]);
		return m_QuadConn[bdryQuad];
	}

	const emInt* getTetConn(const emInt tet) const {
		assert(tet < m_nTets && tet < m_header[eTet]);
		return m_TetConn[tet];
	}

	const emInt* getPyrConn(const emInt pyr) const {
		assert(pyr < m_nPyrs && pyr < m_header[ePyr]);
		return m_PyrConn[pyr];
	}

	const emInt* getPrismConn(const emInt prism) const {
		assert(prism < m_nPrisms && prism < m_header[ePrism]);
		return m_PrismConn[prism];
	}

	const emInt* getHexConn(const emInt hex) const {
		assert(hex < m_nHexes && hex < m_header[eHex]);
		return m_HexConn[hex];
	}

	Mapping::MappingType getDefaultMappingType() const {
		return Mapping::Uniform;
	}

	virtual std::unique_ptr<UMesh> createFineUMesh(const emInt numDivs, Part& P,
			std::vector<CellPartData>& vecCPD, struct RefineStats& RS) const;


	virtual std::unique_ptr<UMesh> extractCoarseMesh(Part& P,	std::vector<CellPartData>& vecCPD, 
	const int numDivs,
			const std::unordered_set<TriFaceVerts> &tris= std::unordered_set<TriFaceVerts>(), 
			const std::unordered_set<QuadFaceVerts> &quads= std::unordered_set<QuadFaceVerts>(), 
			const emInt partID=-1) const;


			
	virtual void TestMPI(const emInt &nDivs, const emInt &nParts);

	void refineForMPI(ParallelTester* const tester, const int numDivs) const;

	void setupCellDataForPartitioning(std::vector<CellPartData>& vecCPD,
			double &xmin, double& ymin, double& zmin, double& xmax, double& ymax,
			double& zmax) const;

	bool writeVTKFile(const char fileName[]);
	bool writeUGridFile(const char fileName[]);

	size_t getFileImageSize() const {
		return m_fileImageSize;
	}

	void incrementVertIndices(emInt* conn, emInt size, int inc);

	// Writing with compression reduces file size by a little over a factor of two,
	// at the expense of making file write slower by two orders of magnitude.
	// So don't do it.
	// bool writeCompressedUGridFile(const char fileName[]);
private:
	void init(const emInt nVerts, const emInt nBdryVerts, const emInt nBdryTris,
			const emInt nBdryQuads, const emInt nTets, const emInt nPyramids,
			const emInt nPrisms, const emInt nHexes);
};




#endif /* SRC_UMESH_H_ */
