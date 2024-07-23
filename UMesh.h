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
			m_nHexes, m_nTrisFromReader, m_nQuadsFromReader;
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
	TableCell2Cell                                                  cell2cell; 
	std::map < std::pair<emInt,emInt>, std::set<std::set<emInt>>>   cell2faces; 
	std::unordered_map<emInt, std::set<std::set<emInt>>>            cell2bdryfaces; 
	UMesh& operator=(const UMesh& inPut);

public:
	void partFaceMatching(
		 std::vector<Part>& parts, const std::vector<CellPartData>& vecCPD,	
		 std::vector<std::unordered_set<TriFaceVerts>>  &tris,
		 std::vector<std::unordered_set<QuadFaceVerts>> &quads, size_t &totalTriSize, size_t &totalQuadSize)const;		 
	UMesh(const emInt nVerts, const emInt nBdryVerts, const emInt nBdryTris,
			const emInt nBdryQuads, const emInt nTets, const emInt nPyramids,
			const emInt nPrisms, const emInt nHexes);
	UMesh(const emInt nVerts, const emInt nBdryVerts, const emInt nBdryTris,
			const emInt nBdryQuads, const emInt nTets, const emInt nPyramids,
			const emInt nPrisms, const emInt nHexes, 
			const std::vector<std::vector<emInt>> &tetConns, const std::vector<emInt> &header,
			const std::vector<std::pair<emInt,emInt>> &cellID2cellTypeLocal,
			const std::vector<std::vector<emInt>> &vTriConns,
			const std::vector<double> &vLengthSclae, 
			const std::vector<std::vector<emInt>> &vQuadConn,
			const std::vector<std::vector<emInt>> &vPyrConn, 
			const std::vector<std::vector<emInt>> &vPrismConn,
			const std::vector<std::vector<emInt>> &vHexConn);		
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
	emInt numBdryTrisFromReader()  const
	{
		return m_nTrisFromReader; 
	} 
	emInt numBdryQuadsFromReader() const
	{
		return m_nQuadsFromReader; 
	}
	emInt numCells() const {
		return numTets() + numPyramids() + numPrisms() 
		+ numHexes() 
		+ numBdryTris() + numBdryQuads();
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

	const double (*getAllCoords() const)[3]
	{
		return m_coords; 
	}
	const emInt  (*getAllTetConn() const)[4]
	{
		return m_TetConn;
	}
	const emInt (*getAllTriConn() const) [3]
	{
		return m_TriConn;
	}
	const emInt (*getAllQuadConn() const) [4]
	{
		return m_QuadConn; 
	}
	const emInt (*getAllPyrConn() const ) [5]
	{
		return m_PyrConn; 
	}
	const emInt (*getAllPrismConn() const) [6]
	{
		return m_PrismConn; 
	}
	const emInt (*getAllHexConn() const) [8]
	{
		return m_HexConn; 
	}
	const emInt* getHeader() const
	{
		return m_header; 
	}
	Mapping::MappingType getDefaultMappingType() const {
		return Mapping::Uniform;
	}

	virtual std::unique_ptr<UMesh> createFineUMesh(const emInt numDivs, Part& P,
			std::vector<CellPartData>& vecCPD, struct RefineStats& RS) const;


	virtual std::unique_ptr<ExaMesh> extractCoarseMesh(Part& P,	std::vector<CellPartData>& vecCPD, 
	const int numDivs,
			const std::unordered_set<TriFaceVerts> &tris= std::unordered_set<TriFaceVerts>(), 
			const std::unordered_set<QuadFaceVerts> &quads= std::unordered_set<QuadFaceVerts>(), 
			const emInt partID=-1) const;

	void setupCellDataForPartitioning(std::vector<CellPartData>& vecCPD,
			double &xmin, double& ymin, double& zmin, double& xmax, double& ymax,
			double& zmax) const;

	bool writeVTKFile(const char fileName[]);
	bool writeUGridFile(const char fileName[]);

	size_t getFileImageSize() const {
		return m_fileImageSize;
	}

	void incrementVertIndices(emInt* conn, emInt size, int inc);
	void calcMemoryRequirements (const UMesh &UMIn, const int nDivs); 
	void buildCell2CellConn(multimpFace2Cell& face2cell, const emInt nCells);
	void buidCell2FacesConn(std::pair<emInt, emInt> cellInfo, emInt v0 , emInt v1, emInt v2); 
	void buidCell2FacesConn(std::pair<emInt, emInt> cellInfo, emInt v0 , emInt v1, emInt v2, emInt v3);
	void testCell2CellConn(emInt nCells); 
	void testCell2FaceConn(emInt nCells);
	std::vector<std::pair<emInt,emInt>> getCellID2CellType () const 
	{
		return cellID2cellTypeLocalID; 
	}
	void setCellId2CellTypeLocal (const std::vector<std::pair<emInt,emInt>>& cellIDs )
	{
		cellID2cellTypeLocalID=cellIDs; 
	}

	std::unique_ptr<UMesh>  
	Extract(const emInt partID, const std::vector<emInt> &partcells , const int numDivs, 
	const std::unordered_set<TriFaceVerts> tris= std::unordered_set<TriFaceVerts>(), 
	const std::unordered_set<QuadFaceVerts> quads= std::unordered_set<QuadFaceVerts>()) const;

	void convertToUmeshFormat(); 
	void 
	partFaceMatching(const std::vector<std::vector<emInt>> &part2cells,	
		 std::vector<std::unordered_set<TriFaceVerts>>  &tris,
		 std::vector<std::unordered_set<QuadFaceVerts>> &quads, size_t &totalTriSize, size_t &totalQuadSize) const;




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
