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
 * examesh.h
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#ifndef SRC_EXAMESH_H_
#define SRC_EXAMESH_H_

#include <limits.h>
#include <assert.h>
#include <memory>

#include "Mapping.h"
#include "Part.h"
#include "exa-defs.h"
#include "FaceVerts.h"
#include "mpiDefs.h"
#include "ParallelTester.h"
#include <set>

class UMesh;
class CubicMesh; 

using vecSharePtrUmesh     = std::vector<std::shared_ptr<UMesh>>; 
using vecSharePtrCubicMesh = std::vector<std::shared_ptr<CubicMesh>>;
struct MeshSize {
	ssize_t nBdryVerts, nVerts, nBdryTris, nBdryQuads, nTets, nPyrs, nPrisms,
			nHexes;
};

class ExaMesh {
protected:
	double *m_lenScale;
	exa_set<QuadFaceVerts> m_partQuads; 
	exa_set<TriFaceVerts>  m_partTris; 
	exa_set<TriFaceVerts>  m_refinedPartTris; 
	exa_set<QuadFaceVerts> m_refinedPartQuads; 

	std::vector<std::vector<emInt>>   vcell2cell;
	std::vector<emInt>                vcellID2type;
	std::vector<std::pair<emInt,emInt>>  cellID2cellTypeLocalID;
	std::map < std::pair<emInt,emInt>, std::set<std::set<emInt>>>   cell2faces;


	void setupLengthScales();

public:
	ExaMesh() :
			m_lenScale(nullptr) {
	}
	virtual ~ExaMesh() {
		if (m_lenScale) delete[] m_lenScale;
	}

	static std::unique_ptr<ExaMesh> readMeshFromFile(const std::string fileName,
			const std::string fileSuffix, const std::string fileInfix = "");

	virtual std::unique_ptr<UMesh>
	subdivideMesh(const emInt nDivs, const emInt partID = 0) const = 0;

	virtual double getX(const emInt vert) const = 0;
	virtual double getY(const emInt vert) const = 0;
	virtual double getZ(const emInt vert) const =0;
	virtual void getCoords(const emInt vert, double coords[3]) const = 0;

	virtual emInt numVerts() const = 0;
	virtual emInt numVertsAllocated() const = 0;
	virtual emInt numBdryVerts() const = 0;
	virtual emInt numBdryTris() const = 0;
	virtual emInt numBdryQuads() const = 0;
	virtual emInt numTets() const = 0;
	virtual emInt numPyramids() const = 0;
	virtual emInt numPrisms() const = 0;
	virtual emInt numHexes() const = 0;
	virtual emInt numCells() const = 0 ; 
	emInt numVolCells() const {
	    return numCells() - numBdryTris() - numBdryQuads();
	}
		
	virtual emInt numVertsToCopy() const {
		return numVerts();
	}

	virtual emInt addVert(const double newCoords[3]) = 0;
	virtual emInt addBdryTri(const emInt verts[]) = 0;
	virtual emInt addBdryQuad(const emInt verts[]) = 0;
	virtual emInt addTet(const emInt verts[]) = 0;
	virtual emInt addPyramid(const emInt verts[]) = 0;
	virtual emInt addPrism(const emInt verts[]) = 0;
	virtual emInt addHex(const emInt verts[]) = 0;


	const virtual emInt* getBdryTriConn(const emInt bdryTri) const=0;
	const virtual emInt* getBdryQuadConn(const emInt bdryQuad) const=0;
	const virtual emInt* getTetConn(const emInt tet) const=0;
	const virtual emInt* getPyrConn(const emInt pyr) const=0;
	const virtual emInt* getPrismConn(const emInt prism) const=0;
	const virtual emInt* getHexConn(const emInt hex) const=0;

	virtual Mapping::MappingType getDefaultMappingType() const = 0;

	std::size_t getCellConnSize (const emInt cellID)
	const
	{
		return vcell2cell[cellID].size();
	}
	emInt getCellConn (const emInt cellID, const emInt neighID)
	const
	{
		return vcell2cell[cellID][neighID];
	}
	emInt getCellType (const emInt cellID)
	const
	{
		return vcellID2type[cellID];
	}
	emInt FastpartFaceMatching(const emInt nParts, const std::vector<std::vector<emInt>> &part2cells,
			const std::vector<emInt> &cell2part,
			vecVecTri &tris, vecVecQuad &quads) const;
	std::vector<std::pair<emInt,emInt>> getCellID2CellType2LocalID() const
	{
		return cellID2cellTypeLocalID;
	}
	void getFaceLists (const emInt ind, const emInt type,
			const emInt partID, const emInt numDivs,
			std::vector<TriFaceVerts> &tris,
			std::vector<QuadFaceVerts> &quads) const;

	void printMeshSizeStats();
	double getLengthScale(const emInt vert) const {
		assert(vert < numVertsAllocated());
		// TODO Would be better to always associate a length scale with the mesh.
		if (m_lenScale) {
			return m_lenScale[vert];
		}
		else {
			return 1;
		}
	}
	const double* getAllLengthScale() const
	{
		return m_lenScale; 
	}
	void setLengthScale(const emInt vert, const double len) const {
		assert(vert < numVertsAllocated());
		assert(len > 0);
		// TODO Would be better to always associate a length scale with the mesh.
		if (m_lenScale) {
			m_lenScale[vert] = len;
		}
	}
	MeshSize computeFineMeshSize(const int nDivs) const;

	void buildCellToCellConnectivity();
	void buildFaceCellConnectivity();
	void buildCell2CellConnFromFace2Cell(multimpFace2Cell& face2cell, const emInt nCells);
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
	virtual emInt getTriType() const = 0;
	virtual emInt getQuadType() const = 0;
	virtual emInt getTetType() const = 0;
	virtual emInt getPyrType() const = 0;
	virtual emInt getPrismType() const = 0;
	virtual emInt getHexType() const = 0;

	virtual void refineForParallel(const emInt numDivs,
			const emInt maxCellsPerPart) const;
	virtual std::unique_ptr<UMesh> createFineUMesh(const emInt numDivs, Part& P,
			std::vector<CellPartData>& vecCPD, struct RefineStats& RS) const = 0;

	void refineForMPI( const int numDivs ,ParallelTester* tester,
	const char MeshType, std::string mshName,FILE* fileAllTimes) const;		

	virtual std::unique_ptr<ExaMesh> extractCoarseMeshPseudoParallel(Part& P,	std::vector<CellPartData>& vecCPD, 
	const int numDivs,
			const std::unordered_set<TriFaceVerts> &tris= std::unordered_set<TriFaceVerts>(), 
			const std::unordered_set<QuadFaceVerts> &quads= std::unordered_set<QuadFaceVerts>(), 
			const emInt partID=-1) const=0;		

	virtual std::unique_ptr<ExaMesh>
	extractCoarseMeshMPI(const emInt partID, const std::vector<emInt> &partcells , const int numDivs,
			const std::unordered_set<TriFaceVerts>& tris,
			const std::unordered_set<QuadFaceVerts>& quads) const = 0;

	void TestMPI(const emInt &nDivs, const emInt &nParts, ParallelTester* tester, 
	const char MeshType); 
	virtual void setupCellDataForPartitioning(std::vector<CellPartData>& vecCPD,
			double &xmin, double& ymin, double& zmin, double& xmax, double& ymax,
			double& zmax) const = 0;
	void prettyPrintCellCount(size_t cells, const char* prefix) const;

	void addPartTritoSet(const TriFaceVerts &obj){
		auto iter = m_partTris.find(obj);
	
		if (iter != m_partTris.end()) {
			m_partTris.erase(iter);
		}
		else {
			m_partTris.insert(obj);
		}

	}
	void addPartQuadtoSet(const QuadFaceVerts &obj){
		auto iter = m_partQuads.find(obj);
	
		if (iter != m_partQuads.end()) {
			m_partQuads.erase(iter);
		}
		else {
			m_partQuads.insert(obj);
		}

	}
	void addRefinedPartTritoSet(const TriFaceVerts &obj){
		auto iter = m_refinedPartTris.find(obj);
	
		if (iter != m_refinedPartTris.end()) {
			m_refinedPartTris.erase(iter);
		}
		else {
			m_refinedPartTris.insert(obj);
		}

	}
	void addRefinedPartQuadtoSet(const QuadFaceVerts &obj){
		auto iter = m_refinedPartQuads.find(obj);
	
		if (iter != m_refinedPartQuads.end()) {
			m_refinedPartQuads.erase(iter);
		}
		else {
			m_refinedPartQuads.insert(obj);
		}

	}

	size_t getSizePartTris()const
	{
		return m_partTris.size();
	}
	size_t getSizePartQuads()const
	{
		return m_partQuads.size();
	}
	const exa_set<QuadFaceVerts>& getTempQuadPart    () const
	{
		return m_partQuads; 
	}
	const exa_set<TriFaceVerts>&  getTempTriPart     () const 
	{
		return m_partTris; 
	}
	const exa_set<TriFaceVerts>&  getRefinedPartTris () const
	{
		return m_refinedPartTris; 
	}
	const exa_set<QuadFaceVerts>& getRefinedPartQuads() const 
	{
		return m_refinedPartQuads; 
	}
	//void refineMPI();
	//virtual void buildCell2CellConn(const std::multimap < std::set<emInt>, std::pair<emInt,emInt>> & face2cell, const emInt nCells)=0; 	 

protected:
	void addCellToPartitionData(const emInt* verts, emInt nPts, emInt ii,
			int type, std::vector<CellPartData>& vecCPD, double& xmin, double& ymin,
			double& zmin, double& xmax, double& ymax, double& zmax) const;
private:
	void findCentroidOfVerts(const emInt* verts, emInt nPts, double& x, double& y,
			double& z) const;
};

template<typename T>
void addUniquely(exa_set<T>& mySet, T& val) {

	typename exa_set<T>::iterator vertIter, VIend = mySet.end();

	auto inserResult = mySet.insert(val);

	if(!inserResult.second)
	{
		mySet.erase(inserResult.first);
	}

	// auto iter = mySet.find(val);
	// if (iter != mySet.end()) {
	// 	mySet.erase(iter);
	// }
	// else {
	// 	mySet.insert(val);
	// }
}
template<typename T>
void addUniquely(std::set<T> &mySet, T& val) {
	auto iter = mySet.find(val);
	
	if (iter != mySet.end()) {
		mySet.erase(iter);
	}
	else {
		mySet.insert(val);
	}
}

bool computeMeshSize(const struct MeshSize& MSIn, const emInt nDivs,
		struct MeshSize& MSOut);

// Defined elsewhere.
emInt subdividePartMesh(const ExaMesh * const pVM_input,
		UMesh * const pVM_output,
		const int nDivs, const emInt partID=-1);

bool partitionCells(const ExaMesh* const pEM, const emInt nPartsToMake,
		std::vector<Part>& parts, std::vector<CellPartData>& vecCPD);

void sortVerts3(const emInt input[3], emInt output[3]);
void sortVerts4(const emInt input[4], emInt output[4]);

#endif /* SRC_EXAMESH_H_ */
