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


	void setupLengthScales();

public:
	ExaMesh() :
			m_lenScale(nullptr) {
	}
	virtual ~ExaMesh() {
		if (m_lenScale) delete[] m_lenScale;
	}
	virtual double getX(const emInt vert) const = 0;
	virtual double getY(const emInt vert) const = 0;
	virtual double getZ(const emInt vert) const =0;
	virtual void getCoords(const emInt vert, double coords[3]) const = 0;

	virtual emInt numVerts() const = 0;
	virtual emInt numBdryVerts() const = 0;
	virtual emInt numBdryTris() const = 0;
	virtual emInt numBdryQuads() const = 0;
	virtual emInt numBdryTrisFromReader()  const= 0 ; 
	virtual emInt numBdryQuadsFromReader() const= 0 ;
	virtual emInt numTets() const = 0;
	virtual emInt numPyramids() const = 0;
	virtual emInt numPrisms() const = 0;
	virtual emInt numHexes() const = 0;
	virtual emInt numCells() const = 0 ; 
		
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



	void printMeshSizeStats();
	double getLengthScale(const emInt vert) const {
		assert(vert < numVerts());
		// TODO Would be better to always associate a length scale with the mesh.
		if (m_lenScale) {
			return m_lenScale[vert];
		}
		else {
			return 1;
		}
	}
	void setLengthScale(const emInt vert, const double len) const {
		assert(vert < numVerts());
		assert(len > 0);
		// TODO Would be better to always associate a length scale with the mesh.
		if (m_lenScale) {
			m_lenScale[vert] = len;
		}
	}
	MeshSize computeFineMeshSize(const int nDivs) const;

	void buildFaceCellConnectivity();

	virtual void refineForParallel(const emInt numDivs,
			const emInt maxCellsPerPart) const;
	virtual std::unique_ptr<UMesh> createFineUMesh(const emInt numDivs, Part& P,
			std::vector<CellPartData>& vecCPD, struct RefineStats& RS) const = 0;

	void refineForMPI( const int numDivs ,ParallelTester* tester,
	const char MeshType, std::string mshName,FILE* fileAllTimes) const;		

	virtual std::unique_ptr<ExaMesh> extractCoarseMesh(Part& P,	std::vector<CellPartData>& vecCPD, 
	const int numDivs,
			const std::unordered_set<TriFaceVerts> &tris= std::unordered_set<TriFaceVerts>(), 
			const std::unordered_set<QuadFaceVerts> &quads= std::unordered_set<QuadFaceVerts>(), 
			const emInt partID=-1) const=0;		
			
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
	virtual void partFaceMatching(
		 std::vector<Part>& parts, const std::vector<CellPartData>& vecCPD,	
		 std::vector<std::unordered_set<TriFaceVerts>>  &tris,
		 std::vector<std::unordered_set<QuadFaceVerts>> &quads, size_t &totalTriSize, size_t &totalQuadSize )const=0;	
	//void refineMPI();
	//virtual void buildCell2CellConn(const std::multimap < std::set<emInt>, std::pair<emInt,emInt>> & face2cell, const emInt nCells)=0; 	 
	virtual std::size_t getCellConnSize (const emInt cellID) const = 0; 
	virtual emInt getCellConn (const emInt cellID, const emInt neighID) const = 0; 

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
	auto iter = mySet.find(val);
	if (iter != mySet.end()) {
		mySet.erase(iter);
	}
	else {
		mySet.insert(val);
	}
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
