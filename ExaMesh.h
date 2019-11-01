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

class UMesh;

struct MeshSize {
	emInt nBdryVerts, nVerts, nBdryTris, nBdryQuads, nTets, nPyrs, nPrisms,
			nHexes;
};

class ExaMesh {
protected:
	double *m_lenScale;
	emInt *m_coarseGlobalIndices;

	void setupLengthScales();

public:
	ExaMesh() :
			m_lenScale(nullptr), m_coarseGlobalIndices(nullptr) {
	}
	virtual ~ExaMesh() {
		if (m_lenScale) delete[] m_lenScale;
		if (m_coarseGlobalIndices) delete[] m_coarseGlobalIndices;
	}
	virtual double getX(const emInt vert) const = 0;
	virtual double getY(const emInt vert) const = 0;
	virtual double getZ(const emInt vert) const =0;
	virtual void getCoords(const emInt vert, double coords[3]) const = 0;

	virtual emInt numVerts() const = 0;
	virtual emInt numBdryVerts() const = 0;
	virtual emInt numBdryTris() const = 0;
	virtual emInt numBdryQuads() const = 0;
	virtual emInt numTets() const = 0;
	virtual emInt numPyramids() const = 0;
	virtual emInt numPrisms() const = 0;
	virtual emInt numHexes() const = 0;
	virtual emInt numVertsToCopy() const {
		return numVerts();
	}

	virtual emInt addVert(const double newCoords[3],
			const emInt coarseGlobalIndex =
			EMINT_MAX) = 0;
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

	virtual void refineForParallel(const emInt numDivs,
			const emInt maxCellsPerPart) const;

	virtual std::unique_ptr<UMesh> createFineUMesh(const emInt numDivs, Part& P,
			std::vector<CellPartData>& vecCPD) const = 0;

	virtual void setupCellDataForPartitioning(std::vector<CellPartData>& vecCPD,
			double &xmin, double& ymin, double& zmin, double& xmax, double& ymax,
			double& zmax) const = 0;

protected:
	void addCellToPartitionData(const emInt* verts, emInt nPts, emInt ii,
			int type, std::vector<CellPartData>& vecCPD, double& xmin, double& ymin,
			double& zmin, double& xmax, double& ymax, double& zmax) const;
private:
	void findCentroidOfVerts(const emInt* verts, emInt nPts, double& x, double& y,
			double& z) const;
};

template<typename T>
void addUniquely(exaSet<T>& mySet, T& val) {
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
		const int nDivs);

bool partitionCells(const ExaMesh* const pEM, const emInt nPartsToMake,
		std::vector<Part>& parts, std::vector<CellPartData>& vecCPD);

void sortVerts3(const emInt input[3], emInt output[3]);
void sortVerts4(const emInt input[4], emInt output[4]);

#endif /* SRC_EXAMESH_H_ */
