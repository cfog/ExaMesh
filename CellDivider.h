/*
 * CellDivider.h
 *
 *  Created on: Jul. 3, 2019
 *      Author: cfog
 */

#ifndef SRC_CELLDIVIDER_H_
#define SRC_CELLDIVIDER_H_

#include <algorithm>
#include <cmath>

#ifdef USE_ORDERED

#include <set>
#include <map>
#define exaSet std::set
#define exaMap std::map

#else

#include <unordered_set>
#include <unordered_map>

#define exaSet std::unordered_set
#define exaMap std::unordered_map

#endif

#include "examesh.h"
#include "Mapping.h"
#include "UMesh.h"

class Edge {
private:
	emInt v0, v1;
public:
	Edge(emInt vA, emInt vB) {
		if (vA < vB) {
			v0 = vA;
			v1 = vB;
		}
		else {
			v0 = vB;
			v1 = vA;
		}
	}
	bool operator<(const Edge &E) const {
		return (v0 < E.v0 || (v0 == E.v0 && v1 < E.v1));
	}
	bool operator==(const Edge &E) const {
		return v0 == E.v0 && v1 == E.v1;
	}
	emInt getV0() const {
		return v0;
	}

	emInt getV1() const {
		return v1;
	}
};

struct EdgeVerts {
	emInt verts[MAX_DIVS + 1];
};

struct TriFaceVerts {
	emInt corners[3], sorted[3];
	emInt intVerts[MAX_DIVS - 2][MAX_DIVS - 2];
	TriFaceVerts() {
	}
	TriFaceVerts(const emInt v0, const emInt v1, const emInt v2);
	void setupSorted();
};

struct QuadFaceVerts {
	emInt corners[4], sorted[4];
	emInt intVerts[MAX_DIVS - 1][MAX_DIVS - 1];
	QuadFaceVerts() {
	}
	QuadFaceVerts(const emInt v0, const emInt v1, const emInt v2, const emInt v3);
	void setupSorted();
};

namespace std {
	template<> struct hash<TriFaceVerts> {
		typedef TriFaceVerts argument_type;
		typedef std::size_t result_type;
		result_type operator()(const argument_type& TFV) const noexcept
		{
			const result_type h0 = TFV.sorted[0];
			const result_type h1 = TFV.sorted[1];
			const result_type h2 = TFV.sorted[2];
			return (h0 ^ (h1 << 1)) ^ (h2 << 2);
		}
	};

	template<> struct hash<QuadFaceVerts> {
		typedef QuadFaceVerts argument_type;
		typedef std::size_t result_type;
		result_type operator()(const argument_type& QFV) const noexcept
		{
			const result_type h0 = QFV.sorted[0];
			const result_type h1 = QFV.sorted[1];
			const result_type h2 = QFV.sorted[2];
			const result_type h3 = QFV.sorted[3];
			return h0 ^ (h1 << 1) ^ (h2 << 2) ^ (h3 << 3);
		}
	};

	template<> struct hash<Edge> {
		typedef Edge argument_type;
		typedef std::size_t result_type;
		result_type operator()(const argument_type& E) const noexcept
		{
			const result_type h0 = E.getV0();
			const result_type h1 = E.getV1();
			return (h0 ^ (h1 << 1));
		}
	};
}


bool operator==(const TriFaceVerts& a, const TriFaceVerts& b);
bool operator<(const TriFaceVerts& a, const TriFaceVerts& b);
bool operator==(const QuadFaceVerts& a, const QuadFaceVerts& b);
bool operator<(const QuadFaceVerts& a, const QuadFaceVerts& b);

class CellDivider {
protected:
	UMesh *m_pMesh;
	Mapping *m_Map;
	emInt (*localVerts)[MAX_DIVS + 1][MAX_DIVS + 1];
	int edgeVertIndices[12][2];
	int faceVertIndices[6][4];
	int numTriFaces, numQuadFaces, numEdges, numVerts;
	int vertIJK[8][3];
	double uvwIJK[8][3];
	emInt cellVerts[8];
	int nDivs;

	// Used by both tets and pyramids.
	int checkOrient3D(const emInt verts[4]) const;
private:
	void getEdgeVerts(exaMap<Edge, EdgeVerts> &vertsOnEdges, const int edge,
			EdgeVerts &EV);

	void getQuadVerts(exaSet<QuadFaceVerts> &vertsOnQuads, const int face,
			QuadFaceVerts &QFV);

	void getTriVerts(exaSet<TriFaceVerts> &vertsOnTris,
			const int face,
			TriFaceVerts& TFV);
public:
	CellDivider(UMesh *pVolMesh, const emInt segmentsPerEdge) :
			m_pMesh(pVolMesh), m_Map(nullptr),
					numTriFaces(0), numQuadFaces(0), numEdges(0),
					numVerts(0), nDivs(segmentsPerEdge) {
		localVerts = new emInt[MAX_DIVS + 1][MAX_DIVS + 1][MAX_DIVS + 1];
//		for (int ii = 0; ii <= MAX_DIVS; ii++) {
//			for (int jj = 0; jj <= MAX_DIVS; jj++) {
//				for (int kk = 0; kk <= MAX_DIVS; kk++) {
//					localVerts[ii][jj][kk] = 100;
//				}
//			}
//		}
	}
	virtual ~CellDivider() {
		delete[] localVerts;
		if (m_Map) delete m_Map;
	}
	void divideEdges(exaMap<Edge, EdgeVerts> &vertsOnEdges);
	void divideFaces(exaSet<TriFaceVerts> &vertsOnTris,
	exaSet<QuadFaceVerts> &vertsOnQuads);
	virtual void divideInterior() = 0;
	virtual void createNewCells() = 0;
	virtual void setupCoordMapping(const emInt verts[]) = 0;
	virtual void getPhysCoordsFromParamCoords(const double uvw[],
			double xyz[]) = 0;
};

#endif /* SRC_CELLDIVIDER_H_ */
