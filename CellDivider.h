/*
 * CellDivider.h
 *
 *  Created on: Jul. 3, 2019
 *      Author: cfog
 */

#ifndef SRC_CELLDIVIDER_H_
#define SRC_CELLDIVIDER_H_

#include <set>

#include "examesh.h"
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
	emInt corners[3];
	emInt intVerts[MAX_DIVS - 2][MAX_DIVS - 2];
	TriFaceVerts() {
	}
	TriFaceVerts(const emInt v0, const emInt v1, const emInt v2);
};

struct QuadFaceVerts {
	emInt corners[4];
	emInt intVerts[MAX_DIVS - 1][MAX_DIVS - 1];
	QuadFaceVerts() {
	}
	QuadFaceVerts(const emInt v0, const emInt v1, const emInt v2, const emInt v3);
};

bool operator<(const TriFaceVerts& a, const TriFaceVerts& b);
bool operator<(const QuadFaceVerts& a, const QuadFaceVerts& b);

class CellDivider {
protected:
	UMesh *m_pMesh;
	emInt (*localVerts)[MAX_DIVS + 1][MAX_DIVS + 1];
	int edgeVertIndices[12][2];
	int faceVertIndices[6][4];
	int numTriFaces, numQuadFaces, numEdges, numVerts;
	int vertIJK[8][3];
	double uvwIJK[8][3];
	emInt cellVerts[8];
	int nDivs;

private:
	void getEdgeVerts(std::map<Edge, EdgeVerts> &vertsOnEdges, const int edge,
			EdgeVerts &EV);

	void getQuadVerts(std::set<QuadFaceVerts> &vertsOnQuads, const int face,
			QuadFaceVerts &QFV);

	void getTriVerts(std::set<TriFaceVerts> &vertsOnTris, const int face,
			TriFaceVerts& TFV);

public:
	CellDivider(UMesh *pVolMesh, const emInt segmentsPerEdge) :
			m_pMesh(pVolMesh), numTriFaces(0), numQuadFaces(0), numEdges(0),
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
	}
	void divideEdges(std::map<Edge, EdgeVerts> &vertsOnEdges);
	void divideFaces(std::set<TriFaceVerts> &vertsOnTris,
			std::set<QuadFaceVerts> &vertsOnQuads);
	virtual void divideInterior() = 0;
	virtual void createNewCells() = 0;
	virtual void setupCoordMapping(const emInt verts[]) = 0;
	virtual void getPhysCoordsFromParamCoords(const double uvw[],
			double xyz[]) = 0;
};

#endif /* SRC_CELLDIVIDER_H_ */
