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


#include "ExaMesh.h"
#include "Mapping.h"
#include "UMesh.h"

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
			const double dihedral, EdgeVerts &EV);

	void getQuadVerts(exaSet<QuadFaceVerts> &vertsOnQuads, const int face,
			QuadFaceVerts &QFV);

	typename exaSet<TriFaceVerts>::iterator getTriVerts(
			exaSet<TriFaceVerts> &vertsOnTris,
			const int face,
			TriFaceVerts& TFV,
			bool& shouldErase);
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
