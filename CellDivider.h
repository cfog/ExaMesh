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
	void getEdgeVerts(exa_map<Edge, EdgeVerts> &vertsOnEdges, const int edge,
			const double dihedral, EdgeVerts &EV);

	void getQuadVerts(exa_set<QuadFaceVerts> &vertsOnQuads, const int face,
			QuadFaceVerts &QFV);

	typename exa_set<TriFaceVerts>::iterator getTriVerts(
			exa_set<TriFaceVerts> &vertsOnTris,
			const int face,
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
	void divideEdges(exa_map<Edge, EdgeVerts> &vertsOnEdges);
	void divideFaces(exa_set<TriFaceVerts> &vertsOnTris,
	exa_set<QuadFaceVerts> &vertsOnQuads);
	virtual void divideInterior() = 0;
	virtual void createNewCells() = 0;
	virtual void setupCoordMapping(const emInt verts[]) = 0;
	virtual void getPhysCoordsFromParamCoords(const double uvw[],
			double xyz[]) = 0;
};

#endif /* SRC_CELLDIVIDER_H_ */
