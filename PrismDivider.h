/*
 * PrismDivider.h
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#ifndef APPS_EXAMESH_PRISMDIVIDER_H_
#define APPS_EXAMESH_PRISMDIVIDER_H_

#include "examesh.h"
#include "CellDivider.h"

class PrismDivider: public CellDivider {
public:
	PrismDivider(UMesh *pVolMesh, const int segmentsPerEdge)
:
			CellDivider(pVolMesh, segmentsPerEdge) {
    vertIJK[0][0] = 0;
    vertIJK[0][1] = 0;
    vertIJK[0][2] = nDivs;

    vertIJK[1][0] = nDivs;
    vertIJK[1][1] = 0;
    vertIJK[1][2] = nDivs;

    vertIJK[2][0] = 0;
    vertIJK[2][1] = nDivs;
    vertIJK[2][2] = nDivs;

    vertIJK[3][0] = 0;
    vertIJK[3][1] = 0;
    vertIJK[3][2] = 0;

    vertIJK[4][0] = nDivs;
    vertIJK[4][1] = 0;
    vertIJK[4][2] = 0;

    vertIJK[5][0] = 0;
    vertIJK[5][1] = nDivs;
    vertIJK[5][2] = 0;

		numVerts = 6;
		numEdges = 9;
		numTriFaces = 2;
		numQuadFaces = 3;
		edgeVertIndices[0][0] = 0;
		edgeVertIndices[0][1] = 1;

		edgeVertIndices[1][0] = 1;
		edgeVertIndices[1][1] = 2;

		edgeVertIndices[2][0] = 0;
		edgeVertIndices[2][1] = 2;

		edgeVertIndices[3][0] = 3;
		edgeVertIndices[3][1] = 4;

		edgeVertIndices[4][0] = 4;
		edgeVertIndices[4][1] = 5;

		edgeVertIndices[5][0] = 3;
		edgeVertIndices[5][1] = 5;

		edgeVertIndices[6][0] = 0;
		edgeVertIndices[6][1] = 3;

		edgeVertIndices[7][0] = 1;
		edgeVertIndices[7][1] = 4;

		edgeVertIndices[8][0] = 2;
		edgeVertIndices[8][1] = 5;

		faceVertIndices[0][0] = 2;
		faceVertIndices[0][1] = 1;
		faceVertIndices[0][2] = 4;
		faceVertIndices[0][3] = 5;

		faceVertIndices[1][0] = 1;
		faceVertIndices[1][1] = 0;
		faceVertIndices[1][2] = 3;
		faceVertIndices[1][3] = 4;

		faceVertIndices[2][0] = 0;
		faceVertIndices[2][1] = 2;
		faceVertIndices[2][2] = 5;
		faceVertIndices[2][3] = 3;

		faceVertIndices[3][0] = 0;
		faceVertIndices[3][1] = 1;
		faceVertIndices[3][2] = 2;

		faceVertIndices[4][0] = 5;
		faceVertIndices[4][1] = 4;
		faceVertIndices[4][2] = 3;
  }
	~PrismDivider() {
	}
	virtual void divideInterior(const emInt verts[]);
	virtual void createNewCells();
};

#endif /* APPS_EXAMESH_PRISMDIVIDER_H_ */
