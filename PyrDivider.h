/*
 * PyrDivider.h
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#ifndef APPS_EXAMESH_PYRDIVIDER_H_
#define APPS_EXAMESH_PYRDIVIDER_H_

#include "examesh.h"
#include "CellDivider.h"

class PyrDivider: public CellDivider {
public:
	PyrDivider(UMesh *pVolMesh, const int segmentsPerEdge)
      :
			CellDivider(pVolMesh, segmentsPerEdge) {
    vertIJK[0][0] = 0;
    vertIJK[0][1] = 0;
    vertIJK[0][2] = nDivs;

    vertIJK[1][0] = nDivs;
    vertIJK[1][1] = 0;
    vertIJK[1][2] = nDivs;

    vertIJK[2][0] = nDivs;
    vertIJK[2][1] = nDivs;
    vertIJK[2][2] = nDivs;

    vertIJK[3][0] = 0;
    vertIJK[3][1] = nDivs;
    vertIJK[3][2] = nDivs;

    vertIJK[4][0] = 0;
    vertIJK[4][1] = 0;
    vertIJK[4][2] = 0;

		numVerts = 5;
		numEdges = 8;
		numTriFaces = 4;
		numQuadFaces = 1;
		edgeVertIndices[0][0] = 0;
		edgeVertIndices[0][1] = 1;

		edgeVertIndices[1][0] = 0;
		edgeVertIndices[1][1] = 3;

		edgeVertIndices[2][0] = 0;
		edgeVertIndices[2][1] = 4;

		edgeVertIndices[3][0] = 1;
		edgeVertIndices[3][1] = 2;

		edgeVertIndices[4][0] = 1;
		edgeVertIndices[4][1] = 4;

		edgeVertIndices[5][0] = 2;
		edgeVertIndices[5][1] = 3;

		edgeVertIndices[6][0] = 2;
		edgeVertIndices[6][1] = 4;

		edgeVertIndices[7][0] = 3;
		edgeVertIndices[7][1] = 4;

		faceVertIndices[0][0] = 0;
		faceVertIndices[0][1] = 1;
		faceVertIndices[0][2] = 2;
		faceVertIndices[0][3] = 3;

		faceVertIndices[1][0] = 0;
		faceVertIndices[1][1] = 4;
		faceVertIndices[1][2] = 1;

		faceVertIndices[2][0] = 1;
		faceVertIndices[2][1] = 4;
		faceVertIndices[2][2] = 2;

		faceVertIndices[3][0] = 2;
		faceVertIndices[3][1] = 4;
		faceVertIndices[3][2] = 3;

		faceVertIndices[4][0] = 3;
		faceVertIndices[4][1] = 4;
		faceVertIndices[4][2] = 0;
  }
	~PyrDivider() {
	}
	void divideInterior(const emInt verts[]);
  void createNewCells();
};

#endif /* APPS_EXAMESH_PYRDIVIDER_H_ */
