/*
 * HexDivider.h
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#ifndef APPS_EXAMESH_HEXDIVIDER_H_
#define APPS_EXAMESH_HEXDIVIDER_H_

#include "examesh.h"
#include "CellDivider.h"

class HexDivider: public CellDivider {
public:
	HexDivider(UMesh *pVolMesh, const int segmentsPerEdge)
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

    vertIJK[5][0] = nDivs;
    vertIJK[5][1] = 0;
    vertIJK[5][2] = 0;

    vertIJK[6][0] = nDivs;
    vertIJK[6][1] = nDivs;
    vertIJK[6][2] = 0;

    vertIJK[7][0] = 0;
    vertIJK[7][1] = nDivs;
    vertIJK[7][2] = 0;

		numVerts = 8;
		numEdges = 12;
		numTriFaces = 0;
		numQuadFaces = 6;
		edgeVertIndices[0][0] = 0;
		edgeVertIndices[0][1] = 1;

		edgeVertIndices[1][0] = 1;
		edgeVertIndices[1][1] = 2;

		edgeVertIndices[2][0] = 2;
		edgeVertIndices[2][1] = 3;

		edgeVertIndices[3][0] = 3;
		edgeVertIndices[3][1] = 0;

		edgeVertIndices[4][0] = 4;
		edgeVertIndices[4][1] = 5;

		edgeVertIndices[5][0] = 5;
		edgeVertIndices[5][1] = 6;

		edgeVertIndices[6][0] = 6;
		edgeVertIndices[6][1] = 7;

		edgeVertIndices[7][0] = 7;
		edgeVertIndices[7][1] = 4;

		edgeVertIndices[8][0] = 0;
		edgeVertIndices[8][1] = 4;

		edgeVertIndices[9][0] = 1;
		edgeVertIndices[9][1] = 5;

		edgeVertIndices[10][0] = 2;
		edgeVertIndices[10][1] = 6;

		edgeVertIndices[11][0] = 3;
		edgeVertIndices[11][1] = 7;

		faceVertIndices[0][0] = 0;
		faceVertIndices[0][1] = 1;
		faceVertIndices[0][2] = 2;
		faceVertIndices[0][3] = 3;

		faceVertIndices[1][0] = 7;
		faceVertIndices[1][1] = 6;
		faceVertIndices[1][2] = 5;
		faceVertIndices[1][3] = 4;

		faceVertIndices[2][0] = 0;
		faceVertIndices[2][1] = 4;
		faceVertIndices[2][2] = 5;
		faceVertIndices[2][3] = 1;

		faceVertIndices[3][0] = 1;
		faceVertIndices[3][1] = 5;
		faceVertIndices[3][2] = 6;
		faceVertIndices[3][3] = 2;

		faceVertIndices[4][0] = 2;
		faceVertIndices[4][1] = 6;
		faceVertIndices[4][2] = 7;
		faceVertIndices[4][3] = 3;

		faceVertIndices[5][0] = 3;
		faceVertIndices[5][1] = 7;
		faceVertIndices[5][2] = 4;
		faceVertIndices[5][3] = 0;
  }
	~HexDivider() {
	}
	void divideInterior(const emInt verts[]);
  void createNewCells();
};

#endif /* APPS_EXAMESH_HEXDIVIDER_H_ */
