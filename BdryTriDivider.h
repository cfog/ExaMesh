/*
 * BdryTriDivider.h
 *
 *  Created on: Jul. 8, 2019
 *      Author: cfog
 */

#ifndef SRC_BDRYTRIDIVIDER_H_
#define SRC_BDRYTRIDIVIDER_H_

#include "CellDivider.h"

class BdryTriDivider: public CellDivider {
public:
	BdryTriDivider(UMesh *pVolMesh, const int segmentsPerEdge) :
			CellDivider(pVolMesh, segmentsPerEdge) {
		vertIJK[0][0] = 0;
		vertIJK[0][1] = 0;
		vertIJK[0][2] = 0;

		vertIJK[1][0] = nDivs;
		vertIJK[1][1] = 0;
		vertIJK[1][2] = 0;

		vertIJK[2][0] = 0;
		vertIJK[2][1] = nDivs;
		vertIJK[2][2] = 0;

		numVerts = 3;
		numEdges = 3;
		numTriFaces = 1;
		numQuadFaces = 0;
		edgeVertIndices[0][0] = 0;
		edgeVertIndices[0][1] = 1;

		edgeVertIndices[1][0] = 1;
		edgeVertIndices[1][1] = 2;

		edgeVertIndices[2][0] = 2;
		edgeVertIndices[2][1] = 0;

		faceVertIndices[0][0] = 0;
		faceVertIndices[0][1] = 1;
		faceVertIndices[0][2] = 2;
	}
	~BdryTriDivider() {
	}
	void divideInterior();
	void createNewCells();
	void setupCoordMapping(const emInt /*verts*/[]) {
	}
	void getPhysCoordsFromParamCoords(const double /*uvw*/[], double /*xyz*/[]) {
	}
};


#endif /* SRC_BDRYTRIDIVIDER_H_ */
