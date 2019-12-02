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
 * BdryQuadDivider.h
 *
 *  Created on: Jul. 8, 2019
 *      Author: cfog
 */

#ifndef SRC_BDRYQUADDIVIDER_H_
#define SRC_BDRYQUADDIVIDER_H_

#include "CellDivider.h"

class BdryQuadDivider: public CellDivider {
public:
	BdryQuadDivider(UMesh *pVolMesh, const int segmentsPerEdge) :
			CellDivider(pVolMesh, segmentsPerEdge) {
		vertIJK[0][0] = 0;
		vertIJK[0][1] = 0;
		vertIJK[0][2] = 0;

		vertIJK[1][0] = nDivs;
		vertIJK[1][1] = 0;
		vertIJK[1][2] = 0;

		vertIJK[2][0] = nDivs;
		vertIJK[2][1] = nDivs;
		vertIJK[2][2] = 0;

		vertIJK[3][0] = 0;
		vertIJK[3][1] = nDivs;
		vertIJK[3][2] = 0;

		numVerts = 4;
		numEdges = 4;
		numTriFaces = 0;
		numQuadFaces = 1;
		edgeVertIndices[0][0] = 0;
		edgeVertIndices[0][1] = 1;

		edgeVertIndices[1][0] = 1;
		edgeVertIndices[1][1] = 2;

		edgeVertIndices[2][0] = 2;
		edgeVertIndices[2][1] = 3;

		edgeVertIndices[3][0] = 3;
		edgeVertIndices[3][1] = 0;

		faceVertIndices[0][0] = 0;
		faceVertIndices[0][1] = 1;
		faceVertIndices[0][2] = 2;
		faceVertIndices[0][3] = 3;
	}
	~BdryQuadDivider() {
	}
	void divideInterior();
	void createNewCells();

	// TODO: Currently, there's no coord mapping set up for bdry faces
	void setupCoordMapping(const emInt verts[]) {
		for (int ii = 0; ii < 4; ii++) {
			cellVerts[ii] = verts[ii];
		}

	}
	void getPhysCoordsFromParamCoords(const double /*uvw*/[], double /*xyz*/[]) {
	}
};


#endif /* SRC_BDRYQUADDIVIDER_H_ */
