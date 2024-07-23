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
 * PrismDivider.h
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#ifndef APPS_EXAMESH_PRISMDIVIDER_H_
#define APPS_EXAMESH_PRISMDIVIDER_H_

#include "CellDivider.h"
#include "ExaMesh.h"

class PrismDivider: public CellDivider {
	double xyzOffsetBot[3], uVecBot[3], vVecBot[3];
	double xyzOffsetTop[3], uVecTop[3], vVecTop[3];
public:
	PrismDivider(UMesh *pVolMesh, const ExaMesh* const pInitMesh,
			const int segmentsPerEdge,
			const Mapping::MappingType type = Mapping::Invalid)
:
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

    vertIJK[3][0] = 0;
    vertIJK[3][1] = 0;
    vertIJK[3][2] = nDivs;

    vertIJK[4][0] = nDivs;
    vertIJK[4][1] = 0;
    vertIJK[4][2] = nDivs;

    vertIJK[5][0] = 0;
    vertIJK[5][1] = nDivs;
    vertIJK[5][2] = nDivs;

		uvwIJK[0][0] = 0;
		uvwIJK[0][1] = 0;
		uvwIJK[0][2] = 0;

		uvwIJK[1][0] = 1;
		uvwIJK[1][1] = 0;
		uvwIJK[1][2] = 0;

		uvwIJK[2][0] = 0;
		uvwIJK[2][1] = 1;
		uvwIJK[2][2] = 0;

		uvwIJK[3][0] = 0;
		uvwIJK[3][1] = 0;
		uvwIJK[3][2] = 1;

		uvwIJK[4][0] = 1;
		uvwIJK[4][1] = 0;
		uvwIJK[4][2] = 1;

		uvwIJK[5][0] = 0;
		uvwIJK[5][1] = 1;
		uvwIJK[5][2] = 1;

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

		faceEdgeIndices[0][0] = 1;
		faceEdgeIndices[0][1] = 7;
		faceEdgeIndices[0][2] = 4;
		faceEdgeIndices[0][3] = 8;

		faceEdgeIndices[1][0] = 0;
		faceEdgeIndices[1][1] = 6;
		faceEdgeIndices[1][2] = 3;
		faceEdgeIndices[1][3] = 7;

		faceEdgeIndices[2][0] = 2;
		faceEdgeIndices[2][1] = 8;
		faceEdgeIndices[2][2] = 5;
		faceEdgeIndices[2][3] = 6;

		faceEdgeIndices[3][0] = 0;
		faceEdgeIndices[3][1] = 1;
		faceEdgeIndices[3][2] = 2;

		faceEdgeIndices[4][0] = 4;
		faceEdgeIndices[4][1] = 3;
		faceEdgeIndices[4][2] = 5;

		Mapping::MappingType myType = type;
		if (myType == Mapping::Invalid) {
			myType = pInitMesh->getDefaultMappingType();
		}
		switch (myType) {
		case Mapping::Lagrange:
			m_Map = new LagrangeCubicPrismMapping(pInitMesh);
			break;
		case Mapping::Uniform:
		default:
			m_Map = new Q1PrismMapping(pInitMesh);
			break;
		}
	}
	~PrismDivider() {
	}
//	virtual void divideInterior();
	virtual void createNewCells();
	void setupCoordMapping(const emInt verts[]);
	void getPhysCoordsFromParamCoords(const double uvw[], double xyz[]);

	virtual emInt maxI(const emInt j, const emInt /*k*/) const {return nDivs - j;}
	virtual emInt maxJ(const emInt i, const emInt /*k*/) const {return nDivs - i;}

};

#endif /* APPS_EXAMESH_PRISMDIVIDER_H_ */
