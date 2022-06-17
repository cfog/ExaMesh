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
 * HexDivider.h
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#ifndef APPS_EXAMESH_HEXDIVIDER_H_
#define APPS_EXAMESH_HEXDIVIDER_H_

#include "CellDivider.h"
#include "ExaMesh.h"

class HexDivider: public CellDivider {
	double xyzOffsetBot[3], uVecBot[3], vVecBot[3], uvVecBot[3];
	double xyzOffsetTop[3], uVecTop[3], vVecTop[3], uvVecTop[3];

public:
	HexDivider(UMesh *pVolMesh, const ExaMesh* const pInitMesh,
			const int segmentsPerEdge,
			const Mapping::MappingType type = Mapping::Uniform)
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

		uvwIJK[0][0] = 0;
		uvwIJK[0][1] = 0;
		uvwIJK[0][2] = 0;

		uvwIJK[1][0] = 1;
		uvwIJK[1][1] = 0;
		uvwIJK[1][2] = 0;

		uvwIJK[2][0] = 1;
		uvwIJK[2][1] = 1;
		uvwIJK[2][2] = 0;

		uvwIJK[3][0] = 0;
		uvwIJK[3][1] = 1;
		uvwIJK[3][2] = 0;

		uvwIJK[4][0] = 0;
		uvwIJK[4][1] = 0;
		uvwIJK[4][2] = 1;

		uvwIJK[5][0] = 1;
		uvwIJK[5][1] = 0;
		uvwIJK[5][2] = 1;

		uvwIJK[6][0] = 1;
		uvwIJK[6][1] = 1;
		uvwIJK[6][2] = 1;

		uvwIJK[7][0] = 0;
		uvwIJK[7][1] = 1;
		uvwIJK[7][2] = 1;

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

		faceEdgeIndices[0][0] = 0;
		faceEdgeIndices[0][1] = 1;
		faceEdgeIndices[0][2] = 2;
		faceEdgeIndices[0][3] = 3;

		faceEdgeIndices[1][0] = 6;
		faceEdgeIndices[1][1] = 5;
		faceEdgeIndices[1][2] = 4;
		faceEdgeIndices[1][3] = 7;

		faceEdgeIndices[2][0] = 8;
		faceEdgeIndices[2][1] = 4;
		faceEdgeIndices[2][2] = 9;
		faceEdgeIndices[2][3] = 0;

		faceEdgeIndices[3][0] = 9;
		faceEdgeIndices[3][1] = 5;
		faceEdgeIndices[3][2] = 10;
		faceEdgeIndices[3][3] = 1;

		faceEdgeIndices[4][0] = 10;
		faceEdgeIndices[4][1] = 6;
		faceEdgeIndices[4][2] = 11;
		faceEdgeIndices[4][3] = 2;
	       
		faceEdgeIndices[5][0] = 11;
		faceEdgeIndices[5][1] = 7;
		faceEdgeIndices[5][2] = 8;
		faceEdgeIndices[5][3] = 3;
	      
		if (type == Mapping::LengthScale) {
			m_Map = new LengthScaleHexMapping(pInitMesh);
		}
		else if (type == Mapping::Lagrange) {
			m_Map = new LagrangeCubicHexMapping(pInitMesh);
		}
		else {
			m_Map = new Q1HexMapping(pInitMesh);
		}
  }
	~HexDivider() {
	}
	void divideInterior();
  void createNewCells();
	void setupCoordMapping(const emInt verts[]);
	void getPhysCoordsFromParamCoords(const double uvw[], double xyz[]);
};

#endif /* APPS_EXAMESH_HEXDIVIDER_H_ */
