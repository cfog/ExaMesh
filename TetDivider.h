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
 * TetDivider.h
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#ifndef APPS_EXAMESH_TETDIVIDER_H_
#define APPS_EXAMESH_TETDIVIDER_H_

#include "CellDivider.h"
#include "ExaMesh.h"
#include "Mapping.h"

class TetDivider: public CellDivider {
public:
	TetDivider(UMesh *pVolMesh, const ExaMesh* const pInitMesh,
			const int segmentsPerEdge,
			const Mapping::MappingType type =
					Mapping::Uniform)
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

		numVerts = 4;
		numEdges = 6;
		numTriFaces = 4;
		numQuadFaces = 0;
		edgeVertIndices[0][0] = 0;
		edgeVertIndices[0][1] = 1;

		edgeVertIndices[1][0] = 0;
		edgeVertIndices[1][1] = 2;

		edgeVertIndices[2][0] = 0;
		edgeVertIndices[2][1] = 3;

		edgeVertIndices[3][0] = 1;
		edgeVertIndices[3][1] = 2;

		edgeVertIndices[4][0] = 1;
		edgeVertIndices[4][1] = 3;

		edgeVertIndices[5][0] = 2;
		edgeVertIndices[5][1] = 3;

		faceVertIndices[0][0] = 0;
		faceVertIndices[0][1] = 1;
		faceVertIndices[0][2] = 2;

		faceVertIndices[1][0] = 0;
		faceVertIndices[1][1] = 3;
		faceVertIndices[1][2] = 1;

		faceVertIndices[2][0] = 1;
		faceVertIndices[2][1] = 3;
		faceVertIndices[2][2] = 2;

		faceVertIndices[3][0] = 2;
		faceVertIndices[3][1] = 3;
		faceVertIndices[3][2] = 0;

		faceEdgeIndices[0][0] = 0;
		faceEdgeIndices[0][1] = 3;
		faceEdgeIndices[0][2] = 1;

		faceEdgeIndices[1][0] = 2;
		faceEdgeIndices[1][1] = 4;
		faceEdgeIndices[1][2] = 0;

		faceEdgeIndices[2][0] = 4;
		faceEdgeIndices[2][1] = 5;
		faceEdgeIndices[2][2] = 3;

		faceEdgeIndices[3][0] = 5;
		faceEdgeIndices[3][1] = 2;
		faceEdgeIndices[3][2] = 1;

		if (type == Mapping::LengthScale) {
			m_Map = new LengthScaleTetMapping(pInitMesh);
		}
		else if (type == Mapping::Lagrange) {
			m_Map = new LagrangeCubicTetMapping(pInitMesh);
		}
		else {
			m_Map = new Q1TetMapping(pInitMesh);
		}
  }
	~TetDivider() {
	}
	void divideInterior();
  void createNewCells();
	void setupCoordMapping(const emInt verts[]);
	void getPhysCoordsFromParamCoords(const double uvw[], double xyz[]);
	void setPolyCoeffs(const double* xyz0, const double* xyz1, const double* xyz2,
			const double* xyz3, double uderiv0[3], double vderiv0[3],
			double wderiv0[3], double uderiv1[3], double vderiv1[3],
			double wderiv1[3], double uderiv2[3], double vderiv2[3],
			double wderiv2[3], double uderiv3[3], double vderiv3[3],
			double wderiv3[3]);
	void stuffTetsIntoOctahedron(emInt vertsNew[][4]);
};

#endif /* APPS_EXAMESH_TETDIVIDER_H_ */
