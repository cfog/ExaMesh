/*
 * PyrDivider.h
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#ifndef APPS_EXAMESH_PYRDIVIDER_H_
#define APPS_EXAMESH_PYRDIVIDER_H_

#include "CellDivider.h"
#include "ExaMesh.h"

class PyrDivider: public CellDivider {
	double xyzOffset[3], uVec[3], vVec[3], uvVec[3], xyzApex[3];
public:
	PyrDivider(UMesh *pInitMesh, const int segmentsPerEdge,
			const Mapping::MappingType type = Mapping::Uniform)
      :
			CellDivider(pInitMesh, segmentsPerEdge) {
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

		if (type == Mapping::LengthScale) {
			// TODO: Must fix this
			m_Map = new UniformPyramidMapping(pInitMesh);
		}
		else if (type == Mapping::Lagrange) {
			m_Map = new LagrangeCubicPyramidMapping(pInitMesh);
		}
		else {
			m_Map = new UniformPyramidMapping(pInitMesh);
		}
  }
	~PyrDivider() {
	}
	void divideInterior();
  void createNewCells();
	void setupCoordMapping(const emInt verts[]);
	void getPhysCoordsFromParamCoords(const double uvw[], double xyz[]);
};

#endif /* APPS_EXAMESH_PYRDIVIDER_H_ */
