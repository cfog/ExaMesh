/*
 * Mapping.cxx
 *
 *  Created on: Oct. 1, 2019
 *      Author: cfog
 */

#include <assert.h>
#include <ExaMesh.h>

#include "Mapping.h"

LagrangeMapping::LagrangeMapping(const ExaMesh* const EM,
		const int nVals) :
		Mapping(EM), m_numValues(nVals) {
	m_nodalValues = new double[nVals][3];
}

LagrangeMapping::~LagrangeMapping() {
	delete[] m_nodalValues;
}

void LagrangeMapping::setupCoordMapping(const emInt verts[]) {
	double coords[m_numValues][3];
	for (int ii = 0; ii < 20; ii++) {
		coords[ii][0] = m_pMesh->getX(verts[ii]);
		coords[ii][1] = m_pMesh->getY(verts[ii]);
		coords[ii][2] = m_pMesh->getZ(verts[ii]);
	}
	setNodalValues(coords);
}



void LagrangeMapping::setNodalValues(double inputValues[][3]) {
	for (int ii = 0; ii < m_numValues; ii++) {
		m_nodalValues[ii][0] = inputValues[ii][0];
		m_nodalValues[ii][1] = inputValues[ii][1];
		m_nodalValues[ii][2] = inputValues[ii][2];
	}
}

void LagrangeMapping::computeTransformedCoords(const double uvw[3], double xyz[3]) const {
	xyz[0] = xyz[1] = xyz[2] = 0;
	for (int ii = 0; ii < m_numValues; ii++) {
		double basis = computeBasisFunction(ii, uvw);
		xyz[0] += basis * m_nodalValues[ii][0];
		xyz[1] += basis * m_nodalValues[ii][1];
		xyz[2] += basis * m_nodalValues[ii][2];
	}
}

double LagrangeCubicTetMapping::computeBasisFunction(const int whichFunc,
		const double uvw[3]) const {
	// These basis functions are in the same order as the nodes, which in
	// turn are ordered according the CGNS numbering system, with base 0
	// instead of 1.
	double b1 = uvw[0];
	double b2 = uvw[1];
	double b3 = uvw[2];
	double b0 = 1 - b3 - b1 - b2;
	switch (whichFunc) {
		// For the nodes
		case 0:
			return 4.5 * b0 * (b0 - 1. / 3) * (b0 - 2. / 3);
		case 1:
			return 4.5 * b1 * (b1 - 1. / 3) * (b1 - 2. / 3);
		case 2:
			return 4.5 * b2 * (b2 - 1. / 3) * (b2 - 2. / 3);
		case 3:
			return 4.5 * b3 * (b3 - 1. / 3) * (b3 - 2. / 3);

			// Points on the edge from 0 to 1
		case 4:
			return 13.5 * b0 * b1 * (b0 - 1. / 3);
		case 5:
			return 13.5 * b0 * b1 * (b1 - 1. / 3);

			// Points on the edge from 1 to 2
		case 6:
			return 13.5 * b1 * b2 * (b1 - 1. / 3);
		case 7:
			return 13.5 * b1 * b2 * (b2 - 1. / 3);

			// Points on the edge from 2 to 0
		case 8:
			return 13.5 * b2 * b0 * (b2 - 1. / 3);
		case 9:
			return 13.5 * b2 * b0 * (b0 - 1. / 3);

			// Points on the edge from 0 to 3
		case 10:
			return 13.5 * b0 * b3 * (b0 - 1. / 3);
		case 11:
			return 13.5 * b0 * b3 * (b3 - 1. / 3);

			// Points on the edge from 1 to 3
		case 12:
			return 13.5 * b1 * b3 * (b1 - 1. / 3);
		case 13:
			return 13.5 * b1 * b3 * (b3 - 1. / 3);

			// Points on the edge from 2 to 3
		case 14:
			return 13.5 * b2 * b3 * (b2 - 1. / 3);
		case 15:
			return 13.5 * b2 * b3 * (b3 - 1. / 3);

		case 16:
			// Point on the face 012
			return 27 * b0 * b1 * b2;

		case 17:
			// Point on the face 013
			return 27 * b0 * b1 * b3;

		case 18:
			// Point on the face 123
			return 27 * b1 * b2 * b3;

		case 19:
			// Point on the face 203
			return 27 * b2 * b0 * b3;

		default:
			assert(0);
			return 0;
	}
}
