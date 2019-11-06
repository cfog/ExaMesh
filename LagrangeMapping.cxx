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
	setModalValues();
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

void LagrangeCubicMapping::computeTransformedCoords(const double uvw[3],
		double xyz[3]) const {
	xyz[0] = xyz[1] = xyz[2] = 0;
	const double& u = uvw[0];
	const double& v = uvw[1];
	const double& w = uvw[2];
	for (int ii = 0; ii < 3; ii++) {
		xyz[ii] = u
				* (Cu[ii] + u * (Cuu[ii] + u * Cuuu[ii] + v * Cuuv[ii] + w * Cuuw[ii])
						+ v * (Cuv[ii] + v * Cuvv[ii] + w * Cuvw[ii])
						+ w * (Cuw[ii] + w * Cuww[ii]))
							+ v * (Cv[ii] + v * (Cvv[ii] + v * Cvvv[ii] + w * Cvvw[ii])
											+ w * (Cvw[ii] + w * Cvww[ii]))
							+ w * (Cw[ii] + w * (Cww[ii] + w * Cwww[ii])) + C[ii];
	}
}

void LagrangeCubicTetMapping::setModalValues() {
	const double (*c)[3] = m_nodalValues;
	for (int ii = 0; ii < 3; ii++) {
		C[ii] = c[0][ii];

		Cu[ii] = -5.5 * c[0][ii] + c[1][ii] + 9 * c[4][ii] - 4.5 * c[5][ii];
		Cv[ii] = -5.5 * c[0][ii] + c[2][ii] + 9 * c[9][ii] - 4.5 * c[8][ii];
		Cw[ii] = -5.5 * c[0][ii] + c[3][ii] + 9 * c[10][ii] - 4.5 * c[11][ii];

		Cuu[ii] = 4.5 * (2 * c[0][ii] - c[1][ii] - 5 * c[4][ii] + 4 * c[5][ii]);
		Cvv[ii] = 4.5 * (2 * c[0][ii] - c[2][ii] - 5 * c[9][ii] + 4 * c[8][ii]);
		Cww[ii] = 4.5 * (2 * c[0][ii] - c[3][ii] - 5 * c[10][ii] + 4 * c[11][ii]);

		Cuv[ii] =
				4.5 * (4 * c[0][ii] - 5 * c[4][ii] + c[5][ii] - c[6][ii] - c[7][ii]
						+ c[8][ii]
								- 5 * c[9][ii]
								+ 6 * c[16][ii]);
		Cvw[ii] =
				4.5 * (4 * c[0][ii] - 5 * c[9][ii] + c[8][ii] - c[14][ii] - c[15][ii]
						+ c[11][ii]
								- 5 * c[10][ii]
								+ 6 * c[19][ii]);
		Cuw[ii] =
				4.5 * (4 * c[0][ii] - 5 * c[4][ii] + c[5][ii] - c[13][ii] - c[12][ii]
						+ c[11][ii]
								- 5 * c[10][ii]
								+ 6 * c[17][ii]);

		Cuuu[ii] = 4.5 * (-c[0][ii] + c[1][ii] + 3 * c[4][ii] - 3 * c[5][ii]);
		Cvvv[ii] = 4.5 * (-c[0][ii] + c[2][ii] + 3 * c[9][ii] - 3 * c[8][ii]);
		Cwww[ii] = 4.5 * (-c[0][ii] + c[3][ii] + 3 * c[10][ii] - 3 * c[11][ii]);

		Cuuv[ii] = 13.5
				* (-c[0][ii] + 2 * c[4][ii] - c[5][ii] + c[6][ii] + c[9][ii] - 2
						* c[16][ii]);
		Cvvw[ii] = 13.5
				* (-c[0][ii] + 2 * c[9][ii] - c[8][ii] + c[10][ii] + c[14][ii] - 2
						* c[19][ii]);
		Cuww[ii] = 13.5
				* (-c[0][ii] + 2 * c[10][ii] - c[11][ii] + c[4][ii] + c[13][ii] - 2
						* c[17][ii]);

		Cuvv[ii] = 13.5
				* (-c[0][ii] + c[4][ii] - c[8][ii] + c[7][ii] + 2 * c[9][ii] - 2
						* c[16][ii]);
		Cvww[ii] = 13.5
				* (-c[0][ii] + c[9][ii] - c[11][ii] + c[15][ii] + 2 * c[10][ii] - 2
						* c[19][ii]);
		Cuuw[ii] = 13.5
				* (-c[0][ii] + c[10][ii] - c[5][ii] + c[12][ii] + 2 * c[4][ii] - 2
						* c[17][ii]);

		Cuvw[ii] = 27
				* (-c[0][ii] + c[4][ii] + c[9][ii] + c[10][ii] - c[16][ii] - c[17][ii]
						- c[19][ii]
						+ c[18][ii]);
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
