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

void LagrangeCubicTetMapping::computeTransformedCoords(const double uvw[3],
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

void LagrangeCubicPyramidMapping::setModalValues() {
	for (int ii = 0; ii < 3; ii++) {
		// Some aliases for these coordinate values, which will greatly shorten the
		// code that computes the modal values.
		double& c0 = m_nodalValues[0][ii];
		double& c1 = m_nodalValues[1][ii];
		double& c2 = m_nodalValues[2][ii];
		double& c3 = m_nodalValues[3][ii];
		double& c4 = m_nodalValues[4][ii];
		double& c5 = m_nodalValues[5][ii];
		double& c6 = m_nodalValues[6][ii];
		double& c7 = m_nodalValues[7][ii];
		double& c8 = m_nodalValues[8][ii];
		double& c9 = m_nodalValues[9][ii];
		double& c10 = m_nodalValues[10][ii];
		double& c11 = m_nodalValues[11][ii];
		double& c12 = m_nodalValues[12][ii];
		double& c13 = m_nodalValues[13][ii];
		double& c14 = m_nodalValues[14][ii];
		double& c15 = m_nodalValues[15][ii];
		double& c16 = m_nodalValues[16][ii];
		double& c17 = m_nodalValues[17][ii];
		double& c18 = m_nodalValues[18][ii];
		double& c19 = m_nodalValues[19][ii];
		double& c20 = m_nodalValues[20][ii];
		double& c21 = m_nodalValues[21][ii];
		double& c22 = m_nodalValues[22][ii];
		double& c23 = m_nodalValues[23][ii];
		double& c24 = m_nodalValues[24][ii];
		double& c25 = m_nodalValues[25][ii];
		double& c26 = m_nodalValues[26][ii];
		double& c27 = m_nodalValues[27][ii];
		double& c28 = m_nodalValues[28][ii];
		double& c29 = m_nodalValues[29][ii];

		Apex[ii] = c4;

		C[ii] =
				1 / 256. * ((c0 + c2 + c1 + c3)
						- 9 * (c11 + c10 + c9 + c8 + c7 + c5 + c12 + c6)
										+ 81 * (c21 + c22 + c23 + c24));

		Cu[ii] = +1 / 256.
				* (27 * (-c6 + c5 - c9 + c10) + 9 * (-c8 - c7 + c12 + c11)
						+ (-c3 + c2 + c1 - c0) + 243 * (-c21 + c22 + c23 - c24));

		Cv[ii] = 1. / 256
				* (9 * (c6 + c5 - c9 - c10) + 27 * (-c8 + c7 + c12 - c11) + c3 + c2 - c1
						- c0
						+ 243 * (-c21 - c22 + c23 + c24));

		Cw[ii] = (+37. / 256 * (-c3 - c2 - c1 - c0)
				+ 45. / 256 * (c6 + c8 + c7 + c5 + c12 + c9 + c11 + c10)
				+ 9. / 16 * (+c13 + c15 + c17 + c19) + c4
				+ 9. / 8 * (-c14 - c16 - c18 - c20)
				+ 405. / 256 * (-c21 - c22 - c23 - c24) + 27. / 4 * c29);

		Cuu[ii] = 9. / 256
				* (c6 + c5 + c9 - c3 - c2 - c1 + c10 - c0 + 9
						* (+c8 + c7 + c12 + c11 - c21 - c22 - c23 - c24));

		Cuv[ii] = 63. / 256 * (+c3 - c2 + c1 - c0)
				+ (9. / 8 * (+4 * (+c13 - c15 + c17 - c19)
						+ 5 * (-c14 + c16 - c18 + c20)))
							- 27. / 256 * (32
									* (c13 - c15 + c17 - c19 - c14 + c16 - c18 + c20)
															+ c12
															- c11
															+ c8
															- c7 - c10 - c6
															+ c5 + c9 + 27 * (-c21 + c22 - c23 + c24));

		Cvv[ii] = 1. / 256
				* (81 * (c6 + c5 + c9 + c10 - c21 - c22 - c23 - c24) + 9
						* (+c12 + c8 + c7 - c3 - c2 - c1 + c11 - c0));

		Cvw[ii] = 9. / 128
				* ((-c3 - c2 + c1 + c0) + 3 * (+c8 - c7 - c12 + c11)
						+ 27 * (c21 + c22 - c23 - c24) + 7 * (c6 + c5 - c10 - c9)
						+ 8 * (-c15 - c13 + c19 + c17) + 16 * (+c14 + c16 - c18 - c20)
						+ 48 * (-c25 + c27));

		Cww[ii] = (63. / 256 * (-c6 - c8 - c7 - c5 - c12 - c9 - c11 - c10)
				+ 135. / 256 * (+c3 + c2 + c1 + c0)
				+ 567. / 256 * (+c21 + c22 + c23 + c24)
				+ 4.5 * (+c14 + c16 + c18 + c20 - c4) + 2.25 * (-c13 - c15 - c17 - c19)
								- 13.5 * c29);

		Cuw[ii] = (+9 / 128.
				* (3 * (c6 - c5 - c10 + c9) + 7 * (-c8 - c7 + c12 + c11)
						+ 27 * (c21 - c22 - c23 + c24))
								+ 9. / 16 * (-c13 + c17 - c19 + c15)
								+ 27. / 8 * (c26 - c28)
								+ 9. / 128 * (c0 + c3 - c1 - c2)
								+ 1.125 * (c14 - c16 - c18 + c20));

		Cuuu[ii] = 9. / 256
				* (3 * (c6 - c5 + c9 - c10) + 9 * (+c8 + c7 - c12 - c11)
						+ 27 * (+c21 - c22 - c23 + c24) + c3
						- c2 - c1
						+ c0);

		Cuuv[ii] = 27. / 256
				* (5 * (c6 + c5 - c9 + c3 + c2 - c1 - c10 - c0) + 9
						* (+c8 - c7 - c12 + c11 + c21 + c22 - c23 - c24)
						+ 16 * (+c13 + c15 - c17 - c19 - 2 * c25 + 2 * c27));

		Cuvv[ii] = (27. / 256
				* (9 * (-c5 - c10 + c6 + c9 + c21 - c22 - c23 + c24) + 5
						* (c12 - c3 + c2 + c1 + c11 - c8 - c7 - c0)
						+ 32 * (-c28 + c26) + 16 * (+c13 - c15 - c17 + c19)));

		Cvvv[ii] = 1 / 256.
				* (+81 * (-c6 - c5 + c9 + c10) + 27 * (c11 - c12 + c8 - c7)
						+ 9 * (c0 + c1 - c2 - c3) + 243 * (c21 + c22 - c23 - c24));

		Cvvw[ii] = 27. / 256
				* ((-c12 - c11 - c8 - c7 + c3 + c2 + c1 + c0) + 9
						* (-c10 + c21 + c22 + c23 + c24 - c6 - c5 - c9)
						+ 32 * (c25 - 2 * c29 + c27));

		Cvww[ii] = 27. / 256
				* (-c8 + c7 + c12 - c11 + 3 * (+c3 + c2 - c1 - c0)
						+ 5 * (+c9 + c10 - c6 - c5) + 9 * (-c21 - c22 + c23 + c24)
						+ 16 * (+c13 + c15 - c19 - c17)
						+ 32 * (-c14 - c16 + c18 + c20 + c25 - c27));

		Cwww[ii] = 9. / 256
				* (3 * (c6 + c8 + c7 + c5 + c12 + c9 + c11 + c10) - 11
						* (c3 + c2 + c1 + c0)
						- 27 * (+c21 + c22 + c23 + c24)
						+ 48 * (+c13 + c15 + c17 + c19)
						+ 96 * (-c20 - c18 - c16 - c14 + 2 * c29) + 128 * c4);

		Cuww[ii] = (27 / 256.
				* (+(-c6 + c5 - c9 + c10) + 3 * (-c3 + c2 + c1 - c0)
						+ 5 * (-c11 + c8 + c7 - c12) + 9 * (-c21 + c22 + c23 - c24))
								+ 27. / 16 * (c13 - c15 - c17 + c19)
								+ 27. / 8 * (-c26 - c14 + c28 + c16 + c18 - c20));

		Cuuw[ii] = 27. / 256
				* (9 * (-c12 - c11 - c8 - c7 + c21 + c22 + c23 + c24) - c10 - c6 - c5
						- c9
						+ c3 + c2 + c1 + c0 + 32 * (+c28 + c26 - 2 * c29));

		Cuvw[ii] =
		(27. / 256
				* (32 * (c13 - c15 + c17 - c19 - c14 + c16 - c18 + c20) + c12
						- c11
											+ c8
											- c7 - c10 - c6
											+ c5 + c9
						+ 27 * (-c21 + c22 - c23 + c24))
		+ 225. / 256 * (c3 - c2 + c1 - c0));

		CuvOverw[ii] =
				(0.25 * (9 * (-c14 + c16 - c18 + c20) + 4.5 * (c13 - c15 + c17 - c19)
									+ c3
									- c2
									+ c1
									- c0));

		Cu2vOverw[ii] =
				9. / 16 * (3 * (+c13 + c15 - c17 - c19 - 2 * c25 + 2 * c27) + c6 + c5
						- c9 - c10
										+ c3 + c2
										- c1 - c0);

		Cuv2Overw[ii] =
				9. / 16 * (3 * (c13 - c15 - c17 + c19) + 6 * (c26 - c28) - c8 - c7 + c12
										+ c11
										- c3
										+ c2 + c1
										- c0);

		Cu3vOverw[ii] = 9. / 256
				* (27 * (-c8 + c7 - c12 + c11) + 81 * (+c21 - c22 + c23 - c24)
						+ 3 * (+c6 - c5 - c9 + c10)
						- c3
						+ c2
						- c1
						+ c0);

		Cu2v2Overw[ii] = (27. / 16
				* (4 * c29 + c13 + c15 + c17 + c19 - 2 * c26 - 2 * c28 - 2 * c25
						- 2 * c27)
							+ 243. / 256 * (+c8 + c7 + c12 + c11 - c21 - c22 - c23 - c24 + c6
															+ c5 + c9 + c10
																			- c3 - c2 - c1 - c0));

		Cuv3Overw[ii] = 9. / 256
				* (3 * (-c12 + c11 - c8 + c7) + 27 * (-c9 + c10 + c6 - c5)
						+ 81 * (c21 - c22 + c23 - c24)
						- c3
						+ c2
						- c1
						+ c0);

		Cu2v2Overw2[ii] = (27. / 16
				* (4 * c29 + c13 + c15 + c17 + c19 - 2 * c26 - 2 * c28 - 2 * c25
						- 2 * c27)
							+ 81. / 128 * (+c6 + c8 + c7 + c5 + c12 + c9 + c11 + c10 - c3 - c2
															- c1 - c0 - c21 - c22 - c23 - c24));

		Cu3v2Overw2[ii] =
				81. / 256 * (3 * (-c6 + c5 - c9 + c10 - c21 + c22 + c23 - c24) - c8 - c7
						+ c12
											- c3
											+ c2 + c1 + c11
											- c0);

		Cu2v3Overw2[ii] = 81. / 256
				* (c6 + c5 - c9 + c3 + c2 - c1 - c10 - c0 + 3
						* (-c8 + c7 + c12 - c11 - c21 - c22 + c23 + c24));

		Cu3v3Overw3[ii] =
				81. / 256 * (3 * (-c6 + c8 - c7 + c5 + c12 + c9 - c11 - c10) + c3 - c2
						+ c1
											- c0
											+ 9 * (-c21 + c22 - c23 + c24));
	}
}

void LagrangeCubicPyramidMapping::computeTransformedCoords(const double uvw[3],
		double xyz[3]) const {
	xyz[0] = xyz[1] = xyz[2] = 0;
	double u = uvw[0];
	double v = uvw[1];
	const double& w = uvw[2];

	// Now we need to re-map u and v.  The rest of examesh assumes that 0 <= u,v, <= (1-w),
	// whereas this routine uses basis functions that assume -(1-w) <= u,v <= (1-w), so
	// we need to recompute u and v.
	u = 2 * u - (1 - w);
	v = 2 * v - (1 - w);
	assert(u >= 0 && v >= 0 && w >= 0 && w <= 1);
	assert(fabs(u) <= (1-w));
	assert(fabs(v) <= (1-w));

	if (w == 1) {
		xyz[0] = Apex[0];
		xyz[1] = Apex[1];
		xyz[2] = Apex[2];
		return;
	}
	double frac = u * v / (w - 1);

	for (int ii = 0; ii < 3; ii++) {
		// First the polynomial terms.
		xyz[ii] = u
				* (Cu[ii] + u * (Cuu[ii] + u * Cuuu[ii] + v * Cuuv[ii] + w * Cuuw[ii])
						+ v * (Cuv[ii] + v * Cuvv[ii] + w * Cuvw[ii])
						+ w * (Cuw[ii] + w * Cuww[ii]))
							+ v * (Cv[ii] + v * (Cvv[ii] + v * Cvvv[ii] + w * Cvvw[ii])
											+ w * (Cvw[ii] + w * Cvww[ii]))
							+ w * (Cw[ii] + w * (Cww[ii] + w * Cwww[ii])) + C[ii];
		// Now the rational terms
		xyz[ii] += frac
				* (CuvOverw[ii] + u
						* (Cu2vOverw[ii] + Cu3vOverw[ii] * u + Cu2v2Overw[ii] * v)
						+ v * (Cuv2Overw[ii] + Cuv3Overw[ii] * v)
						+ frac * (Cu2v2Overw2[ii] + Cu3v2Overw2[ii] * u
											+ Cu2v3Overw2[ii] * v + Cu3v3Overw3[ii] * frac));
	}
}

// Function is deprecated.
double LagrangeCubicPyramidMapping::computeBasisFunction(const int whichFunc,
		const double uvw[3]) const {
	return 0;
}


