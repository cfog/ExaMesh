/*
 * Mapping.cxx
 *
 *  Created on: Oct. 1, 2019
 *      Author: cfog
 */

#include <assert.h>
#include <ExaMesh.h>

#include "Mapping.h"

LagrangeMapping::LagrangeMapping(const ExaMesh* const EM, const int nVals) :
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
				(27. / 256 * (32 * (c13 - c15 + c17 - c19 - c14 + c16 - c18 + c20) + c12
						- c11
											+ c8
											- c7 - c10 - c6
											+ c5 + c9 + 27 * (-c21 + c22 - c23 + c24))
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
											+ 243. / 256 * (+c8 + c7 + c12 + c11 - c21 - c22 - c23
																			- c24
																			+ c6 + c5 + c9 + c10
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
												+ 81. / 128 * (+c6 + c8 + c7 + c5 + c12 + c9 + c11 + c10
														- c3 - c2 - c1 - c0 - c21 - c22 - c23 - c24));

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

void LagrangeCubicPrismMapping::setModalValues() {
	for (int ii = 0; ii < 3; ii++) {
		{
			// Bottom
			double& c0 = m_nodalValues[0][ii];
			double& c1 = m_nodalValues[1][ii];
			double& c2 = m_nodalValues[2][ii];
			double& c3 = m_nodalValues[6][ii];
			double& c4 = m_nodalValues[7][ii];
			double& c5 = m_nodalValues[8][ii];
			double& c6 = m_nodalValues[9][ii];
			double& c7 = m_nodalValues[10][ii];
			double& c8 = m_nodalValues[11][ii];
			double& c9 = m_nodalValues[24][ii];

			C0[ii] = +c0;
			Cu0[ii] = +(c1 + 9 * c3 - 4.5 * c4 - 5.5 * c0);
			Cv0[ii] = +(c2 - 4.5 * c7 + 9 * c8 - 5.5 * c0);
			Cuu0[ii] = +4.5 * (2 * c0 - c1 - 5 * c3 + 4 * c4);
			Cuv0[ii] = +4.5 * (4 * c0 + c4 + c7 - c5 - c6 + 6 * c9 - 5 * c3 - 5 * c8);
			Cvv0[ii] = +4.5 * (2 * c0 - c2 + 4 * c7 - 5 * c8);
			Cuuu0[ii] = +4.5 * (-c0 + c1 + 3 * c3 - 3 * c4);
			Cuuv0[ii] = +13.5 * (c5 - c0 - 2 * c9 + 2 * c3 + c8 - c4);
			Cuvv0[ii] = +13.5 * (-2 * c9 - c0 + 2 * c8 - c7 + c6 + c3);
			Cvvv0[ii] = +4.5 * (-c0 + c2 - 3 * c7 + 3 * c8);
		}
		{
			// Lower interior layer
			double& c0 = m_nodalValues[12][ii];
			double& c1 = m_nodalValues[14][ii];
			double& c2 = m_nodalValues[16][ii];
			double& c3 = m_nodalValues[25][ii];
			double& c4 = m_nodalValues[26][ii];
			double& c5 = m_nodalValues[29][ii];
			double& c6 = m_nodalValues[30][ii];
			double& c7 = m_nodalValues[33][ii];
			double& c8 = m_nodalValues[34][ii];
			double& c9 = m_nodalValues[38][ii];

			C1[ii] = +c0;
			Cu1[ii] = +(c1 + 9 * c3 - 4.5 * c4 - 5.5 * c0);
			Cv1[ii] = +(c2 - 4.5 * c7 + 9 * c8 - 5.5 * c0);
			Cuu1[ii] = +4.5 * (2 * c0 - c1 - 5 * c3 + 4 * c4);
			Cuv1[ii] = +4.5 * (4 * c0 + c4 + c7 - c5 - c6 + 6 * c9 - 5 * c3 - 5 * c8);
			Cvv1[ii] = +4.5 * (2 * c0 - c2 + 4 * c7 - 5 * c8);
			Cuuu1[ii] = +4.5 * (-c0 + c1 + 3 * c3 - 3 * c4);
			Cuuv1[ii] = +13.5 * (c5 - c0 - 2 * c9 + 2 * c3 + c8 - c4);
			Cuvv1[ii] = +13.5 * (-2 * c9 - c0 + 2 * c8 - c7 + c6 + c3);
			Cvvv1[ii] = +4.5 * (-c0 + c2 - 3 * c7 + 3 * c8);
		}
		{
			// Upper interior layer
			double& c0 = m_nodalValues[13][ii];
			double& c1 = m_nodalValues[15][ii];
			double& c2 = m_nodalValues[17][ii];
			double& c3 = m_nodalValues[28][ii];
			double& c4 = m_nodalValues[27][ii];
			double& c5 = m_nodalValues[32][ii];
			double& c6 = m_nodalValues[31][ii];
			double& c7 = m_nodalValues[36][ii];
			double& c8 = m_nodalValues[35][ii];
			double& c9 = m_nodalValues[39][ii];

			C2[ii] = +c0;
			Cu2[ii] = +(c1 + 9 * c3 - 4.5 * c4 - 5.5 * c0);
			Cv2[ii] = +(c2 - 4.5 * c7 + 9 * c8 - 5.5 * c0);
			Cuu2[ii] = +4.5 * (2 * c0 - c1 - 5 * c3 + 4 * c4);
			Cuv2[ii] = +4.5 * (4 * c0 + c4 + c7 - c5 - c6 + 6 * c9 - 5 * c3 - 5 * c8);
			Cvv2[ii] = +4.5 * (2 * c0 - c2 + 4 * c7 - 5 * c8);
			Cuuu2[ii] = +4.5 * (-c0 + c1 + 3 * c3 - 3 * c4);
			Cuuv2[ii] = +13.5 * (c5 - c0 - 2 * c9 + 2 * c3 + c8 - c4);
			Cuvv2[ii] = +13.5 * (-2 * c9 - c0 + 2 * c8 - c7 + c6 + c3);
			Cvvv2[ii] = +4.5 * (-c0 + c2 - 3 * c7 + 3 * c8);
		}
		{
			// Top
			double& c0 = m_nodalValues[3][ii];
			double& c1 = m_nodalValues[4][ii];
			double& c2 = m_nodalValues[5][ii];
			double& c3 = m_nodalValues[18][ii];
			double& c4 = m_nodalValues[19][ii];
			double& c5 = m_nodalValues[20][ii];
			double& c6 = m_nodalValues[21][ii];
			double& c7 = m_nodalValues[22][ii];
			double& c8 = m_nodalValues[23][ii];
			double& c9 = m_nodalValues[37][ii];

			C3[ii] = +c0;
			Cu3[ii] = +(c1 + 9 * c3 - 4.5 * c4 - 5.5 * c0);
			Cv3[ii] = +(c2 - 4.5 * c7 + 9 * c8 - 5.5 * c0);
			Cuu3[ii] = +4.5 * (2 * c0 - c1 - 5 * c3 + 4 * c4);
			Cuv3[ii] = +4.5 * (4 * c0 + c4 + c7 - c5 - c6 + 6 * c9 - 5 * c3 - 5 * c8);
			Cvv3[ii] = +4.5 * (2 * c0 - c2 + 4 * c7 - 5 * c8);
			Cuuu3[ii] = +4.5 * (-c0 + c1 + 3 * c3 - 3 * c4);
			Cuuv3[ii] = +13.5 * (c5 - c0 - 2 * c9 + 2 * c3 + c8 - c4);
			Cuvv3[ii] = +13.5 * (-2 * c9 - c0 + 2 * c8 - c7 + c6 + c3);
			Cvvv3[ii] = +4.5 * (-c0 + c2 - 3 * c7 + 3 * c8);
		}
	}
}

void LagrangeCubicPrismMapping::computeTransformedCoords(const double uvw[3],
		double xyz[3]) const {
	const double& u = uvw[0];
	const double& v = uvw[1];
	const double& w = uvw[2];

	double segment0 = -4.5 * (w - 1) * (w - 2. / 3) * (w - 1. / 3);
	double segment1 = 13.5 * (w - 1) * (w - 2. / 3) * w;
	double segment2 = -13.5 * (w - 1) * (w - 1. / 3) * w;
	double segment3 = 4.5 * (w - 2. / 3) * (w - 1. / 3) * w;

	for (int ii = 0; ii < 3; ii++) {
		double val0 = C0[ii]
				+ u * (Cu0[ii] + u * (Cuu0[ii] + u * Cuuu0[ii] + v * Cuuv0[ii])
								+ v * (Cuv0[ii] + v * Cuvv0[ii]))
				+ v * (Cv0[ii] + v * (Cvv0[ii] + v * Cvvv0[ii]));
		double val1 = C1[ii]
				+ u * (Cu1[ii] + u * (Cuu1[ii] + u * Cuuu1[ii] + v * Cuuv1[ii])
								+ v * (Cuv1[ii] + v * Cuvv1[ii]))
				+ v * (Cv1[ii] + v * (Cvv1[ii] + v * Cvvv1[ii]));
		double val2 = C2[ii]
				+ u * (Cu2[ii] + u * (Cuu2[ii] + u * Cuuu2[ii] + v * Cuuv2[ii])
								+ v * (Cuv2[ii] + v * Cuvv2[ii]))
				+ v * (Cv2[ii] + v * (Cvv2[ii] + v * Cvvv2[ii]));
		double val3 = C3[ii]
				+ u * (Cu3[ii] + u * (Cuu3[ii] + u * Cuuu3[ii] + v * Cuuv3[ii])
								+ v * (Cuv3[ii] + v * Cuvv3[ii]))
				+ v * (Cv3[ii] + v * (Cvv3[ii] + v * Cvvv3[ii]));
		xyz[ii] = segment0 * val0 + segment1 * val1 + segment2 * val2
							+ segment3 * val3;
	}
}

//void LagrangeCubicHexMapping::setModalValues() {
//	for (int ii = 0; ii < 3; ii++) {
//		{
//			// Bottom
//			double& c00 = m_nodalValues[0][ii];
//			double& c10 = m_nodalValues[8][ii];
//			double& c20 = m_nodalValues[9][ii];
//			double& c30 = m_nodalValues[1][ii];
//			double& c01 = m_nodalValues[15][ii];
//			double& c11 = m_nodalValues[32][ii];
//			double& c21 = m_nodalValues[33][ii];
//			double& c31 = m_nodalValues[10][ii];
//			double& c02 = m_nodalValues[14][ii];
//			double& c12 = m_nodalValues[35][ii];
//			double& c22 = m_nodalValues[34][ii];
//			double& c32 = m_nodalValues[11][ii];
//			double& c03 = m_nodalValues[3][ii];
//			double& c13 = m_nodalValues[13][ii];
//			double& c23 = m_nodalValues[12][ii];
//			double& c33 = m_nodalValues[2][ii];
//
//			C0[ii] = c00;
//
//			Cu0[ii] = c30 - 11 / 2 * c00 + 9 * c10 - 9 / 2 * c20;
//			Cv0[ii] = 9 * c01 - 9 / 2 * c02 + c03 - 11 / 2 * c00;
//
//			Cuu0[ii] = 4.5 * (-c30 + 2 * c00 - 5 * c10 + 4 * c20);
//			Cuv0[ii] = 99. / 4 * (-2 * c01 - 2 * c10 + c02 + c20)
//					+ 81. / 4 * (4 * c11 - 2 * c12 - 2 * c21 + c22)
//					+ 9 * (c13 + c31 - 0.5 * (c23 + c32)) + c33
//									- 5.5 * (c03 + c30)
//									+ 121. / 4 * c00;
//			Cvv0[ii] = 4.5 * (-5 * c01 + 4 * c02 - c03 + 2 * c00);
//
//			Cuuu0[ii] = 4.5 * (c30 - c00 + 3 * (c10 - c20));
//			Cuuv0[ii] = 81 * c01 - 81 * c22 - 81. / 2 * (c02 + c31) + 81. / 4 * c32
//									+ 18 * c23 + 9 * c03
//									- 9. / 2 * c33
//									+ 495. / 4 * c10
//									- 405. / 2 * c11
//									+ 405. / 4 * c12
//									- 45. / 2 * c13
//									+ 162 * c21
//									- 99 * c20
//									+ 99. / 4 * c30
//									- 99. / 2 * c00;
//			Cuvv0[ii] = 162 * c12 + 81 * c10 - 405. / 2 * c11 - 81. / 2 * c13
//									- 81. / 2 * c20
//									+ 405. / 4 * c21
//									- 81 * c22
//									+ 81. / 4 * c23 + 9 * c30
//									- 45 / 2 * c31
//									+ 18 * c32
//									- 9. / 2 * c33
//									+ 495. / 4 * c01
//									- 99 * c02
//									+ 99. / 4 * c03
//									- 99. / 2 * c00;
//			Cvvv0[ii] = 9. / 2 * (3 * (c01 - c02) + c03 - c00);
//
//			Cuuuv0[ii] = -81. / 2 * c01 + 81. / 4 * c02 - 9. / 2 * c03
//										- 297. / 4 * c10
//										+ 243. / 2 * c11
//										- 243. / 4 * c12
//										+ 27. / 2 * c13
//										- 243. / 2 * c21
//										+ 243. / 4 * c22
//										- 27. / 2 * c23
//										+ 81. / 2 * c31
//										- 81. / 4 * c32
//										+ 9. / 2 * c33 + 297 / 4 * c20
//										- 99. / 4 * c30
//										+ 99. / 4 * c00;
//			Cuuvv0[ii] = -405 * c12 - 405. / 2 * c10 + 2025. / 4 * c11
//										+ 405. / 4 * c13 + 162 * c20
//										- 405 * c21
//										+ 324 * c22
//										- 81 * c23 - 81. / 2 * c30
//										+ 405. / 4 * c31
//										- 81 * c32
//										+ 81. / 4 * c33
//										- 405. / 2 * c01
//										+ 162 * c02
//										- 81. / 2 * c03
//										+ 81 * c00;
//			Cuvvv0[ii] = -243. / 2 * c12 + 81. / 2 * c13 + 81. / 4 * c20
//					- 243. / 4 * c21
//										+ 243. / 4 * c22
//										- 81. / 4 * c23 - 81. / 2 * c10
//										+ 243. / 2 * c11
//										- 9. / 2 * c30
//										+ 27. / 2 * c31
//										- 27. / 2 * c32
//										+ 9. / 2 * c33
//										- 297. / 4 * c01
//										+ 297. / 4 * c02
//										- 99. / 4 * c03
//										+ 99. / 4 * c00;
//
//			Cuuuvv0[ii] = 243 * c12 + 243 / 2 * c10 - 1215 / 4 * c11 - 243 / 4 * c13
//										- 243 / 2 * c20
//										+ 1215 / 4 * c21
//										- (243 * c22)
//										+ 243 / 4 * c23 + 81 / 2 * c30
//										- 405 / 4 * c31
//										+ 81 * c32
//										- 81 / 4 * c33
//										+ 405 / 4 * c01
//										- 81 * c02
//										+ 81 / 4 * c03
//										- 81 / 2 * c00;
//			Cuuvvv0[ii] = 1215 / 4 * c12 - 405 / 4 * c13 - 81 * c20 + (243 * c21)
//					- (243 * c22)
//										+ (81 * c23) + 405 / 4 * c10
//										- 1215 / 4 * c11
//										+ 81 / 4 * c30
//										- 243 / 4 * c31
//										+ 243 / 4 * c32
//										- 81 / 4 * c33
//										+ 243 / 2 * c01
//										- 243 / 2 * c02
//										+ 81 / 2 * c03
//										- 81 / 2 * c00;
//
//			Cuuuvvv0[ii] = -729 / 4 * c12 + 243 / 4 * c13 + 243 / 4 * c20
//					- 729 / 4 * c21
//											+ 729 / 4 * c22
//											- 243 / 4 * c23 - 243 / 4 * c10
//											+ 729 / 4 * c11
//											- 81 / 4 * c30
//											+ 243 / 4 * c31
//											- 243 / 4 * c32
//											+ 81 / 4 * c33
//											- 243 / 4 * c01
//											+ 243 / 4 * c02
//											- 81 / 4 * c03
//											+ 81 / 4 * c00;
//		}
//
//		{
//			// Lower interior layer
//			double& c00 = m_nodalValues[16][ii];
//			double& c10 = m_nodalValues[36][ii];
//			double& c20 = m_nodalValues[37][ii];
//			double& c30 = m_nodalValues[18][ii];
//			double& c01 = m_nodalValues[49][ii];
//			double& c11 = m_nodalValues[56][ii];
//			double& c21 = m_nodalValues[57][ii];
//			double& c31 = m_nodalValues[40][ii];
//			double& c02 = m_nodalValues[48][ii];
//			double& c12 = m_nodalValues[59][ii];
//			double& c22 = m_nodalValues[48][ii];
//			double& c32 = m_nodalValues[41][ii];
//			double& c03 = m_nodalValues[22][ii];
//			double& c13 = m_nodalValues[45][ii];
//			double& c23 = m_nodalValues[44][ii];
//			double& c33 = m_nodalValues[20][ii];
//
//			C0[ii] = c00;
//
//			Cu0[ii] = c30 - 11 / 2 * c00 + 9 * c10 - 9 / 2 * c20;
//			Cv0[ii] = 9 * c01 - 9 / 2 * c02 + c03 - 11 / 2 * c00;
//
//			Cuu0[ii] = 4.5 * (-c30 + 2 * c00 - 5 * c10 + 4 * c20);
//			Cuv0[ii] = 99. / 4 * (-2 * c01 - 2 * c10 + c02 + c20)
//					+ 81. / 4 * (4 * c11 - 2 * c12 - 2 * c21 + c22)
//					+ 9 * (c13 + c31 - 0.5 * (c23 + c32)) + c33
//									- 5.5 * (c03 + c30)
//									+ 121. / 4 * c00;
//			Cvv0[ii] = 4.5 * (-5 * c01 + 4 * c02 - c03 + 2 * c00);
//
//			Cuuu0[ii] = 4.5 * (c30 - c00 + 3 * (c10 - c20));
//			Cuuv0[ii] = 81 * c01 - 81 * c22 - 81. / 2 * (c02 + c31) + 81. / 4 * c32
//									+ 18 * c23 + 9 * c03
//									- 9. / 2 * c33
//									+ 495. / 4 * c10
//									- 405. / 2 * c11
//									+ 405. / 4 * c12
//									- 45. / 2 * c13
//									+ 162 * c21
//									- 99 * c20
//									+ 99. / 4 * c30
//									- 99. / 2 * c00;
//			Cuvv0[ii] = 162 * c12 + 81 * c10 - 405. / 2 * c11 - 81. / 2 * c13
//									- 81. / 2 * c20
//									+ 405. / 4 * c21
//									- 81 * c22
//									+ 81. / 4 * c23 + 9 * c30
//									- 45 / 2 * c31
//									+ 18 * c32
//									- 9. / 2 * c33
//									+ 495. / 4 * c01
//									- 99 * c02
//									+ 99. / 4 * c03
//									- 99. / 2 * c00;
//			Cvvv0[ii] = 9. / 2 * (3 * (c01 - c02) + c03 - c00);
//
//			Cuuuv0[ii] = -81. / 2 * c01 + 81. / 4 * c02 - 9. / 2 * c03
//										- 297. / 4 * c10
//										+ 243. / 2 * c11
//										- 243. / 4 * c12
//										+ 27. / 2 * c13
//										- 243. / 2 * c21
//										+ 243. / 4 * c22
//										- 27. / 2 * c23
//										+ 81. / 2 * c31
//										- 81. / 4 * c32
//										+ 9. / 2 * c33 + 297 / 4 * c20
//										- 99. / 4 * c30
//										+ 99. / 4 * c00;
//			Cuuvv0[ii] = -405 * c12 - 405. / 2 * c10 + 2025. / 4 * c11
//										+ 405. / 4 * c13 + 162 * c20
//										- 405 * c21
//										+ 324 * c22
//										- 81 * c23 - 81. / 2 * c30
//										+ 405. / 4 * c31
//										- 81 * c32
//										+ 81. / 4 * c33
//										- 405. / 2 * c01
//										+ 162 * c02
//										- 81. / 2 * c03
//										+ 81 * c00;
//			Cuvvv0[ii] = -243. / 2 * c12 + 81. / 2 * c13 + 81. / 4 * c20
//					- 243. / 4 * c21
//										+ 243. / 4 * c22
//										- 81. / 4 * c23 - 81. / 2 * c10
//										+ 243. / 2 * c11
//										- 9. / 2 * c30
//										+ 27. / 2 * c31
//										- 27. / 2 * c32
//										+ 9. / 2 * c33
//										- 297. / 4 * c01
//										+ 297. / 4 * c02
//										- 99. / 4 * c03
//										+ 99. / 4 * c00;
//
//			Cuuuvv0[ii] = 243 * c12 + 243 / 2 * c10 - 1215 / 4 * c11 - 243 / 4 * c13
//										- 243 / 2 * c20
//										+ 1215 / 4 * c21
//										- (243 * c22)
//										+ 243 / 4 * c23 + 81 / 2 * c30
//										- 405 / 4 * c31
//										+ 81 * c32
//										- 81 / 4 * c33
//										+ 405 / 4 * c01
//										- 81 * c02
//										+ 81 / 4 * c03
//										- 81 / 2 * c00;
//			Cuuvvv0[ii] = 1215 / 4 * c12 - 405 / 4 * c13 - 81 * c20 + (243 * c21)
//					- (243 * c22)
//										+ (81 * c23) + 405 / 4 * c10
//										- 1215 / 4 * c11
//										+ 81 / 4 * c30
//										- 243 / 4 * c31
//										+ 243 / 4 * c32
//										- 81 / 4 * c33
//										+ 243 / 2 * c01
//										- 243 / 2 * c02
//										+ 81 / 2 * c03
//										- 81 / 2 * c00;
//
//			Cuuuvvv0[ii] = -729 / 4 * c12 + 243 / 4 * c13 + 243 / 4 * c20
//					- 729 / 4 * c21
//											+ 729 / 4 * c22
//											- 243 / 4 * c23 - 243 / 4 * c10
//											+ 729 / 4 * c11
//											- 81 / 4 * c30
//											+ 243 / 4 * c31
//											- 243 / 4 * c32
//											+ 81 / 4 * c33
//											- 243 / 4 * c01
//											+ 243 / 4 * c02
//											- 81 / 4 * c03
//											+ 81 / 4 * c00;
//		}
//
//		{
//			// Upper interior layer
//			double& c00 = m_nodalValues[17][ii];
//			double& c10 = m_nodalValues[39][ii];
//			double& c20 = m_nodalValues[38][ii];
//			double& c30 = m_nodalValues[19][ii];
//			double& c01 = m_nodalValues[50][ii];
//			double& c11 = m_nodalValues[60][ii];
//			double& c21 = m_nodalValues[61][ii];
//			double& c31 = m_nodalValues[43][ii];
//			double& c02 = m_nodalValues[51][ii];
//			double& c12 = m_nodalValues[63][ii];
//			double& c22 = m_nodalValues[62][ii];
//			double& c32 = m_nodalValues[42][ii];
//			double& c03 = m_nodalValues[23][ii];
//			double& c13 = m_nodalValues[46][ii];
//			double& c23 = m_nodalValues[47][ii];
//			double& c33 = m_nodalValues[21][ii];
//
//			C0[ii] = c00;
//
//			Cu0[ii] = c30 - 11 / 2 * c00 + 9 * c10 - 9 / 2 * c20;
//			Cv0[ii] = 9 * c01 - 9 / 2 * c02 + c03 - 11 / 2 * c00;
//
//			Cuu0[ii] = 4.5 * (-c30 + 2 * c00 - 5 * c10 + 4 * c20);
//			Cuv0[ii] = 99. / 4 * (-2 * c01 - 2 * c10 + c02 + c20)
//					+ 81. / 4 * (4 * c11 - 2 * c12 - 2 * c21 + c22)
//					+ 9 * (c13 + c31 - 0.5 * (c23 + c32)) + c33
//									- 5.5 * (c03 + c30)
//									+ 121. / 4 * c00;
//			Cvv0[ii] = 4.5 * (-5 * c01 + 4 * c02 - c03 + 2 * c00);
//
//			Cuuu0[ii] = 4.5 * (c30 - c00 + 3 * (c10 - c20));
//			Cuuv0[ii] = 81 * c01 - 81 * c22 - 81. / 2 * (c02 + c31) + 81. / 4 * c32
//									+ 18 * c23 + 9 * c03
//									- 9. / 2 * c33
//									+ 495. / 4 * c10
//									- 405. / 2 * c11
//									+ 405. / 4 * c12
//									- 45. / 2 * c13
//									+ 162 * c21
//									- 99 * c20
//									+ 99. / 4 * c30
//									- 99. / 2 * c00;
//			Cuvv0[ii] = 162 * c12 + 81 * c10 - 405. / 2 * c11 - 81. / 2 * c13
//									- 81. / 2 * c20
//									+ 405. / 4 * c21
//									- 81 * c22
//									+ 81. / 4 * c23 + 9 * c30
//									- 45 / 2 * c31
//									+ 18 * c32
//									- 9. / 2 * c33
//									+ 495. / 4 * c01
//									- 99 * c02
//									+ 99. / 4 * c03
//									- 99. / 2 * c00;
//			Cvvv0[ii] = 9. / 2 * (3 * (c01 - c02) + c03 - c00);
//
//			Cuuuv0[ii] = -81. / 2 * c01 + 81. / 4 * c02 - 9. / 2 * c03
//										- 297. / 4 * c10
//										+ 243. / 2 * c11
//										- 243. / 4 * c12
//										+ 27. / 2 * c13
//										- 243. / 2 * c21
//										+ 243. / 4 * c22
//										- 27. / 2 * c23
//										+ 81. / 2 * c31
//										- 81. / 4 * c32
//										+ 9. / 2 * c33 + 297 / 4 * c20
//										- 99. / 4 * c30
//										+ 99. / 4 * c00;
//			Cuuvv0[ii] = -405 * c12 - 405. / 2 * c10 + 2025. / 4 * c11
//										+ 405. / 4 * c13 + 162 * c20
//										- 405 * c21
//										+ 324 * c22
//										- 81 * c23 - 81. / 2 * c30
//										+ 405. / 4 * c31
//										- 81 * c32
//										+ 81. / 4 * c33
//										- 405. / 2 * c01
//										+ 162 * c02
//										- 81. / 2 * c03
//										+ 81 * c00;
//			Cuvvv0[ii] = -243. / 2 * c12 + 81. / 2 * c13 + 81. / 4 * c20
//					- 243. / 4 * c21
//										+ 243. / 4 * c22
//										- 81. / 4 * c23 - 81. / 2 * c10
//										+ 243. / 2 * c11
//										- 9. / 2 * c30
//										+ 27. / 2 * c31
//										- 27. / 2 * c32
//										+ 9. / 2 * c33
//										- 297. / 4 * c01
//										+ 297. / 4 * c02
//										- 99. / 4 * c03
//										+ 99. / 4 * c00;
//
//			Cuuuvv0[ii] = 243 * c12 + 243 / 2 * c10 - 1215 / 4 * c11 - 243 / 4 * c13
//										- 243 / 2 * c20
//										+ 1215 / 4 * c21
//										- (243 * c22)
//										+ 243 / 4 * c23 + 81 / 2 * c30
//										- 405 / 4 * c31
//										+ 81 * c32
//										- 81 / 4 * c33
//										+ 405 / 4 * c01
//										- 81 * c02
//										+ 81 / 4 * c03
//										- 81 / 2 * c00;
//			Cuuvvv0[ii] = 1215 / 4 * c12 - 405 / 4 * c13 - 81 * c20 + (243 * c21)
//					- (243 * c22)
//										+ (81 * c23) + 405 / 4 * c10
//										- 1215 / 4 * c11
//										+ 81 / 4 * c30
//										- 243 / 4 * c31
//										+ 243 / 4 * c32
//										- 81 / 4 * c33
//										+ 243 / 2 * c01
//										- 243 / 2 * c02
//										+ 81 / 2 * c03
//										- 81 / 2 * c00;
//
//			Cuuuvvv0[ii] = -729 / 4 * c12 + 243 / 4 * c13 + 243 / 4 * c20
//					- 729 / 4 * c21
//											+ 729 / 4 * c22
//											- 243 / 4 * c23 - 243 / 4 * c10
//											+ 729 / 4 * c11
//											- 81 / 4 * c30
//											+ 243 / 4 * c31
//											- 243 / 4 * c32
//											+ 81 / 4 * c33
//											- 243 / 4 * c01
//											+ 243 / 4 * c02
//											- 81 / 4 * c03
//											+ 81 / 4 * c00;
//		}
//
//		{
//			// Top
//			double& c00 = m_nodalValues[4][ii];
//			double& c10 = m_nodalValues[24][ii];
//			double& c20 = m_nodalValues[25][ii];
//			double& c30 = m_nodalValues[5][ii];
//			double& c01 = m_nodalValues[31][ii];
//			double& c11 = m_nodalValues[52][ii];
//			double& c21 = m_nodalValues[53][ii];
//			double& c31 = m_nodalValues[26][ii];
//			double& c02 = m_nodalValues[30][ii];
//			double& c12 = m_nodalValues[55][ii];
//			double& c22 = m_nodalValues[54][ii];
//			double& c32 = m_nodalValues[27][ii];
//			double& c03 = m_nodalValues[7][ii];
//			double& c13 = m_nodalValues[29][ii];
//			double& c23 = m_nodalValues[28][ii];
//			double& c33 = m_nodalValues[6][ii];
//
//			C0[ii] = c00;
//
//			Cu0[ii] = c30 - 11 / 2 * c00 + 9 * c10 - 9 / 2 * c20;
//			Cv0[ii] = 9 * c01 - 9 / 2 * c02 + c03 - 11 / 2 * c00;
//
//			Cuu0[ii] = 4.5 * (-c30 + 2 * c00 - 5 * c10 + 4 * c20);
//			Cuv0[ii] = 99. / 4 * (-2 * c01 - 2 * c10 + c02 + c20)
//					+ 81. / 4 * (4 * c11 - 2 * c12 - 2 * c21 + c22)
//					+ 9 * (c13 + c31 - 0.5 * (c23 + c32)) + c33
//									- 5.5 * (c03 + c30)
//									+ 121. / 4 * c00;
//			Cvv0[ii] = 4.5 * (-5 * c01 + 4 * c02 - c03 + 2 * c00);
//
//			Cuuu0[ii] = 4.5 * (c30 - c00 + 3 * (c10 - c20));
//			Cuuv0[ii] = 81 * c01 - 81 * c22 - 81. / 2 * (c02 + c31) + 81. / 4 * c32
//									+ 18 * c23 + 9 * c03
//									- 9. / 2 * c33
//									+ 495. / 4 * c10
//									- 405. / 2 * c11
//									+ 405. / 4 * c12
//									- 45. / 2 * c13
//									+ 162 * c21
//									- 99 * c20
//									+ 99. / 4 * c30
//									- 99. / 2 * c00;
//			Cuvv0[ii] = 162 * c12 + 81 * c10 - 405. / 2 * c11 - 81. / 2 * c13
//									- 81. / 2 * c20
//									+ 405. / 4 * c21
//									- 81 * c22
//									+ 81. / 4 * c23 + 9 * c30
//									- 45 / 2 * c31
//									+ 18 * c32
//									- 9. / 2 * c33
//									+ 495. / 4 * c01
//									- 99 * c02
//									+ 99. / 4 * c03
//									- 99. / 2 * c00;
//			Cvvv0[ii] = 9. / 2 * (3 * (c01 - c02) + c03 - c00);
//
//			Cuuuv0[ii] = -81. / 2 * c01 + 81. / 4 * c02 - 9. / 2 * c03
//										- 297. / 4 * c10
//										+ 243. / 2 * c11
//										- 243. / 4 * c12
//										+ 27. / 2 * c13
//										- 243. / 2 * c21
//										+ 243. / 4 * c22
//										- 27. / 2 * c23
//										+ 81. / 2 * c31
//										- 81. / 4 * c32
//										+ 9. / 2 * c33 + 297 / 4 * c20
//										- 99. / 4 * c30
//										+ 99. / 4 * c00;
//			Cuuvv0[ii] = -405 * c12 - 405. / 2 * c10 + 2025. / 4 * c11
//										+ 405. / 4 * c13 + 162 * c20
//										- 405 * c21
//										+ 324 * c22
//										- 81 * c23 - 81. / 2 * c30
//										+ 405. / 4 * c31
//										- 81 * c32
//										+ 81. / 4 * c33
//										- 405. / 2 * c01
//										+ 162 * c02
//										- 81. / 2 * c03
//										+ 81 * c00;
//			Cuvvv0[ii] = -243. / 2 * c12 + 81. / 2 * c13 + 81. / 4 * c20
//					- 243. / 4 * c21
//										+ 243. / 4 * c22
//										- 81. / 4 * c23 - 81. / 2 * c10
//										+ 243. / 2 * c11
//										- 9. / 2 * c30
//										+ 27. / 2 * c31
//										- 27. / 2 * c32
//										+ 9. / 2 * c33
//										- 297. / 4 * c01
//										+ 297. / 4 * c02
//										- 99. / 4 * c03
//										+ 99. / 4 * c00;
//
//			Cuuuvv0[ii] = 243 * c12 + 243 / 2 * c10 - 1215 / 4 * c11 - 243 / 4 * c13
//										- 243 / 2 * c20
//										+ 1215 / 4 * c21
//										- (243 * c22)
//										+ 243 / 4 * c23 + 81 / 2 * c30
//										- 405 / 4 * c31
//										+ 81 * c32
//										- 81 / 4 * c33
//										+ 405 / 4 * c01
//										- 81 * c02
//										+ 81 / 4 * c03
//										- 81 / 2 * c00;
//			Cuuvvv0[ii] = 1215 / 4 * c12 - 405 / 4 * c13 - 81 * c20 + (243 * c21)
//					- (243 * c22)
//										+ (81 * c23) + 405 / 4 * c10
//										- 1215 / 4 * c11
//										+ 81 / 4 * c30
//										- 243 / 4 * c31
//										+ 243 / 4 * c32
//										- 81 / 4 * c33
//										+ 243 / 2 * c01
//										- 243 / 2 * c02
//										+ 81 / 2 * c03
//										- 81 / 2 * c00;
//
//			Cuuuvvv0[ii] = -729 / 4 * c12 + 243 / 4 * c13 + 243 / 4 * c20
//					- 729 / 4 * c21
//											+ 729 / 4 * c22
//											- 243 / 4 * c23 - 243 / 4 * c10
//											+ 729 / 4 * c11
//											- 81 / 4 * c30
//											+ 243 / 4 * c31
//											- 243 / 4 * c32
//											+ 81 / 4 * c33
//											- 243 / 4 * c01
//											+ 243 / 4 * c02
//											- 81 / 4 * c03
//											+ 81 / 4 * c00;
//		}
//	}
//}
//

void LagrangeCubicHexMapping::setModalValues() {
	for (int ii = 0; ii < 3; ii++) {
		// Bottom
		double& c000 = m_nodalValues[0][ii];
		double& c100 = m_nodalValues[8][ii];
		double& c200 = m_nodalValues[9][ii];
		double& c300 = m_nodalValues[1][ii];
		double& c010 = m_nodalValues[15][ii];
		double& c110 = m_nodalValues[32][ii];
		double& c210 = m_nodalValues[33][ii];
		double& c310 = m_nodalValues[10][ii];
		double& c020 = m_nodalValues[14][ii];
		double& c120 = m_nodalValues[35][ii];
		double& c220 = m_nodalValues[34][ii];
		double& c320 = m_nodalValues[11][ii];
		double& c030 = m_nodalValues[3][ii];
		double& c130 = m_nodalValues[13][ii];
		double& c230 = m_nodalValues[12][ii];
		double& c330 = m_nodalValues[2][ii];
			// Lower interior layer
		double& c001 = m_nodalValues[16][ii];
		double& c101 = m_nodalValues[36][ii];
		double& c201 = m_nodalValues[37][ii];
		double& c301 = m_nodalValues[18][ii];
		double& c011 = m_nodalValues[49][ii];
		double& c111 = m_nodalValues[56][ii];
		double& c211 = m_nodalValues[57][ii];
		double& c311 = m_nodalValues[40][ii];
		double& c021 = m_nodalValues[48][ii];
		double& c121 = m_nodalValues[59][ii];
		double& c221 = m_nodalValues[58][ii];
		double& c321 = m_nodalValues[41][ii];
		double& c031 = m_nodalValues[22][ii];
		double& c131 = m_nodalValues[45][ii];
		double& c231 = m_nodalValues[44][ii];
		double& c331 = m_nodalValues[20][ii];

			// Upper interior layer
		double& c002 = m_nodalValues[17][ii];
		double& c102 = m_nodalValues[39][ii];
		double& c202 = m_nodalValues[38][ii];
		double& c302 = m_nodalValues[19][ii];
		double& c012 = m_nodalValues[50][ii];
		double& c112 = m_nodalValues[60][ii];
		double& c212 = m_nodalValues[61][ii];
		double& c312 = m_nodalValues[43][ii];
		double& c022 = m_nodalValues[51][ii];
		double& c122 = m_nodalValues[63][ii];
		double& c222 = m_nodalValues[62][ii];
		double& c322 = m_nodalValues[42][ii];
		double& c032 = m_nodalValues[23][ii];
		double& c132 = m_nodalValues[46][ii];
		double& c232 = m_nodalValues[47][ii];
		double& c332 = m_nodalValues[21][ii];

			// Top
		double& c003 = m_nodalValues[4][ii];
		double& c103 = m_nodalValues[24][ii];
		double& c203 = m_nodalValues[25][ii];
		double& c303 = m_nodalValues[5][ii];
		double& c013 = m_nodalValues[31][ii];
		double& c113 = m_nodalValues[52][ii];
		double& c213 = m_nodalValues[53][ii];
		double& c313 = m_nodalValues[26][ii];
		double& c023 = m_nodalValues[30][ii];
		double& c123 = m_nodalValues[55][ii];
		double& c223 = m_nodalValues[54][ii];
		double& c323 = m_nodalValues[27][ii];
		double& c033 = m_nodalValues[7][ii];
		double& c133 = m_nodalValues[29][ii];
		double& c233 = m_nodalValues[28][ii];
		double& c333 = m_nodalValues[6][ii];
		// Constant term

		C[ii] = c000;

		// Linear terms (3)

		Cu[ii] = +(-11. / 2. * c000 + 9. * c100 + c300 - 9. / 2. * c200);
		Cv[ii] = +(-11. / 2. * c000 + c030 + 9. * c010 - 9. / 2. * c020);
		Cw[ii] = +(-11. / 2. * c000 + 9. * c001 - 9. / 2. * c002 + c003);

		// Quadratic terms (6) 2-0-0 (3), 1-1-0 (3)

		Cu2[ii] = +(-9. / 2. * c300 + 18. * c200 - 45. / 2. * c100 + 9. * c000);
		Cv2[ii] = +(-9. / 2. * c030 + 18. * c020 - 45. / 2. * c010 + 9. * c000);
		Cw2[ii] = +(18. * c002 - 9. / 2. * c003 - 45. / 2. * c001 + 9. * c000);

		Cuv[ii] = +(121. / 4. * c000 + 99. / 4. * c020 - 99. / 2. * c010
						- 11. / 2. * c030 - 9. / 2. * c320
						+ c330 + 9. * c310
						- 81. / 2. * c120 - 81. / 2. * c210
						+ 81. / 4. * c220
						- 9. / 2. * c230
						+ 9. * c130 + 81. * c110
						- 99. / 2. * c100 - 11. / 2. * c300
								+ 99. / 4. * c200);

		Cuw[ii] = +(121. / 4. * c000 - 11. / 2. * c003 + 99. / 4. * c002
				- 99. / 2. * c001
						+ 81. * c101 + c303
						- 9. / 2. * c203 - 81. / 2. * c102
						+ 81. / 4. * c202 + 9. * c301
						- 9. / 2. * c302 - 81. / 2. * c201
						+ 9. * c103
						- 99. / 2. * c100 - 11. / 2. * c300
								+ 99. / 4. * c200);

		Cvw[ii] = +(121. / 4. * c000 - 11. / 2. * c003 + 99. / 4. * c002
				- 99. / 2. * c001
						+ 99. / 4. * c020
						- 99. / 2. * c010 - 11. / 2. * c030
						+ 81. / 4. * c022
						- 9. / 2. * c023
						+ 9. * c031
						- 9. / 2. * c032
						+ c033 + 81. * c011
						- 81. / 2. * c012
						+ 9. * c013
								- 81. / 2. * c021);

		// Cubic terms (10) 3-0-0 (3), 2-1-0 (6), 1-1-1 (1)

		Cu3[ii] = +(9. / 2. * c300 - 27. / 2. * c200 + 27. / 2. * c100
				- 9. / 2. * c000);

		Cv3[ii] = +(9. / 2. * c030 - 27. / 2. * c020 + 27. / 2. * c010
				- 9. / 2. * c000);

		Cw3[ii] = +(-27. / 2. * c002 + 9. / 2. * c003 + 27. / 2. * c001
				- 9. / 2. * c000);

		Cu2v[ii] = +(-9. / 2. * c330 - 81. / 2. * c310 + 81. / 4. * c320
				- 99. / 2. * c000 - 81. / 2. * c020
							+ 81. * c010 + 9. * c030 + 162. * c210 + 495. / 4. * c100
							- 81. * c220
							+ 18. * c230 + 99. / 4. * c300
							- 99. * c200 - 45. / 2. * c130 - 405. / 2. * c110
									+ 405. / 4. * c120);

		Cuv2[ii] = +(-9. / 2. * c330 - 45. / 2. * c310 + 18. * c320
				- 99. / 2. * c000
							- 99. * c020
							+ 495. / 4. * c010
							+ 99. / 4. * c030
							+ 405. / 4. * c210
							+ 81. * c100
							- 81. * c220
							+ 81. / 4. * c230 + 9. * c300
							- 81. / 2. * c200 - 81. / 2. * c130 - 405. / 2. * c110
									+ 162. * c120);

		Cu2w[ii] = +(18. * c203 - 9. / 2. * c303 - 405. / 2. * c101
									- 99. / 2. * c000
							- 81. / 2. * c002
							+ 9. * c003 + 81. * c001 + 495. / 4. * c100 + 405. / 4. * c102
							- 81. * c202
							+ 99. / 4. * c300
							- 81. / 2. * c301
							+ 81. / 4. * c302 + 162. * c201
									- 99. * c200 - 45. / 2. * c103);

		Cuw2[ii] = +(81. / 4. * c203 - 9. / 2. * c303 - 405. / 2. * c101
							- 99. / 2. * c000 - 99. * c002
							+ 99. / 4. * c003 + 495. / 4. * c001 + 81. * c100 + (162 * c102)
							- 81. * c202
							+ 9. * c300
							- 45. / 2. * c301
							+ (18 * c302) + 405. / 4. * c201
									- 81. / 2. * c200 - 81. / 2. * c103);

		Cv2w[ii] = +(-99. / 2. * c000 - 81. / 2. * c002 + 9. * c003 + 81. * c001
							+ 162. * c021
							- 81. * c022 - 99. * c020
							+ 405. / 4. * c012
							- 45. / 2. * c013 - 405. / 2. * c011
							+ 495. / 4. * c010 + 81. / 4. * c032
							- 9. / 2. * c033
							+ 99. / 4. * c030
							- 81. / 2. * c031
									+ 18. * c023);

		Cvw2[ii] = +(-99. / 2. * c000 + 99. / 4. * c003 - 99. * c002
				+ 495. / 4. * c001
							+ 405. / 4. * c021
							- 81. * c022 - 81. / 2. * c020
							+ 162. * c012
							- 81. / 2. * c013 - 405. / 2. * c011
							+ 81. * c010 + 18. * c032
							- 9. / 2. * c033
							+ 9. * c030
							- 45. / 2. * c031
									+ 81. / 4. * c023);

		Cuvw[ii] = +(-81. / 2. * c231 + 99. / 4. * c230 + 81. * c311
									+ 1089. / 4. * c010
							+ 81. / 4. * c232
							- 1331. / 8. * c000 - 729. / 2. * c112
							+ 729. / 4. * c221
							- 81. / 2. * c123
							+ 891. / 4. * c102
							+ 9. * c331
							+ 891. / 4. * c120
							+ 891. / 4. * c021
							- 81. / 2. * c321
							- 1089. / 8. * c020
							- 99. / 2. * c103
							- 99. / 2. * c310
							+ 729. * c111
							- 891. / 8. * c202 - 99. / 2. * c013 - 9. / 2. * c323
							+ 729. / 4. * c122
							- 81. / 2. * c312
							+ 9. * c313 + 99. / 4. * c023
							- 729. / 2. * c121 - 9. / 2. * c233
							+ 121. / 4. * c300
							- 9. / 2. * c332
							+ 99. / 4. * c302 + 81. * c131
							- 99. / 2. * c301
							+ c333 + 81. * c113
							- 729. / 2. * c211
							+ 81. / 4. * c322 + 9. * c133
							- 729. / 8. * c222
							+ 81. / 4. * c223 + 99. / 4. * c203
							- 1089. / 8. * c200
							+ 891. / 4. * c201
							- 11. / 2. * c033 - 99. / 2. * c031
							+ 99. / 4. * c032
							- 99. / 2. * c130
							+ 1089. / 4. * c100
							- 81. / 2. * c213
							+ 729. / 4. * c212
							- 891. / 8. * c022 - 11. / 2. * c303
							+ 121. / 4. * c030 + 121. / 4. * c003 + 891. / 4. * c210
							- 1089. / 8. * c002
							+ 1089. / 4. * c001 + 891. / 4. * c012
							- 81. / 2. * c132
							+ 99. / 4. * c320
							- 11. / 2. * c330
							- 891. / 2. * c110
							- 891. / 2. * c011
									- 891. / 8. * c220 - 891. / 2. * c101);

		// Quartic terms (12?) 2-2-0 (3), 2-1-1 (3), 3-1-0 (6)

		Cu2v2[ii] = +(81. / 4. * c330 + 405. / 4. * c310 - 81. * c320 + 81. * c000
							+ 162. * c020
							- 405. / 2. * c010
							- 81. / 2. * c030
							- 405. * c210
							- 405. / 2. * c100
							+ 324. * c220
							- 81. * c230 - 81. / 2. * c300
							+ 162. * c200 + 405. / 4. * c130 + 2025. / 4. * c110
									- 405. * c120);

		Cu2w2[ii] = +(-81. * c203 + 81. / 4. * c303 + 2025. / 4. * c101 + 81. * c000
				- 81. / 2. * c003
							+ 162. * c002
							- 405. / 2. * c001 - 405. / 2. * c100 - (405 * c102)
							+ 324. * c202
							- 81. / 2. * c300
							+ 405. / 4. * c301
							- (81 * c302) - 405. * c201
									+ 162. * c200 + 405. / 4. * c103);

		Cv2w2[ii] = +(81. * c000 - 81. / 2. * c003 + 162. * c002 - 405. / 2. * c001
							- 405. * c021
							+ 324. * c022 + 162. * c020
							- 405. * c012
							+ 405. / 4. * c013 + 2025. / 4. * c011
							- 405. / 2. * c010 - 81. * c032
							+ 81. / 4. * c033
							- 81. / 2. * c030
							+ 405. / 4. * c031
									- 81. * c023);
		;

		Cu2vw[ii] = +(162. * c231 - 99. * c230 - 729. / 2. * c311 - 891. / 2. * c010
							- 81. * c232
							+ 1089. / 4. * c000 + 3645. / 4. * c112
							- 729. * c221
							+ 405. / 4. * c123
							- 0.4455e4 / 8. * c102
							- 81. / 2. * c331
							- 0.4455e4 / 8. * c120
							- 729. / 2. * c021
							+ 729. / 4. * c321
							+ 891. / 4. * c020
							+ 495. / 4. * c103
							+ 891. / 4. * c310
							- 3645. / 2. * c111
							+ 891. / 2. * c202 + 81. * c013 + 81. / 4. * c323
							- 3645. / 8. * c122
							+ 729. / 4. * c312
							- 81. / 2. * c313 - 81. / 2. * c023
							+ 3645. / 4. * c121 + 18. * c233
							- 1089. / 8. * c300
							+ 81. / 4. * c332
							- 891. / 8. * c302 - 405. / 2. * c131
							+ 891. / 4. * c301
							- 9. / 2. * c333 - 405. / 2. * c113
							+ 0.1458e4 * c211
							- 729. / 8. * c322 - 45. / 2. * c133
							+ 729. / 2. * c222
							- 81. * c223 - 99. * c203
							+ 1089. / 2. * c200
							- 891. * c201
							+ 9. * c033 + 81. * c031
							- 81. / 2. * c032
							+ 495. / 4. * c130
							- 0.5445e4 / 8. * c100
							+ 162. * c213
							- 729. * c212
							+ 729. / 4. * c022 + 99. / 4. * c303
							- 99. / 2. * c030 - 99. / 2. * c003 - 891. * c210
							+ 891. / 4. * c002
							- 891. / 2. * c001 - 729. / 2. * c012
							+ 405. / 4. * c132
							- 891. / 8. * c320
							+ 99. / 4. * c330
							+ 0.4455e4 / 4. * c110
							+ 729. * c011
									+ 891. / 2. * c220 + 0.4455e4 / 4. * c101);

		Cuv2w[ii] = +(729. / 4. * c231 - 891. / 8. * c230 - 405. / 2. * c311
							- 0.5445e4 / 8. * c010 - 729. / 8. * c232
							+ 1089. / 4. * c000 + 3645. / 4. * c112
							- 729. * c221
							+ 162. * c123
							- 729. / 2. * c102 - 81. / 2. * c331 - 891. * c120 - 891. * c021
							+ 162. * c321 + 1089. / 2. * c020 + 81. * c103 + 495. / 4. * c310
							- 3645. / 2. * c111
							+ 729. / 4. * c202 + 495. / 4. * c013 + 18. * c323
							- (729 * c122)
							+ 405. / 4. * c312
							- 0.45e2 / 2. * c313 - 99. * c023
							+ 0.1458e4 * c121 + 81. / 4. * c233
							- 99. / 2. * c300
							+ 81. / 4. * c332
							- 81. / 2. * c302 - 729. / 2. * c131
							+ 81. * c301
							- 9. / 2. * c333 - 405. / 2. * c113
							+ 3645. / 4. * c211
							- (81 * c322) - 81. / 2. * c133
							+ 729. / 2. * c222
							- 81. * c223 - 81. / 2. * c203
							+ 891. / 4. * c200
							- 729. / 2. * c201
							+ 99. / 4. * c033 + 891. / 4. * c031
							- 891. / 8. * c032
							+ 891. / 4. * c130
							- 891. / 2. * c100
							+ 405. / 4. * c213
							- 3645. / 8. * c212
							+ 891. / 2. * c022 + 9. * c303
							- 1089. / 8. * c030 - 99. / 2. * c003 - 0.4455e4 / 8. * c210
							+ 891. / 4. * c002
							- 891. / 2. * c001 - 0.4455e4 / 8. * c012
							+ 729. / 4. * c132
							- 99. * c320
							+ 99. / 4. * c330
							+ 0.4455e4 / 4. * c110
							+ 0.4455e4 / 4. * c011
							+ 891. / 2. * c220 + 729. * c101);

		Cuvw2[ii] = +(405. / 4. * c231 - 81. / 2. * c230 - 405. / 2. * c311
							- 891. / 2. * c010 - 81. * c232
							+ 1089. / 4. * c000 + (1458 * c112)
							- 3645. / 8. * c221
							+ 729. / 4. * c123
							- (891 * c102)
							- 0.45e2 / 2. * c331
							- 729. / 2. * c120
							- 0.4455e4 / 8. * c021
							+ 405. / 4. * c321
							+ 891. / 4. * c020
							+ 891. / 4. * c103
							+ 81. * c310
							- 3645. / 2. * c111
							+ 891. / 2. * c202 + 891. / 4. * c013 + 81. / 4. * c323
							- (729 * c122)
							+ (162 * c312)
							- 81. / 2. * c313 - 891. / 8. * c023
							+ 3645. / 4. * c121 + 81. / 4. * c233
							- 99. / 2. * c300
							+ (18 * c332)
							- (99 * c302) - 405. / 2. * c131
							+ 495. / 4. * c301
							- 9. / 2. * c333 - 729. / 2. * c113
							+ 3645. / 4. * c211
							- (81 * c322) - 81. / 2. * c133
							+ 729. / 2. * c222
							- 729. / 8. * c223 - 891. / 8. * c203
							+ 891. / 4. * c200
							- 0.4455e4 / 8. * c201
							+ 99. / 4. * c033 + 495. / 4. * c031
							- 99. * c032
							+ 81. * c130
							- 891. / 2. * c100
							+ 729. / 4. * c213
							- 729. * c212
							+ 891. / 2. * c022 + 99. / 4. * c303
							- 99. / 2. * c030 - 1089. / 8. * c003 - 729. / 2. * c210
							+ 1089. / 2. * c002
							- 0.5445e4 / 8. * c001 - 891. * c012
							+ (162 * c132)
							- 81. / 2. * c320
							+ 9. * c330
							+ 729. * c110
							+ 0.4455e4 / 4. * c011
									+ 729. / 4. * c220 + 0.4455e4 / 4. * c101);

		Cu3v[ii] = +(9. / 2. * c330 + 81. / 2. * c310 - 81. / 4. * c320
				+ 99. / 4. * c000 + 81. / 4. * c020
							- 81. / 2. * c010
							- 9. / 2. * c030
							- 0.243e3 / 2. * c210
							- 0.297e3 / 4. * c100
							+ 0.243e3 / 4. * c220
							- 0.27e2 / 2. * c230 - 99. / 4. * c300
							+ 0.297e3 / 4. * c200 + 0.27e2 / 2. * c130 + 0.243e3 / 2. * c110
									- 0.243e3 / 4. * c120);

		Cuv3[ii] = +(9. / 2. * c330 + 0.27e2 / 2. * c310 - 0.27e2 / 2. * c320
				+ 99. / 4. * c000 + 0.297e3 / 4. * c020
							- 0.297e3 / 4. * c010
							- 99. / 4. * c030
							- 0.243e3 / 4. * c210
							- 81. / 2. * c100
							+ 0.243e3 / 4. * c220
							- 81. / 4. * c230 - 9. / 2. * c300
							+ 81. / 4. * c200 + 81. / 2. * c130 + 0.243e3 / 2. * c110
									- 0.243e3 / 2. * c120);

		Cv3w[ii] = +(99. / 4. * c000 + 81. / 4. * c002 - 9. / 2. * c003
							- 81. / 2. * c001 - 0.243e3 / 2. * c021
							+ 0.243e3 / 4. * c022
							+ 0.297e3 / 4. * c020
							+ 0.27e2 / 2. * c013
							+ 0.243e3 / 2. * c011
							- 0.243e3 / 4. * c012 - 0.297e3 / 4. * c010 - 81. / 4. * c032
							+ 9. / 2. * c033
							- 99. / 4. * c030
							+ 81. / 2. * c031
									- 0.27e2 / 2. * c023);

		Cvw3[ii] = +(99. / 4. * c000 - 99. / 4. * c003 + 0.297e3 / 4. * c002
				- 0.297e3 / 4. * c001 - 0.243e3 / 4. * c021
							+ 0.243e3 / 4. * c022
							+ 81. / 4. * c020
							+ 81. / 2. * c013
							+ 0.243e3 / 2. * c011
							- 0.243e3 / 2. * c012 - 81. / 2. * c010 - 0.27e2 / 2. * c032
							+ 9. / 2. * c033
							- 9. / 2. * c030
							+ 0.27e2 / 2. * c031
									- 81. / 4. * c023);

		Cuw3[ii] = +(-81. / 4. * c203 + 9. / 2. * c303 + 0.243e3 / 2. * c101
							+ 99. / 4. * c000
							- 99. / 4. * c003
							+ 0.297e3 / 4. * c002
							- 0.297e3 / 4. * c001 - 81. / 2. * c100 - 0.243e3 / 2. * c102
							+ 0.243e3 / 4. * c202
							- 9. / 2. * c300
							+ 0.27e2 / 2. * c301
							- 0.27e2 / 2. * c302 - 0.243e3 / 4. * c201
									+ 81. / 4. * c200 + 81. / 2. * c103);

		Cu3w[ii] = +(-0.27e2 / 2. * c203 + 9. / 2. * c303 + 0.243e3 / 2. * c101
							+ 99. / 4. * c000 + 81. / 4. * c002
							- 9. / 2. * c003
							- 81. / 2. * c001
							- 0.297e3 / 4. * c100
							- 0.243e3 / 4. * c102
							+ 0.243e3 / 4. * c202
							- 99. / 4. * c300
							+ 81. / 2. * c301
							- 81. / 4. * c302 - 0.243e3 / 2. * c201
									+ 0.297e3 / 4. * c200 + 0.27e2 / 2. * c103);

		// Quintic terms (12?) 3-1-1 (3), 3-2-0 (6), 2-2-1 (3)

		Cuvw3[ii] = +(-0.243e3 / 4. * c231 + 81. / 4. * c230 + 0.243e3 / 2. * c311
							+ 891. / 4. * c010 + 0.243e3 / 4. * c232
							- 1089. / 8. * c000 - 0.2187e4 / 2. * c112
							+ 0.2187e4 / 8. * c221
							- 729. / 4. * c123
							+ 0.2673e4 / 4. * c102
							+ 0.27e2 / 2. * c331
							+ 729. / 4. * c120
							+ 0.2673e4 / 8. * c021
							- 0.243e3 / 4. * c321
							- 891. / 8. * c020
							- 891. / 4. * c103
							- 81. / 2. * c310
							+ 0.2187e4 / 2. * c111
							- 0.2673e4 / 8. * c202 - 891. / 4. * c013 - 81. / 4. * c323
							+ 0.2187e4 / 4. * c122
							- 0.243e3 / 2. * c312
							+ 81. / 2. * c313 + 891. / 8. * c023
							- 0.2187e4 / 4. * c121 - 81. / 4. * c233
							+ 99. / 4. * c300
							- 0.27e2 / 2. * c332
							+ 0.297e3 / 4. * c302 + 0.243e3 / 2. * c131
							- 0.297e3 / 4. * c301
							+ 9. / 2. * c333 + 729. / 2. * c113
							- 0.2187e4 / 4. * c211
							+ 0.243e3 / 4. * c322 + 81. / 2. * c133
							- 0.2187e4 / 8. * c222
							+ 729. / 8. * c223 + 891. / 8. * c203
							- 891. / 8. * c200
							+ 0.2673e4 / 8. * c201
							- 99. / 4. * c033 - 0.297e3 / 4. * c031
							+ 0.297e3 / 4. * c032
							- 81. / 2. * c130
							+ 891. / 4. * c100
							- 729. / 4. * c213
							+ 0.2187e4 / 4. * c212
							- 0.2673e4 / 8. * c022 - 99. / 4. * c303
							+ 99. / 4. * c030 + 1089. / 8. * c003 + 729. / 4. * c210
							- 0.3267e4 / 8. * c002
							+ 0.3267e4 / 8. * c001 + 0.2673e4 / 4. * c012
							- 0.243e3 / 2. * c132
							+ 81. / 4. * c320
							- 9. / 2. * c330
							- 729. / 2. * c110
							- 0.2673e4 / 4. * c011
									- 729. / 8. * c220 - 0.2673e4 / 4. * c101);

		Cu3vw[ii] = +(-0.243e3 / 2. * c231 + 0.297e3 / 4. * c230 + 729. / 2. * c311
							+ 891. / 4. * c010 + 0.243e3 / 4. * c232
							- 1089. / 8. * c000 - 0.2187e4 / 4. * c112
							+ 0.2187e4 / 4. * c221
							- 0.243e3 / 4. * c123
							+ 0.2673e4 / 8. * c102
							+ 81. / 2. * c331
							+ 0.2673e4 / 8. * c120
							+ 729. / 4. * c021
							- 729. / 4. * c321
							- 891. / 8. * c020
							- 0.297e3 / 4. * c103
							- 891. / 4. * c310
							+ 0.2187e4 / 2. * c111
							- 0.2673e4 / 8. * c202 - 81. / 2. * c013 - 81. / 4. * c323
							+ 0.2187e4 / 8. * c122
							- 729. / 4. * c312
							+ 81. / 2. * c313 + 81. / 4. * c023
							- 0.2187e4 / 4. * c121 - 0.27e2 / 2. * c233
							+ 1089. / 8. * c300
							- 81. / 4. * c332
							+ 891. / 8. * c302 + 0.243e3 / 2. * c131
							- 891. / 4. * c301
							+ 9. / 2. * c333 + 0.243e3 / 2. * c113
							- 0.2187e4 / 2. * c211
							+ 729. / 8. * c322 + 0.27e2 / 2. * c133
							- 0.2187e4 / 8. * c222
							+ 0.243e3 / 4. * c223 + 0.297e3 / 4. * c203
							- 0.3267e4 / 8. * c200
							+ 0.2673e4 / 4. * c201
							- 9. / 2. * c033 - 81. / 2. * c031
							+ 81. / 4. * c032
							- 0.297e3 / 4. * c130
							+ 0.3267e4 / 8. * c100
							- 0.243e3 / 2. * c213
							+ 0.2187e4 / 4. * c212
							- 729. / 8. * c022 - 99. / 4. * c303
							+ 99. / 4. * c030 + 99. / 4. * c003 + 0.2673e4 / 4. * c210
							- 891. / 8. * c002
							+ 891. / 4. * c001 + 729. / 4. * c012
							- 0.243e3 / 4. * c132
							+ 891. / 8. * c320
							- 99. / 4. * c330
							- 0.2673e4 / 4. * c110
							- 729. / 2. * c011
									- 0.2673e4 / 8. * c220 - 0.2673e4 / 4. * c101);

		Cuv3w[ii] = +(-729. / 4. * c231 + 891. / 8. * c230 + 0.243e3 / 2. * c311
							+ 0.3267e4 / 8. * c010 + 729. / 8. * c232
							- 1089. / 8. * c000 - 0.2187e4 / 4. * c112
							+ 0.2187e4 / 4. * c221
							- 0.243e3 / 2. * c123
							+ 729. / 4. * c102
							+ 81. / 2. * c331
							+ 0.2673e4 / 4. * c120
							+ 0.2673e4 / 4. * c021
							- 0.243e3 / 2. * c321
							- 0.3267e4 / 8. * c020
							- 81. / 2. * c103
							- 0.297e3 / 4. * c310
							+ 0.2187e4 / 2. * c111
							- 729. / 8. * c202 - 0.297e3 / 4. * c013 - 0.27e2 / 2. * c323
							+ 0.2187e4 / 4. * c122
							- 0.243e3 / 4. * c312
							+ 0.27e2 / 2. * c313 + 0.297e3 / 4. * c023
							- 0.2187e4 / 2. * c121 - 81. / 4. * c233
							+ 99. / 4. * c300
							- 81. / 4. * c332
							+ 81. / 4. * c302 + 729. / 2. * c131
							- 81. / 2. * c301
							+ 9. / 2. * c333 + 0.243e3 / 2. * c113
							- 0.2187e4 / 4. * c211
							+ 0.243e3 / 4. * c322 + 81. / 2. * c133
							- 0.2187e4 / 8. * c222
							+ 0.243e3 / 4. * c223 + 81. / 4. * c203
							- 891. / 8. * c200
							+ 729. / 4. * c201
							- 99. / 4. * c033 - 891. / 4. * c031
							+ 891. / 8. * c032
							- 891. / 4. * c130
							+ 891. / 4. * c100
							- 0.243e3 / 4. * c213
							+ 0.2187e4 / 8. * c212
							- 0.2673e4 / 8. * c022 - 9. / 2. * c303
							+ 1089. / 8. * c030 + 99. / 4. * c003 + 0.2673e4 / 8. * c210
							- 891. / 8. * c002
							+ 891. / 4. * c001 + 0.2673e4 / 8. * c012
							- 729. / 4. * c132
							+ 0.297e3 / 4. * c320
							- 99. / 4. * c330
							- 0.2673e4 / 4. * c110
							- 0.2673e4 / 4. * c011
									- 0.2673e4 / 8. * c220
									- 729. / 2. * c101);

		Cv3w2[ii] = +(-81. / 2. * c000 + 81. / 4. * c003 - 81. * c002
				+ 405. / 4. * c001
							+ 0.1215e4 / 4. * c021
							- 0.243e3 * c022 - 0.243e3 / 2. * c020
							+ 0.243e3 * c012
							- 0.243e3 / 4. * c013 - 0.1215e4 / 4. * c011
							+ 0.243e3 / 2. * c010 + 81. * c032
							- 81. / 4. * c033
							+ 81. / 2. * c030
							- 405. / 4. * c031
									+ 0.243e3 / 4. * c023);

		Cv2w3[ii] = +(-81. / 2. * c000 + 81. / 2. * c003 - 0.243e3 / 2. * c002
				+ 0.243e3 / 2. * c001 + 0.243e3 * c021
							- 0.243e3 * c022 - 81. * c020
							+ 0.1215e4 / 4. * c012
							- 405. / 4. * c013 - 0.1215e4 / 4. * c011
							+ 405. / 4. * c010 + 0.243e3 / 4. * c032
							- 81. / 4. * c033
							+ 81. / 4. * c030
							- 0.243e3 / 4. * c031
									+ 81. * c023);

		Cu3v2[ii] = +(-81. / 4. * c330 - 405. / 4. * c310 + 81. * c320
				- 81. / 2. * c000
							- 81. * c020
							+ 405. / 4. * c010
							+ 81. / 4. * c030
							+ 0.1215e4 / 4. * c210
							+ 0.243e3 / 2. * c100
							- 0.243e3 * c220
							+ 0.243e3 / 4. * c230 + 81. / 2. * c300
							- 0.243e3 / 2. * c200 - 0.243e3 / 4. * c130 - 0.1215e4 / 4. * c110
									+ 0.243e3 * c120);

		Cu2v3[ii] = +(-81. / 4. * c330 - 0.243e3 / 4. * c310 + 0.243e3 / 4. * c320
				- 81. / 2. * c000 - 0.243e3 / 2. * c020
							+ 0.243e3 / 2. * c010
							+ 81. / 2. * c030
							+ 0.243e3 * c210
							+ 405. / 4. * c100
							- 0.243e3 * c220
							+ 81. * c230 + 81. / 4. * c300
							- 81. * c200 - 405. / 4. * c130 - 0.1215e4 / 4. * c110
									+ 0.1215e4 / 4. * c120);

		Cu2w3[ii] = +(81. * c203 - 81. / 4. * c303 - 0.1215e4 / 4. * c101
							- 81. / 2. * c000
							+ 81. / 2. * c003
							- 0.243e3 / 2. * c002
							+ 0.243e3 / 2. * c001 + 405. / 4. * c100 + 0.1215e4 / 4. * c102
							- 0.243e3 * c202
							+ 81. / 4. * c300
							- 0.243e3 / 4. * c301
							+ 0.243e3 / 4. * c302 + 0.243e3 * c201
									- 81. * c200 - 405. / 4. * c103);

		Cu3w2[ii] = +(0.243e3 / 4. * c203 - 81. / 4. * c303 - 0.1215e4 / 4. * c101
							- 81. / 2. * c000
							+ 81. / 4. * c003
							- 81. * c002
							+ 405. / 4. * c001 + 0.243e3 / 2. * c100 + (243 * c102)
							- 0.243e3 * c202
							+ 81. / 2. * c300
							- 405. / 4. * c301
							+ (81 * c302) + 0.1215e4 / 4. * c201
									- 0.243e3 / 2. * c200 - 0.243e3 / 4. * c103);

		Cu2v2w[ii] = +(-729. * c231 + 891. / 2. * c230 + 3645. / 4. * c311
								+ 0.4455e4 / 4. * c010 + 729. / 2. * c232
								- 891. / 2. * c000 - 0.18225e5 / 8. * c112
								+ 0.2916e4 * c221
								- 405. * c123
								+ 3645. / 4. * c102
								+ 729. / 4. * c331
								+ 0.4455e4 / 2. * c120
								+ 0.1458e4 * c021
								- 729. * c321
								- 891. * c020
								- 405. / 2. * c103
								- 0.4455e4 / 8. * c310
								+ 0.18225e5 / 4. * c111
								- 729. * c202 - 405. / 2. * c013 - 81. * c323
								+ 3645. / 2. * c122
								- 3645. / 8. * c312
								+ 405. / 4. * c313 + 162. * c023
								- 3645. * c121 - 81. * c233
								+ 891. / 4. * c300
								- 729. / 8. * c332
								+ 729. / 4. * c302 + 3645. / 4. * c131
								- 729. / 2. * c301
								+ 81. / 4. * c333 + 2025. / 4. * c113
								- 3645. * c211
								+ 729. / 2. * c322 + 405. / 4. * c133
								- 0.1458e4 * c222
								+ 324. * c223 + 162. * c203
								- 891. * c200
								+ 0.1458e4 * c201
								- 81. / 2. * c033 - 729. / 2. * c031
								+ 729. / 4. * c032
								- 0.4455e4 / 8. * c130
								+ 0.4455e4 / 4. * c100
								- 405. * c213
								+ 3645. / 2. * c212
								- 729. * c022 - 81. / 2. * c303
								+ 891. / 4. * c030 + 81. * c003 + 0.4455e4 / 2. * c210
								- 729. / 2. * c002
								+ 729. * c001 + 3645. / 4. * c012
								- 3645. / 8. * c132
								+ 891. / 2. * c320
								- 891. / 8. * c330
								- 0.22275e5 / 8. * c110
								- 3645. / 2. * c011
								- 0.1782e4 * c220 - 3645. / 2. * c101);

		Cu2vw2[ii] = +(-405. * c231 + 162. * c230 + 3645. / 4. * c311 + 729. * c010
								+ 324. * c232
								- 891. / 2. * c000 - (3645 * c112)
								+ 3645. / 2. * c221
								- 3645. / 8. * c123
								+ 0.4455e4 / 2. * c102
								+ 405. / 4. * c331
								+ 3645. / 4. * c120
								+ 3645. / 4. * c021
								- 3645. / 8. * c321
								- 729. / 2. * c020
								- 0.4455e4 / 8. * c103
								- 729. / 2. * c310
								+ 0.18225e5 / 4. * c111
								- 0.1782e4 * c202 - 729. / 2. * c013 - 729. / 8. * c323
								+ 3645. / 2. * c122
								- (729 * c312)
								+ 729. / 4. * c313 + 729. / 4. * c023
								- 0.18225e5 / 8. * c121 - 81. * c233
								+ 891. / 4. * c300
								- (81 * c332)
								+ 891. / 2. * c302 + 2025. / 4. * c131
								- 0.4455e4 / 8. * c301
								+ 81. / 4. * c333 + 3645. / 4. * c113
								- 3645. * c211
								+ 729. / 2. * c322 + 405. / 4. * c133
								- 0.1458e4 * c222
								+ 729. / 2. * c223 + 891. / 2. * c203
								- 891. * c200
								+ 0.4455e4 / 2. * c201
								- 81. / 2. * c033 - 405. / 2. * c031
								+ 162. * c032
								- 405. / 2. * c130
								+ 0.4455e4 / 4. * c100
								- 729. * c213
								+ 0.2916e4 * c212
								- 729. * c022 - 891. / 8. * c303
								+ 81. * c030 + 891. / 4. * c003 + 0.1458e4 * c210
								- 891. * c002
								+ 0.4455e4 / 4. * c001 + 0.1458e4 * c012
								- (405 * c132)
								+ 729. / 4. * c320
								- 81. / 2. * c330
								- 3645. / 2. * c110
								- 3645. / 2. * c011
										- 729. * c220 - 0.22275e5 / 8. * c101);

		Cuv2w2[ii] =
				+(-3645. / 8. * c231 + 729. / 4. * c230 + 2025. / 4. * c311
								+ 0.4455e4 / 4. * c010 + 729. / 2. * c232
								- 891. / 2. * c000 - (3645 * c112)
								+ 0.3645e4 / 2. * c221
								- 729. * c123
								+ (1458 * c102)
								+ 405. / 4. * c331
								+ 0.1458e4 * c120
								+ 0.4455e4 / 2. * c021
								- 405. * c321
								- 891. * c020
								- 729. / 2. * c103
								- 405. / 2. * c310
								+ 0.18225e5 / 4. * c111
								- 729. * c202 - 0.4455e4 / 8. * c013 - 81. * c323
								+ (2916 * c122)
								- (405 * c312)
								+ 405. / 4. * c313 + 891. / 2. * c023
								- 0.3645e4 * c121 - 729. / 8. * c233
								+ 81. * c300
								- (81 * c332)
								+ (162 * c302) + 0.3645e4 / 4. * c131
								- 405. / 2. * c301
								+ 81. / 4. * c333 + 0.3645e4 / 4. * c113
								- 0.18225e5 / 8. * c211
								+ (324 * c322) + 729. / 4. * c133
								- 0.1458e4 * c222
								+ 729. / 2. * c223 + 729. / 4. * c203
								- 729. / 2. * c200
								+ 0.3645e4 / 4. * c201
								- 891. / 8. * c033 - 0.4455e4 / 8. * c031
								+ 891. / 2. * c032
								- 729. / 2. * c130
								+ 729. * c100
								- 0.3645e4 / 8. * c213
								+ 0.3645e4 / 2. * c212
								- 0.1782e4 * c022 - 81. / 2. * c303
								+ 891. / 4. * c030 + 891. / 4. * c003 + 0.3645e4 / 4. * c210
								- 891. * c002
								+ 0.4455e4 / 4. * c001 + 0.4455e4 / 2. * c012
								- (729 * c132)
								+ 162. * c320
								- 81. / 2. * c330
								- 0.3645e4 / 2. * c110
								- 0.22275e5 / 8. * c011
					- 729. * c220 - 0.3645e4 / 2. * c101);

		// Sixth power terms (10) 3-2-1 (6), 3-3-0 (3), 2-2-2 (1)

		Cuv3w2[ii] = +(0.3645e4 / 8. * c231 - 729. / 4. * c230
										- 0.1215e4 / 4. * c311
								- 0.2673e4 / 4. * c010 - 729. / 2. * c232
								+ 891. / 4. * c000 + (2187 * c112)
								- 0.10935e5 / 8. * c221
								+ 0.2187e4 / 4. * c123
								- (729 * c102)
								- 405. / 4. * c331
								- 0.2187e4 / 2. * c120
								- 0.13365e5 / 8. * c021
								+ 0.1215e4 / 4. * c321
								+ 0.2673e4 / 4. * c020
								+ 729. / 4. * c103
								+ 0.243e3 / 2. * c310
								- 0.10935e5 / 4. * c111
								+ 729. / 2. * c202 + 0.2673e4 / 8. * c013 + 0.243e3 / 4. * c323
								- (2187 * c122)
								+ (243 * c312)
								- 0.243e3 / 4. * c313 - 0.2673e4 / 8. * c023
								+ 0.10935e5 / 4. * c121 + 729. / 8. * c233
								- 81. / 2. * c300
								+ (81 * c332)
								- (81 * c302) - 0.3645e4 / 4. * c131
								+ 405. / 4. * c301
								- 81. / 4. * c333 - 0.2187e4 / 4. * c113
								+ 0.10935e5 / 8. * c211
								- (243 * c322) - 729. / 4. * c133
								+ 0.2187e4 / 2. * c222
								- 0.2187e4 / 8. * c223 - 729. / 8. * c203
								+ 729. / 4. * c200
								- 0.3645e4 / 8. * c201
								+ 891. / 8. * c033 + 0.4455e4 / 8. * c031
								- 891. / 2. * c032
								+ 729. / 2. * c130
								- 729. / 2. * c100
								+ 0.2187e4 / 8. * c213
								- 0.2187e4 / 2. * c212
								+ 0.2673e4 / 2. * c022 + 81. / 4. * c303
								- 891. / 4. * c030 - 891. / 8. * c003 - 0.2187e4 / 4. * c210
								+ 891. / 2. * c002
								- 0.4455e4 / 8. * c001 - 0.2673e4 / 2. * c012
								+ (729 * c132)
								- 0.243e3 / 2. * c320
								+ 81. / 2. * c330
								+ 0.2187e4 / 2. * c110
								+ 0.13365e5 / 8. * c011
										+ 0.2187e4 / 4. * c220
										+ 0.3645e4 / 4. * c101);

		Cuv2w3[ii] = +(0.2187e4 / 8. * c231 - 729. / 8. * c230
										- 0.1215e4 / 4. * c311
								- 0.4455e4 / 8. * c010 - 0.2187e4 / 8. * c232
								+ 891. / 4. * c000 + 0.10935e5 / 4. * c112
								- 0.2187e4 / 2. * c221
								+ 729. * c123
								- 0.2187e4 / 2. * c102
								- 0.243e3 / 4. * c331
								- 729. * c120
								- 0.2673e4 / 2. * c021
								+ 0.243e3 * c321
								+ 891. / 2. * c020
								+ 729. / 2. * c103
								+ 405. / 4. * c310
								- 0.10935e5 / 4. * c111
								+ 0.2187e4 / 4. * c202 + 0.4455e4 / 8. * c013 + 81. * c323
								- (2187 * c122)
								+ 0.1215e4 / 4. * c312
								- 405. / 4. * c313 - 891. / 2. * c023
								+ 0.2187e4 * c121 + 729. / 8. * c233
								- 81. / 2. * c300
								+ 0.243e3 / 4. * c332
								- 0.243e3 / 2. * c302 - 0.2187e4 / 4. * c131
								+ 0.243e3 / 2. * c301
								- 81. / 4. * c333 - 0.3645e4 / 4. * c113
								+ 0.10935e5 / 8. * c211
								- (243 * c322) - 729. / 4. * c133
								+ 0.2187e4 / 2. * c222
								- 729. / 2. * c223 - 729. / 4. * c203
								+ 729. / 4. * c200
								- 0.2187e4 / 4. * c201
								+ 891. / 8. * c033 + 0.2673e4 / 8. * c031
								- 0.2673e4 / 8. * c032
								+ 729. / 4. * c130
								- 729. / 2. * c100
								+ 0.3645e4 / 8. * c213
								- 0.10935e5 / 8. * c212
								+ 0.2673e4 / 2. * c022 + 81. / 2. * c303
								- 891. / 8. * c030 - 891. / 4. * c003 - 0.3645e4 / 8. * c210
								+ 0.2673e4 / 4. * c002
								- 0.2673e4 / 4. * c001 - 0.13365e5 / 8. * c012
								+ 0.2187e4 / 4. * c132
								- 81. * c320
								+ 81. / 4. * c330
								+ 0.3645e4 / 4. * c110
								+ 0.13365e5 / 8. * c011
										+ 729. / 2. * c220
										+ 0.2187e4 / 2. * c101);

		Cu2v3w[ii] = +(729. * c231 - 891. / 2. * c230 - 0.2187e4 / 4. * c311
								- 0.2673e4 / 4. * c010 - 729. / 2. * c232
								+ 891. / 4. * c000 + 0.10935e5 / 8. * c112
								- 0.2187e4 * c221
								+ 0.1215e4 / 4. * c123
								- 0.3645e4 / 8. * c102
								- 729. / 4. * c331
								- 0.13365e5 / 8. * c120 - 0.2187e4 / 2. * c021
								+ 0.2187e4 / 4. * c321
								+ 0.2673e4 / 4. * c020
								+ 405. / 4. * c103
								+ 0.2673e4 / 8. * c310
								- 0.10935e5 / 4. * c111
								+ 729. / 2. * c202 + 0.243e3 / 2. * c013 + 0.243e3 / 4. * c323
								- 0.10935e5 / 8. * c122
								+ 0.2187e4 / 8. * c312
								- 0.243e3 / 4. * c313 - 0.243e3 / 2. * c023
								+ 0.10935e5 / 4. * c121 + 81. * c233
								- 891. / 8. * c300
								+ 729. / 8. * c332
								- 729. / 8. * c302 - 0.3645e4 / 4. * c131
								+ 729. / 4. * c301
								- 81. / 4. * c333 - 0.1215e4 / 4. * c113
								+ 0.2187e4 * c211
								- 0.2187e4 / 8. * c322 - 405. / 4. * c133
								+ 0.2187e4 / 2. * c222
								- 0.243e3 * c223 - 81. * c203
								+ 891. / 2. * c200
								- 729. * c201
								+ 81. / 2. * c033 + 729. / 2. * c031
								- 729. / 4. * c032
								+ 0.4455e4 / 8. * c130
								- 0.4455e4 / 8. * c100
								+ 0.243e3 * c213
								- 0.2187e4 / 2. * c212
								+ 0.2187e4 / 4. * c022 + 81. / 4. * c303
								- 891. / 4. * c030 - 81. / 2. * c003 - 0.2673e4 / 2. * c210
								+ 729. / 4. * c002
								- 729. / 2. * c001 - 0.2187e4 / 4. * c012
								+ 0.3645e4 / 8. * c132
								- 0.2673e4 / 8. * c320
								+ 891. / 8. * c330
								+ 0.13365e5 / 8. * c110
								+ 0.2187e4 / 2. * c011
								+ 0.2673e4 / 2. * c220
										+ 0.3645e4 / 4. * c101);

		Cu2vw3[ii] = +(0.243e3 * c231 - 81. * c230 - 0.2187e4 / 4. * c311
								- 729. / 2. * c010 - 0.243e3 * c232
								+ 891. / 4. * c000 + 0.10935e5 / 4. * c112
								- 0.2187e4 / 2. * c221
								+ 0.3645e4 / 8. * c123
								- 0.13365e5 / 8. * c102
								- 0.243e3 / 4. * c331
								- 0.3645e4 / 8. * c120 - 0.2187e4 / 4. * c021
								+ 0.2187e4 / 8. * c321
								+ 729. / 4. * c020
								+ 0.4455e4 / 8. * c103
								+ 729. / 4. * c310
								- 0.10935e5 / 4. * c111
								+ 0.2673e4 / 2. * c202 + 729. / 2. * c013 + 729. / 8. * c323
								- 0.10935e5 / 8. * c122
								+ 0.2187e4 / 4. * c312
								- 729. / 4. * c313 - 729. / 4. * c023
								+ 0.10935e5 / 8. * c121 + 81. * c233
								- 891. / 8. * c300
								+ 0.243e3 / 4. * c332
								- 0.2673e4 / 8. * c302 - 0.1215e4 / 4. * c131
								+ 0.2673e4 / 8. * c301
								- 81. / 4. * c333 - 0.3645e4 / 4. * c113
								+ 0.2187e4 * c211
								- 0.2187e4 / 8. * c322 - 405. / 4. * c133
								+ 0.2187e4 / 2. * c222
								- 729. / 2. * c223 - 891. / 2. * c203
								+ 891. / 2. * c200
								- 0.2673e4 / 2. * c201
								+ 81. / 2. * c033 + 0.243e3 / 2. * c031
								- 0.243e3 / 2. * c032
								+ 405. / 4. * c130
								- 0.4455e4 / 8. * c100
								+ 729. * c213
								- 0.2187e4 * c212
								+ 0.2187e4 / 4. * c022 + 891. / 8. * c303
								- 81. / 2. * c030 - 891. / 4. * c003 - 729. * c210
								+ 0.2673e4 / 4. * c002
								- 0.2673e4 / 4. * c001 - 0.2187e4 / 2. * c012
								+ 0.1215e4 / 4. * c132
								- 729. / 8. * c320
								+ 81. / 4. * c330
								+ 0.3645e4 / 4. * c110
								+ 0.2187e4 / 2. * c011
										+ 729. / 2. * c220
										+ 0.13365e5 / 8. * c101);

		Cu3vw2[ii] = +(0.1215e4 / 4. * c231 - 0.243e3 / 2. * c230
										- 0.3645e4 / 4. * c311
								- 729. / 2. * c010 - 0.243e3 * c232
								+ 891. / 4. * c000 + (2187 * c112)
								- 0.10935e5 / 8. * c221
								+ 0.2187e4 / 8. * c123
								- 0.2673e4 / 2. * c102
								- 405. / 4. * c331
								- 0.2187e4 / 4. * c120
								- 0.3645e4 / 8. * c021
								+ 0.3645e4 / 8. * c321
								+ 729. / 4. * c020
								+ 0.2673e4 / 8. * c103
								+ 729. / 2. * c310
								- 0.10935e5 / 4. * c111
								+ 0.2673e4 / 2. * c202 + 729. / 4. * c013 + 729. / 8. * c323
								- 0.2187e4 / 2. * c122
								+ (729 * c312)
								- 729. / 4. * c313 - 729. / 8. * c023
								+ 0.10935e5 / 8. * c121 + 0.243e3 / 4. * c233
								- 891. / 4. * c300
								+ (81 * c332)
								- 891. / 2. * c302 - 0.1215e4 / 4. * c131
								+ 0.4455e4 / 8. * c301
								- 81. / 4. * c333 - 0.2187e4 / 4. * c113
								+ 0.10935e5 / 4. * c211
								- 729. / 2. * c322 - 0.243e3 / 4. * c133
								+ 0.2187e4 / 2. * c222
								- 0.2187e4 / 8. * c223 - 0.2673e4 / 8. * c203
								+ 0.2673e4 / 4. * c200
								- 0.13365e5 / 8. * c201
								+ 81. / 4. * c033 + 405. / 4. * c031
								- 81. * c032
								+ 0.243e3 / 2. * c130
								- 0.2673e4 / 4. * c100
								+ 0.2187e4 / 4. * c213
								- 0.2187e4 * c212
								+ 729. / 2. * c022 + 891. / 8. * c303
								- 81. / 2. * c030 - 891. / 8. * c003 - 0.2187e4 / 2. * c210
								+ 891. / 2. * c002
								- 0.4455e4 / 8. * c001 - 729. * c012
								+ (243 * c132)
								- 729. / 4. * c320
								+ 81. / 2. * c330
								+ 0.2187e4 / 2. * c110
								+ 0.3645e4 / 4. * c011
										+ 0.2187e4 / 4. * c220
										+ 0.13365e5 / 8. * c101);

		Cu3v2w[ii] = +(0.2187e4 / 4. * c231 - 0.2673e4 / 8. * c230
								- 0.3645e4 / 4. * c311 - 0.4455e4 / 8. * c010
								- 0.2187e4 / 8. * c232
								+ 891. / 4. * c000 + 0.10935e5 / 8. * c112
								- 0.2187e4 * c221
								+ 0.243e3 * c123
								- 0.2187e4 / 4. * c102
								- 729. / 4. * c331
								- 0.2673e4 / 2. * c120
								- 729. * c021
								+ 729. * c321
								+ 891. / 2. * c020
								+ 0.243e3 / 2. * c103
								+ 0.4455e4 / 8. * c310
								- 0.10935e5 / 4. * c111
								+ 0.2187e4 / 4. * c202 + 405. / 4. * c013 + 81. * c323
								- 0.2187e4 / 2. * c122
								+ 0.3645e4 / 8. * c312
								- 405. / 4. * c313 - 81. * c023
								+ 0.2187e4 * c121 + 0.243e3 / 4. * c233
								- 891. / 4. * c300
								+ 729. / 8. * c332
								- 729. / 4. * c302 - 0.2187e4 / 4. * c131
								+ 729. / 2. * c301
								- 81. / 4. * c333 - 0.1215e4 / 4. * c113
								+ 0.10935e5 / 4. * c211
								- 729. / 2. * c322 - 0.243e3 / 4. * c133
								+ 0.2187e4 / 2. * c222
								- 0.243e3 * c223 - 0.243e3 / 2. * c203
								+ 0.2673e4 / 4. * c200
								- 0.2187e4 / 2. * c201
								+ 81. / 4. * c033 + 729. / 4. * c031
								- 729. / 8. * c032
								+ 0.2673e4 / 8. * c130
								- 0.2673e4 / 4. * c100
								+ 0.1215e4 / 4. * c213
								- 0.10935e5 / 8. * c212
								+ 729. / 2. * c022 + 81. / 2. * c303
								- 891. / 8. * c030 - 81. / 2. * c003 - 0.13365e5 / 8. * c210
								+ 729. / 4. * c002
								- 729. / 2. * c001 - 0.3645e4 / 8. * c012
								+ 0.2187e4 / 8. * c132
								- 891. / 2. * c320
								+ 891. / 8. * c330
								+ 0.13365e5 / 8. * c110
								+ 0.3645e4 / 4. * c011
								+ 0.2673e4 / 2. * c220
										+ 0.2187e4 / 2. * c101);

		Cu3v3[ii] = +(81. / 4. * c330 + 0.243e3 / 4. * c310 - 0.243e3 / 4. * c320
				+ 81. / 4. * c000 + 0.243e3 / 4. * c020
							- 0.243e3 / 4. * c010
							- 81. / 4. * c030
							- 729. / 4. * c210
							- 0.243e3 / 4. * c100
							+ 729. / 4. * c220
							- 0.243e3 / 4. * c230 - 81. / 4. * c300
							+ 0.243e3 / 4. * c200 + 0.243e3 / 4. * c130 + 729. / 4. * c110
									- 729. / 4. * c120);

		Cv3w3[ii] = +(81. / 4. * c000 - 81. / 4. * c003 + 0.243e3 / 4. * c002
				- 0.243e3 / 4. * c001 - 729. / 4. * c021
							+ 729. / 4. * c022 + 0.243e3 / 4. * c020
							- 729. / 4. * c012
							+ 0.243e3 / 4. * c013 + 729. / 4. * c011
							- 0.243e3 / 4. * c010 - 0.243e3 / 4. * c032
							+ 81. / 4. * c033
							- 81. / 4. * c030
							+ 0.243e3 / 4. * c031
									- 0.243e3 / 4. * c023);

		Cu3w3[ii] = +(-0.243e3 / 4. * c203 + 81. / 4. * c303 + 729. / 4. * c101
							+ 81. / 4. * c000
							- 81. / 4. * c003
							+ 0.243e3 / 4. * c002
							- 0.243e3 / 4. * c001 - 0.243e3 / 4. * c100 - 729. / 4. * c102
							+ 729. / 4. * c202
							- 81. / 4. * c300
							+ 0.243e3 / 4. * c301
							- 0.243e3 / 4. * c302 - 729. / 4. * c201
									+ 0.243e3 / 4. * c200 + 0.243e3 / 4. * c103);

		Cu2v2w2[ii] = +(0.3645e4 / 2. * c231 - 729. * c230 - 0.18225e5 / 8. * c311
								- 0.3645e4 / 2. * c010 - 0.1458e4 * c232
								+ 729. * c000 + 0.18225e5 / 2. * c112
								- 0.7290e4 * c221
								+ 0.3645e4 / 2. * c123
								- (3645 * c102)
								- 0.3645e4 / 8. * c331
								- 0.3645e4 * c120
								- 0.3645e4 * c021
								+ 0.3645e4 / 2. * c321
								+ 0.1458e4 * c020
								+ 0.3645e4 / 4. * c103
								+ 0.3645e4 / 4. * c310
								- 0.91125e5 / 8. * c111
								+ 0.2916e4 * c202 + 0.3645e4 / 4. * c013 + 729. / 2. * c323
								- (7290 * c122)
								+ 0.3645e4 / 2. * c312
								- 0.3645e4 / 8. * c313 - 729. * c023
								+ 0.18225e5 / 2. * c121 + 729. / 2. * c233
								- 729. / 2. * c300
								+ 729. / 2. * c332
								- (729 * c302) - 0.18225e5 / 8. * c131
								+ 0.3645e4 / 4. * c301
								- 729. / 8. * c333 - 0.18225e5 / 8. * c113
								+ 0.18225e5 / 2. * c211
								- (1458 * c322) - 0.3645e4 / 8. * c133
								+ 0.5832e4 * c222
								- 0.1458e4 * c223 - 729. * c203
								+ 0.1458e4 * c200
								- 0.3645e4 * c201
								+ 729. / 4. * c033 + 0.3645e4 / 4. * c031
								- 729. * c032
								+ 0.3645e4 / 4. * c130
								- 0.3645e4 / 2. * c100
								+ 0.3645e4 / 2. * c213
								- 0.7290e4 * c212
								+ 0.2916e4 * c022 + 729. / 4. * c303
								- 729. / 2. * c030 - 729. / 2. * c003 - 0.3645e4 * c210
								+ 0.1458e4 * c002
								- 0.3645e4 / 2. * c001 - 0.3645e4 * c012
								+ 0.3645e4 / 2. * c132
								- 729. * c320
								+ 729. / 4. * c330
								+ 0.18225e5 / 4. * c110
								+ 0.18225e5 / 4. * c011
								+ 0.2916e4 * c220
										+ 0.18225e5 / 4. * c101);

		// Seventh power terms (6) 3-3-1 (3), 3-2-2 (3)

		Cu3vw3[ii] = +(-729. / 4. * c231 + 0.243e3 / 4. * c230
										+ 0.2187e4 / 4. * c311
								+ 729. / 4. * c010 + 729. / 4. * c232
								- 891. / 8. * c000 - 0.6561e4 / 4. * c112
								+ 0.6561e4 / 8. * c221
								- 0.2187e4 / 8. * c123
								+ 0.8019e4 / 8. * c102
								+ 0.243e3 / 4. * c331
								+ 0.2187e4 / 8. * c120 + 0.2187e4 / 8. * c021
								- 0.2187e4 / 8. * c321
								- 729. / 8. * c020
								- 0.2673e4 / 8. * c103
								- 729. / 4. * c310
								+ 0.6561e4 / 4. * c111
								- 0.8019e4 / 8. * c202 - 729. / 4. * c013 - 729. / 8. * c323
								+ 0.6561e4 / 8. * c122
								- 0.2187e4 / 4. * c312
								+ 729. / 4. * c313 + 729. / 8. * c023
								- 0.6561e4 / 8. * c121 - 0.243e3 / 4. * c233
								+ 891. / 8. * c300
								- 0.243e3 / 4. * c332
								+ 0.2673e4 / 8. * c302 + 729. / 4. * c131
								- 0.2673e4 / 8. * c301
								+ 81. / 4. * c333 + 0.2187e4 / 4. * c113
								- 0.6561e4 / 4. * c211
								+ 0.2187e4 / 8. * c322 + 0.243e3 / 4. * c133
								- 0.6561e4 / 8. * c222
								+ 0.2187e4 / 8. * c223 + 0.2673e4 / 8. * c203
								- 0.2673e4 / 8. * c200
								+ 0.8019e4 / 8. * c201
								- 81. / 4. * c033 - 0.243e3 / 4. * c031
								+ 0.243e3 / 4. * c032
								- 0.243e3 / 4. * c130
								+ 0.2673e4 / 8. * c100
								- 0.2187e4 / 4. * c213
								+ 0.6561e4 / 4. * c212
								- 0.2187e4 / 8. * c022 - 891. / 8. * c303
								+ 81. / 4. * c030 + 891. / 8. * c003 + 0.2187e4 / 4. * c210
								- 0.2673e4 / 8. * c002
								+ 0.2673e4 / 8. * c001 + 0.2187e4 / 4. * c012
								- 729. / 4. * c132
								+ 729. / 8. * c320
								- 81. / 4. * c330
								- 0.2187e4 / 4. * c110
								- 0.2187e4 / 4. * c011
										- 0.2187e4 / 8. * c220
										- 0.8019e4 / 8. * c101);

		Cu3v3w[ii] = +(-0.2187e4 / 4. * c231 + 0.2673e4 / 8. * c230
								+ 0.2187e4 / 4. * c311 + 0.2673e4 / 8. * c010
								+ 0.2187e4 / 8. * c232
								- 891. / 8. * c000 - 0.6561e4 / 8. * c112
								+ 0.6561e4 / 4. * c221
								- 729. / 4. * c123
								+ 0.2187e4 / 8. * c102
								+ 729. / 4. * c331
								+ 0.8019e4 / 8. * c120
								+ 0.2187e4 / 4. * c021
								- 0.2187e4 / 4. * c321
								- 0.2673e4 / 8. * c020
								- 0.243e3 / 4. * c103 - 0.2673e4 / 8. * c310
								+ 0.6561e4 / 4. * c111
								- 0.2187e4 / 8. * c202
								- 0.243e3 / 4. * c013
								- 0.243e3 / 4. * c323
								+ 0.6561e4 / 8. * c122
								- 0.2187e4 / 8. * c312
								+ 0.243e3 / 4. * c313 + 0.243e3 / 4. * c023
								- 0.6561e4 / 4. * c121 - 0.243e3 / 4. * c233
								+ 891. / 8. * c300
								- 729. / 8. * c332
								+ 729. / 8. * c302 + 0.2187e4 / 4. * c131
								- 729. / 4. * c301
								+ 81. / 4. * c333 + 729. / 4. * c113
								- 0.6561e4 / 4. * c211
								+ 0.2187e4 / 8. * c322 + 0.243e3 / 4. * c133
								- 0.6561e4 / 8. * c222
								+ 729. / 4. * c223 + 0.243e3 / 4. * c203
								- 0.2673e4 / 8. * c200
								+ 0.2187e4 / 4. * c201
								- 81. / 4. * c033 - 729. / 4. * c031
								+ 729. / 8. * c032
								- 0.2673e4 / 8. * c130
								+ 0.2673e4 / 8. * c100
								- 729. / 4. * c213
								+ 0.6561e4 / 8. * c212
								- 0.2187e4 / 8. * c022 - 81. / 4. * c303
								+ 891. / 8. * c030 + 81. / 4. * c003 + 0.8019e4 / 8. * c210
								- 729. / 8. * c002
								+ 729. / 4. * c001 + 0.2187e4 / 8. * c012
								- 0.2187e4 / 8. * c132
								+ 0.2673e4 / 8. * c320
								- 891. / 8. * c330
								- 0.8019e4 / 8. * c110
								- 0.2187e4 / 4. * c011
										- 0.8019e4 / 8. * c220
										- 0.2187e4 / 4. * c101);

		Cuv3w3[ii] = +(-0.2187e4 / 8. * c231 + 729. / 8. * c230 + 729. / 4. * c311
								+ 0.2673e4 / 8. * c010 + 0.2187e4 / 8. * c232
								- 891. / 8. * c000 - 0.6561e4 / 4. * c112
								+ 0.6561e4 / 8. * c221
								- 0.2187e4 / 4. * c123
								+ 0.2187e4 / 4. * c102
								+ 0.243e3 / 4. * c331
								+ 0.2187e4 / 4. * c120 + 0.8019e4 / 8. * c021
								- 729. / 4. * c321
								- 0.2673e4 / 8. * c020
								- 729. / 4. * c103
								- 0.243e3 / 4. * c310
								+ 0.6561e4 / 4. * c111
								- 0.2187e4 / 8. * c202
								- 0.2673e4 / 8. * c013
								- 0.243e3 / 4. * c323
								+ 0.6561e4 / 4. * c122
								- 729. / 4. * c312
								+ 0.243e3 / 4. * c313 + 0.2673e4 / 8. * c023
								- 0.6561e4 / 4. * c121 - 729. / 8. * c233
								+ 81. / 4. * c300
								- 0.243e3 / 4. * c332
								+ 0.243e3 / 4. * c302 + 0.2187e4 / 4. * c131
								- 0.243e3 / 4. * c301
								+ 81. / 4. * c333 + 0.2187e4 / 4. * c113
								- 0.6561e4 / 8. * c211
								+ 729. / 4. * c322 + 729. / 4. * c133
								- 0.6561e4 / 8. * c222
								+ 0.2187e4 / 8. * c223 + 729. / 8. * c203
								- 729. / 8. * c200
								+ 0.2187e4 / 8. * c201
								- 891. / 8. * c033 - 0.2673e4 / 8. * c031
								+ 0.2673e4 / 8. * c032
								- 729. / 4. * c130
								+ 729. / 4. * c100
								- 0.2187e4 / 8. * c213
								+ 0.6561e4 / 8. * c212
								- 0.8019e4 / 8. * c022 - 81. / 4. * c303
								+ 891. / 8. * c030 + 891. / 8. * c003 + 0.2187e4 / 8. * c210
								- 0.2673e4 / 8. * c002
								+ 0.2673e4 / 8. * c001 + 0.8019e4 / 8. * c012
								- 0.2187e4 / 0.4e1 * c132
								+ 0.243e3 / 0.4e1 * c320
								- 0.81e2 / 0.4e1 * c330
								- 0.2187e4 / 0.4e1 * c110
								- 0.8019e4 / 8. * c011
								- 0.2187e4 / 8. * c220
										- 0.2187e4 / 0.4e1 * c101);

		Cu3v2w2[ii] = +(-0.10935e5 / 8. * c231 + 0.2187e4 / 0.4e1 * c230
								+ 0.18225e5 / 8. * c311 + 0.3645e4 / 0.4e1 * c010
								+ 0.2187e4 / 2. * c232
								- 729. / 2. * c000 - 0.10935e5 / 2. * c112
								+ 0.10935e5 / 2. * c221
								- 0.2187e4 / 2. * c123
								+ (2187 * c102)
								+ 0.3645e4 / 8. * c331
								+ 0.2187e4 * c120
								+ 0.3645e4 / 2. * c021
								- 0.3645e4 / 2. * c321
								- 729. * c020
								- 0.2187e4 / 0.4e1 * c103
								- 0.3645e4 / 0.4e1 * c310
								+ 0.54675e5 / 8. * c111
								- 0.2187e4 * c202 - 0.3645e4 / 8. * c013 - 729. / 2. * c323
								+ (4374 * c122)
								- 0.3645e4 / 2. * c312
								+ 0.3645e4 / 8. * c313 + 729. / 2. * c023
								- 0.10935e5 / 2. * c121 - 0.2187e4 / 8. * c233
								+ 729. / 2. * c300
								- 729. / 2. * c332
								+ (729 * c302) + 0.10935e5 / 8. * c131
								- 0.3645e4 / 0.4e1 * c301
								+ 729. / 8. * c333 + 0.10935e5 / 8. * c113
								- 0.54675e5 / 8. * c211
								+ (1458 * c322) + 0.2187e4 / 8. * c133
								- 0.4374e4 * c222
								+ 0.2187e4 / 2. * c223 + 0.2187e4 / 0.4e1 * c203
								- 0.2187e4 / 2. * c200
								+ 0.10935e5 / 0.4e1 * c201
								- 729. / 8. * c033 - 0.3645e4 / 8. * c031
								+ 729. / 2. * c032
								- 0.2187e4 / 0.4e1 * c130
								+ 0.2187e4 / 2. * c100
								- 0.10935e5 / 8. * c213
								+ 0.10935e5 / 2. * c212
								- 0.1458e4 * c022 - 729. / 0.4e1 * c303
								+ 729. / 0.4e1 * c030
								+ 729. / 0.4e1 * c003
								+ 0.10935e5 / 0.4e1 * c210
								- 729. * c002
								+ 0.3645e4 / 0.4e1 * c001 + 0.3645e4 / 2. * c012
								- 0.2187e4 / 2. * c132
								+ 729. * c320
								- 729. / 0.4e1 * c330
								- 0.10935e5 / 0.4e1 * c110
								- 0.18225e5 / 8. * c011
								- 0.2187e4 * c220
										- 0.10935e5 / 0.4e1 * c101);

		Cu2v3w2[ii] = +(-0.3645e4 / 2. * c231 + 729. * c230 + 0.10935e5 / 8. * c311
								+ 0.2187e4 / 2. * c010 + 0.1458e4 * c232
								- 729. / 2. * c000 - 0.10935e5 / 2. * c112
								+ 0.10935e5 / 2. * c221
								- 0.10935e5 / 8. * c123
								+ 0.3645e4 / 2. * c102
								+ 0.3645e4 / 8. * c331
								+ 0.10935e5 / 0.4e1 * c120 + 0.10935e5 / 0.4e1 * c021
								- 0.10935e5 / 8. * c321
								- 0.2187e4 / 2. * c020
								- 0.3645e4 / 8. * c103 - 0.2187e4 / 0.4e1 * c310
								+ 0.54675e5 / 8. * c111
								- 0.1458e4 * c202
								- 0.2187e4 / 0.4e1 * c013
								- 0.2187e4 / 8. * c323
								+ 0.10935e5 / 2. * c122
								- 0.2187e4 / 2. * c312
								+ 0.2187e4 / 8. * c313 + 0.2187e4 / 0.4e1 * c023
								- 0.54675e5 / 8. * c121 - 729. / 2. * c233
								+ 729. / 0.4e1 * c300
								- 729. / 2. * c332
								+ 729. / 2. * c302 + 0.18225e5 / 8. * c131
								- 0.3645e4 / 8. * c301
								+ 729. / 8. * c333 + 0.10935e5 / 8. * c113
								- 0.10935e5 / 2. * c211
								+ 0.2187e4 / 2. * c322 + 0.3645e4 / 8. * c133
								- 0.4374e4 * c222
								+ 0.2187e4 / 2. * c223 + 729. / 2. * c203
								- 729. * c200
								+ 0.3645e4 / 2. * c201
								- 729. / 0.4e1 * c033 - 0.3645e4 / 0.4e1 * c031
								+ 729. * c032
								- 0.3645e4 / 0.4e1 * c130
								+ 0.3645e4 / 0.4e1 * c100
								- 0.2187e4 / 2. * c213
								+ 0.4374e4 * c212
								- 0.2187e4 * c022 - 729. / 8. * c303
								+ 729. / 2. * c030 + 729. / 0.4e1 * c003 + 0.2187e4 * c210
								- 729. * c002
								+ 0.3645e4 / 0.4e1 * c001 + 0.2187e4 * c012
								- 0.3645e4 / 2. * c132
								+ 0.2187e4 / 0.4e1 * c320
								- 729. / 0.4e1 * c330
								- 0.10935e5 / 0.4e1 * c110
								- 0.10935e5 / 0.4e1 * c011
								- 0.2187e4 * c220
										- 0.18225e5 / 8. * c101);

		Cu2v2w3[ii] = +(-0.2187e4 / 2. * c231 + 729. / 2. * c230
										+ 0.10935e5 / 8. * c311
								+ 0.3645e4 / 0.4e1 * c010 + 0.2187e4 / 2. * c232
								- 729. / 2. * c000 - 0.54675e5 / 8. * c112
								+ 0.4374e4 * c221
								- 0.3645e4 / 2. * c123
								+ 0.10935e5 / 0.4e1 * c102
								+ 0.2187e4 / 8. * c331
								+ 0.3645e4 / 2. * c120 + 0.2187e4 * c021
								- 0.2187e4 / 2. * c321
								- 729. * c020
								- 0.3645e4 / 0.4e1 * c103
								- 0.3645e4 / 8. * c310
								+ 0.54675e5 / 8. * c111
								- 0.2187e4 * c202 - 0.3645e4 / 0.4e1 * c013 - 729. / 2. * c323
								+ 0.10935e5 / 2. * c122
								- 0.10935e5 / 8. * c312
								+ 0.3645e4 / 8. * c313 + 729. * c023
								- 0.10935e5 / 2. * c121 - 729. / 2. * c233
								+ 729. / 0.4e1 * c300
								- 0.2187e4 / 8. * c332
								+ 0.2187e4 / 0.4e1 * c302 + 0.10935e5 / 8. * c131
								- 0.2187e4 / 0.4e1 * c301
								+ 729. / 8. * c333 + 0.18225e5 / 8. * c113
								- 0.10935e5 / 2. * c211
								+ 0.2187e4 / 2. * c322 + 0.3645e4 / 8. * c133
								- 0.4374e4 * c222
								+ 0.1458e4 * c223 + 729. * c203
								- 729. * c200
								+ 0.2187e4 * c201
								- 729. / 0.4e1 * c033 - 0.2187e4 / 0.4e1 * c031
								+ 0.2187e4 / 0.4e1 * c032
								- 0.3645e4 / 8. * c130
								+ 0.3645e4 / 0.4e1 * c100
								- 0.3645e4 / 2. * c213
								+ 0.10935e5 / 2. * c212
								- 0.2187e4 * c022 - 729. / 0.4e1 * c303
								+ 729. / 0.4e1 * c030 + 729. / 2. * c003 + 0.3645e4 / 2. * c210
								- 0.2187e4 / 2. * c002
								+ 0.2187e4 / 2. * c001 + 0.10935e5 / 0.4e1 * c012
								- 0.10935e5 / 8. * c132
								+ 729. / 2. * c320
								- 729. / 8. * c330
								- 0.18225e5 / 8. * c110
								- 0.10935e5 / 0.4e1 * c011
								- 0.1458e4 * c220
										- 0.10935e5 / 0.4e1 * c101);

		// Eighth power terms (3) 3-3-2

		Cu3v3w2[ii] = +(0.10935e5 / 8. * c231 - 0.2187e4 / 0.4e1 * c230
								- 0.10935e5 / 8. * c311 - 0.2187e4 / 0.4e1 * c010
								- 0.2187e4 / 2. * c232
								+ 729. / 0.4e1 * c000 + 0.6561e4 / 2. * c112
								- 0.32805e5 / 8. * c221
								+ 0.6561e4 / 8. * c123
								- 0.2187e4 / 2. * c102
								- 0.3645e4 / 8. * c331
								- 0.6561e4 / 0.4e1 * c120 - 0.10935e5 / 8. * c021
								+ 0.10935e5 / 8. * c321
								+ 0.2187e4 / 0.4e1 * c020
								+ 0.2187e4 / 8. * c103 + 0.2187e4 / 0.4e1 * c310
								- 0.32805e5 / 8. * c111
								+ 0.2187e4 / 2. * c202
								+ 0.2187e4 / 8. * c013
								+ 0.2187e4 / 8. * c323
								- 0.6561e4 / 2. * c122
								+ 0.2187e4 / 2. * c312
								- 0.2187e4 / 8. * c313 - 0.2187e4 / 8. * c023
								+ 0.32805e5 / 8. * c121 + 0.2187e4 / 8. * c233
								- 729. / 0.4e1 * c300
								+ 729. / 2. * c332
								- 729. / 2. * c302 - 0.10935e5 / 8. * c131
								+ 0.3645e4 / 8. * c301
								- 729. / 8. * c333 - 0.6561e4 / 8. * c113
								+ 0.32805e5 / 8. * c211
								- 0.2187e4 / 2. * c322 - 0.2187e4 / 8. * c133
								+ 0.6561e4 / 2. * c222
								- 0.6561e4 / 8. * c223 - 0.2187e4 / 8. * c203
								+ 0.2187e4 / 0.4e1 * c200
								- 0.10935e5 / 8. * c201
								+ 729. / 8. * c033 + 0.3645e4 / 8. * c031
								- 729. / 2. * c032
								+ 0.2187e4 / 0.4e1 * c130
								- 0.2187e4 / 0.4e1 * c100
								+ 0.6561e4 / 8. * c213
								- 0.6561e4 / 2. * c212
								+ 0.2187e4 / 2. * c022 + 729. / 8. * c303
								- 729. / 0.4e1 * c030
								- 729. / 8. * c003
								- 0.6561e4 / 0.4e1 * c210
								+ 729. / 2. * c002
								- 0.3645e4 / 8. * c001 - 0.2187e4 / 2. * c012
								+ 0.2187e4 / 2. * c132
								- 0.2187e4 / 0.4e1 * c320
								+ 729. / 0.4e1 * c330
								+ 0.6561e4 / 0.4e1 * c110
								+ 0.10935e5 / 8. * c011
								+ 0.6561e4 / 0.4e1 * c220
										+ 0.10935e5 / 8. * c101);

		Cu3v2w3[ii] = +(0.6561e4 / 8. * c231 - 0.2187e4 / 8. * c230
								- 0.10935e5 / 8. * c311 - 0.3645e4 / 8. * c010
								- 0.6561e4 / 8. * c232
								+ 729. / 0.4e1 * c000 + 0.32805e5 / 8. * c112
								- 0.6561e4 / 2. * c221
								+ 0.2187e4 / 2. * c123
								- 0.6561e4 / 0.4e1 * c102
								- 0.2187e4 / 8. * c331
								- 0.2187e4 / 2. * c120 - 0.2187e4 / 2. * c021
								+ 0.2187e4 / 2. * c321
								+ 729. / 2. * c020
								+ 0.2187e4 / 0.4e1 * c103 + 0.3645e4 / 8. * c310
								- 0.32805e5 / 8. * c111
								+ 0.6561e4 / 0.4e1 * c202
								+ 0.3645e4 / 8. * c013
								+ 729. / 2. * c323
								- 0.6561e4 / 2. * c122
								+ 0.10935e5 / 8. * c312
								- 0.3645e4 / 8. * c313 - 729. / 2. * c023
								+ 0.6561e4 / 2. * c121 + 0.2187e4 / 8. * c233
								- 729. / 0.4e1 * c300
								+ 0.2187e4 / 8. * c332
								- 0.2187e4 / 0.4e1 * c302 - 0.6561e4 / 8. * c131
								+ 0.2187e4 / 0.4e1 * c301
								- 729. / 8. * c333 - 0.10935e5 / 8. * c113
								+ 0.32805e5 / 8. * c211
								- 0.2187e4 / 2. * c322 - 0.2187e4 / 8. * c133
								+ 0.6561e4 / 2. * c222
								- 0.2187e4 / 2. * c223 - 0.2187e4 / 0.4e1 * c203
								+ 0.2187e4 / 0.4e1 * c200
								- 0.6561e4 / 0.4e1 * c201
								+ 729. / 8. * c033 + 0.2187e4 / 8. * c031
								- 0.2187e4 / 8. * c032
								+ 0.2187e4 / 8. * c130
								- 0.2187e4 / 0.4e1 * c100
								+ 0.10935e5 / 8. * c213
								- 0.32805e5 / 8. * c212
								+ 0.2187e4 / 2. * c022 + 729. / 0.4e1 * c303
								- 729. / 8. * c030 - 729. / 0.4e1 * c003 - 0.10935e5 / 8. * c210
								+ 0.2187e4 / 0.4e1 * c002
								- 0.2187e4 / 0.4e1 * c001 - 0.10935e5 / 8. * c012
								+ 0.6561e4 / 8. * c132
								- 729. / 2. * c320
								+ 729. / 8. * c330
								+ 0.10935e5 / 8. * c110
								+ 0.10935e5 / 8. * c011
								+ 0.2187e4 / 2. * c220
										+ 0.6561e4 / 0.4e1 * c101);

		Cu2v3w3[ii] = +(0.2187e4 / 2. * c231 - 729. / 2. * c230
										- 0.6561e4 / 8. * c311
								- 0.2187e4 / 0.4e1 * c010 - 0.2187e4 / 2. * c232
								+ 729. / 0.4e1 * c000 + 0.32805e5 / 8. * c112
								- 0.6561e4 / 2. * c221
								+ 0.10935e5 / 8. * c123
								- 0.10935e5 / 8. * c102
								- 0.2187e4 / 8. * c331
								- 0.10935e5 / 8. * c120 - 0.6561e4 / 0.4e1 * c021
								+ 0.6561e4 / 8. * c321
								+ 0.2187e4 / 0.4e1 * c020
								+ 0.3645e4 / 8. * c103 + 0.2187e4 / 8. * c310
								- 0.32805e5 / 8. * c111
								+ 0.2187e4 / 2. * c202
								+ 0.2187e4 / 0.4e1 * c013
								+ 0.2187e4 / 8. * c323
								- 0.32805e5 / 8. * c122
								+ 0.6561e4 / 8. * c312
								- 0.2187e4 / 8. * c313 - 0.2187e4 / 0.4e1 * c023
								+ 0.32805e5 / 8. * c121 + 729. / 2. * c233
								- 729. / 8. * c300
								+ 0.2187e4 / 8. * c332
								- 0.2187e4 / 8. * c302 - 0.10935e5 / 8. * c131
								+ 0.2187e4 / 8. * c301
								- 729. / 8. * c333 - 0.10935e5 / 8. * c113
								+ 0.6561e4 / 2. * c211
								- 0.6561e4 / 8. * c322 - 0.3645e4 / 8. * c133
								+ 0.6561e4 / 2. * c222
								- 0.2187e4 / 2. * c223 - 729. / 2. * c203
								+ 729. / 2. * c200
								- 0.2187e4 / 2. * c201
								+ 729. / 0.4e1 * c033 + 0.2187e4 / 0.4e1 * c031
								- 0.2187e4 / 0.4e1 * c032
								+ 0.3645e4 / 8. * c130
								- 0.3645e4 / 8. * c100
								+ 0.2187e4 / 2. * c213
								- 0.6561e4 / 2. * c212
								+ 0.6561e4 / 0.4e1 * c022 + 729. / 8. * c303
								- 729. / 0.4e1 * c030
								- 729. / 0.4e1 * c003
								- 0.2187e4 / 2. * c210
								+ 0.2187e4 / 0.4e1 * c002
								- 0.2187e4 / 0.4e1 * c001 - 0.6561e4 / 0.4e1 * c012
								+ 0.10935e5 / 8. * c132
								- 0.2187e4 / 8. * c320
								+ 729. / 8. * c330
								+ 0.10935e5 / 8. * c110
								+ 0.6561e4 / 0.4e1 * c011
								+ 0.2187e4 / 2. * c220
										+ 0.10935e5 / 8. * c101);

		// Ninth power terms (1) 3-3-3
		Cu3v3w3[ii] = +(-0.6561e4 / 8. * c231 + 0.2187e4 / 8. * c230
								+ 0.6561e4 / 8. * c311 + 0.2187e4 / 8. * c010
								+ 0.6561e4 / 8. * c232
								- 729. / 8. * c000 - 0.19683e5 / 8. * c112
								+ 0.19683e5 / 8. * c221
								- 0.6561e4 / 8. * c123
								+ 0.6561e4 / 8. * c102
								+ 0.2187e4 / 8. * c331
								+ 0.6561e4 / 8. * c120 + 0.6561e4 / 8. * c021
								- 0.6561e4 / 8. * c321
								- 0.2187e4 / 8. * c020
								- 0.2187e4 / 8. * c103 - 0.2187e4 / 8. * c310
								+ 0.19683e5 / 8. * c111
								- 0.6561e4 / 8. * c202
								- 0.2187e4 / 8. * c013
								- 0.2187e4 / 8. * c323
								+ 0.19683e5 / 8. * c122
								- 0.6561e4 / 8. * c312
								+ 0.2187e4 / 8. * c313 + 0.2187e4 / 8. * c023
								- 0.19683e5 / 8. * c121 - 0.2187e4 / 8. * c233
								+ 729. / 8. * c300
								- 0.2187e4 / 8. * c332
								+ 0.2187e4 / 8. * c302 + 0.6561e4 / 8. * c131
								- 0.2187e4 / 8. * c301
								+ 729. / 8. * c333 + 0.6561e4 / 8. * c113
								- 0.19683e5 / 8. * c211
								+ 0.6561e4 / 8. * c322 + 0.2187e4 / 8. * c133
								- 0.19683e5 / 8. * c222
								+ 0.6561e4 / 8. * c223 + 0.2187e4 / 8. * c203
								- 0.2187e4 / 8. * c200
								+ 0.6561e4 / 8. * c201
								- 729. / 8. * c033 - 0.2187e4 / 8. * c031
								+ 0.2187e4 / 8. * c032
								- 0.2187e4 / 8. * c130
								+ 0.2187e4 / 8. * c100
								- 0.6561e4 / 8. * c213
								+ 0.19683e5 / 8. * c212
								- 0.6561e4 / 8. * c022 - 729. / 8. * c303
								+ 729. / 8. * c030 + 729. / 8. * c003 + 0.6561e4 / 8. * c210
								- 0.2187e4 / 8. * c002
								+ 0.2187e4 / 8. * c001 + 0.6561e4 / 8. * c012
								- 0.6561e4 / 8. * c132
								+ 0.2187e4 / 8. * c320
								- 729. / 8. * c330
								- 0.6561e4 / 8. * c110
								- 0.6561e4 / 8. * c011
										- 0.6561e4 / 8. * c220
										- 0.6561e4 / 8. * c101);

	}
}


//void LagrangeCubicHexMapping::computeTransformedCoords(const double uvw[3],
//		double xyz[3]) const {
//	const double& u = uvw[0];
//	const double& v = uvw[1];
//	const double& w = uvw[2];
//
//	double segment0 = -4.5 * (w - 1) * (w - 2. / 3) * (w - 1. / 3);
//	double segment1 = 13.5 * (w - 1) * (w - 2. / 3) * w;
//	double segment2 = -13.5 * (w - 1) * (w - 1. / 3) * w;
//	double segment3 = 4.5 * (w - 2. / 3) * (w - 1. / 3) * w;
//
//	for (int ii = 0; ii < 3; ii++) {
//		double val0 = C0[ii]
//				+ u * (Cu0[ii]
//						+ u * (Cuu0[ii]
//								+ u * (Cuuu0[ii]
//										+ v * (Cuuuv0[ii] + v * (Cuuuvv0[ii] + v * Cuuuvvv0[ii])))
//								+ v * (Cuuv0[ii] + v * (Cuuvv0[ii] + v * Cuuvvv0[ii])))
//						+ v * (Cuv0[ii] + v * (Cuvv0[ii] + v * Cuvvv0[ii])))
//				+ v * (Cv0[ii] + v * (Cvv0[ii] + v * Cvvv0[ii]));
//		double val1 = C1[ii]
//				+ u * (Cu1[ii]
//						+ u * (Cuu1[ii]
//								+ u * (Cuuu1[ii]
//										+ v * (Cuuuv1[ii] + v * (Cuuuvv1[ii] + v * Cuuuvvv1[ii])))
//								+ v * (Cuuv1[ii] + v * (Cuuvv1[ii] + v * Cuuvvv1[ii])))
//						+ v * (Cuv1[ii] + v * (Cuvv1[ii] + v * Cuvvv1[ii])))
//				+ v * (Cv1[ii] + v * (Cvv1[ii] + v * Cvvv1[ii]));
//		double val2 = C2[ii]
//				+ u * (Cu2[ii]
//						+ u * (Cuu2[ii]
//								+ u * (Cuuu2[ii]
//										+ v * (Cuuuv2[ii] + v * (Cuuuvv2[ii] + v * Cuuuvvv2[ii])))
//								+ v * (Cuuv2[ii] + v * (Cuuvv2[ii] + v * Cuuvvv2[ii])))
//						+ v * (Cuv2[ii] + v * (Cuvv2[ii] + v * Cuvvv2[ii])))
//				+ v * (Cv2[ii] + v * (Cvv2[ii] + v * Cvvv2[ii]));
//		double val3 = C3[ii]
//				+ u * (Cu3[ii]
//						+ u * (Cuu3[ii]
//								+ u * (Cuuu3[ii]
//										+ v * (Cuuuv3[ii] + v * (Cuuuvv3[ii] + v * Cuuuvvv3[ii])))
//								+ v * (Cuuv3[ii] + v * (Cuuvv3[ii] + v * Cuuvvv3[ii])))
//						+ v * (Cuv3[ii] + v * (Cuvv3[ii] + v * Cuvvv3[ii])))
//				+ v * (Cv3[ii] + v * (Cvv3[ii] + v * Cvvv3[ii]));
//		xyz[ii] = segment0 * val0 + segment1 * val1 + segment2 * val2
//							+ segment3 * val3;
//	}
//}

void LagrangeCubicHexMapping::computeTransformedCoords(
		const double uvwCoords[3], double xyz[3]) const {
	const double& u = uvwCoords[0];
	const double& v = uvwCoords[1];
	const double& w = uvwCoords[2];

	// Construct all the monomials
	const double u2 = u * u;
	const double v2 = v * v;
	const double w2 = w * w;
	const double uv = u * v;
	const double vw = v * w;
	const double uw = u * w;

	const double u3 = u2 * u;
	const double u2v = u2 * v;
	const double uv2 = v2 * u;
	const double v3 = v2 * v;
	const double v2w = v2 * w;
	const double vw2 = w2 * v;
	const double w3 = w2 * w;
	const double uw2 = w2 * u;
	const double u2w = u2 * w;
	const double uvw = uv * w;

	const double u2v2 = u2 * v2;
	const double v2w2 = v2 * w2;
	const double u2w2 = u2 * w2;
	const double u2vw = u2 * vw;
	const double uv2w = v2 * uw;
	const double uvw2 = w2 * uv;
	const double u3v = u3 * v;
	const double u3w = u3 * w;
	const double uv3 = v3 * u;
	const double v3w = v3 * w;
	const double uw3 = w3 * u;
	const double vw3 = w3 * v;

	const double u3vw = u3 * vw;
	const double uv3w = v3 * uw;
	const double uvw3 = w3 * uv;
	const double u3v2 = u3 * v2;
	const double u3w2 = u3 * w2;
	const double u2v3 = v3 * u2;
	const double v3w2 = v3 * w2;
	const double u2w3 = w3 * u2;
	const double v2w3 = w3 * v2;
	const double u2v2w = u2v2 * w;
	const double u2vw2 = u2w2 * v;
	const double uv2w2 = v2w2 * u;

	const double u3v2w = u3v2 * w;
	const double u3vw2 = u3w2 * v;
	const double u2v3w = u2v3 * w;
	const double uv3w2 = v3w2 * u;
	const double u2vw3 = u2w3 * v;
	const double uv2w3 = v2w3 * u;
	const double u3v3 = u3v2 * v;
	const double v3w3 = v3w2 * w;
	const double u3w3 = u3w2 * w;
	const double u2v2w2 = u2v2w * w;

	const double u3v3w = u3v3 * w;
	const double u3vw3 = u3vw2 * w;
	const double uv3w3 = uv3w2 * w;
	const double u3v2w2 = u3v2w * w;
	const double u2v3w2 = u2v3w * w;
	const double u2v2w3 = u2v2w2 * w;

	const double u3v3w2 = u3v3w * w;
	const double u3v2w3 = u3v2w2 * 2;
	const double u2v3w3 = u2v3w2 * w;

	const double u3v3w3 = u3v3w2 * w;

	// Yikes!  64 monomials all rolled into one gigantic polynomial!!
	for (int ii = 0; ii < 3; ii++) {
		xyz[ii] = C[ii]
				+ (Cu[ii] * u + Cv[ii] * v + Cw[ii] * w)
				+ (Cu2[ii] * u2 + Cuv[ii] * uv + Cv2[ii] * v2 + Cvw[ii] * vw
						+ Cw2[ii] * w2 + Cuw[ii] * uw)
				+ (Cu3[ii] * u3 + Cv3[ii] * v3 + Cw3[ii] * w3 + Cu2v[ii] * u2v
						+ Cuv2[ii] * uv2 + Cv2w[ii] * v2w + Cvw2[ii] * vw2 + Cu2w[ii] * u2w
						+ Cuw2[ii] * uw2 + Cuvw[ii] * uvw)
				+ (Cu2v2[ii] * u2v2 + Cv2w2[ii] * v2w2 + Cu2w2[ii] * u2w2
						+ Cu2vw[ii] * u2vw + Cuv2w[ii] * uv2w + Cuvw2[ii] * uvw2
						+ Cu3v[ii] * u3v + Cu3w[ii] * u3w + Cuv3[ii] * uv3 + Cv3w[ii] * v3w
						+ Cuw3[ii] * uw3 + Cvw3[ii] * vw3)
				+ (Cu3vw[ii] * u3vw + Cuv3w[ii] * uv3w + Cuvw3[ii] * uvw3
						+ Cu3v2[ii] * u3v2 + Cu3w2[ii] * u3w2 + Cu2v3[ii] * u2v3
						+ Cv3w2[ii] * v3w2 + Cu2w3[ii] * u2w3 + Cv2w3[ii] * v2w3
						+ Cu2v2w[ii] * u2v2w + Cu2vw2[ii] * u2vw2 + Cuv2w2[ii] * uv2w2)
				+ (Cu3v2w[ii] * u3v2w + Cu3vw2[ii] * u3vw2 + Cu2v3w[ii] * u2v3w
						+ Cuv3w2[ii] * uv3w2 + Cu2vw3[ii] * u2vw3 + Cuv2w3[ii] * uv2w3
						+ Cu3v3[ii] * u3v3 + Cv3w3[ii] * v3w3 + Cu3w3[ii] * u3w3
						+ Cu2v2w2[ii] * u2v2w2)
				+ (Cu3v3w[ii] * u3v3w + Cu3vw3[ii] * u3vw3 + Cuv3w3[ii] * uv3w3
						+ Cu3v2w2[ii] * u3v2w2 + Cu2v3w2[ii] * u2v3w2 + Cu2v2w3[ii] * u2v2w3)
				+ (Cu3v3w2[ii] * u3v3w2 + Cu3v2w3[ii] * u3v2w3 + Cu2v3w3[ii] * u2v3w3)
				+ Cu3v3w3[ii] * u3v3w3;
	}
}

