#include "Mapping.h"

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

