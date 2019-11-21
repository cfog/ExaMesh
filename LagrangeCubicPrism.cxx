#include "Mapping.h"

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

