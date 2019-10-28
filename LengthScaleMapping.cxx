/*
 * LengthScaleMapping.cxx
 *
 *  Created on: Oct. 1, 2019
 *      Author: cfog
 */

#include <ExaMesh.h>
#include <algorithm>

#include "Mapping.h"

double LengthScaleMapping::getIsoLengthScale(const emInt vert) {
	return m_pMesh->getLengthScale(vert);
}


void TetLengthScaleMapping::setPolyCoeffs(const double* xyz0,
		const double* xyz1, const double* xyz2, const double* xyz3,
		double uderiv0[3], double vderiv0[3], double wderiv0[3], double uderiv1[3],
		double vderiv1[3], double wderiv1[3], double uderiv2[3], double vderiv2[3],
		double wderiv2[3], double uderiv3[3], double vderiv3[3],
		double wderiv3[3]) {
	for (int ii = 0; ii < 3; ii++) {
		//		xyzOffset[ii] = xyz0[ii];
		//		uVec[ii] = xyz1[ii] - xyz0[ii];
		//		vVec[ii] = xyz2[ii] - xyz0[ii];
		//		wVec[ii] = xyz3[ii] - xyz0[ii];
		//
		A[ii] = -2 * xyz1[ii] + 2 * xyz0[ii] + uderiv0[ii] + uderiv1[ii];
		B[ii] = vderiv1[ii] - vderiv0[ii];
		C[ii] = uderiv2[ii] - uderiv0[ii];
		E[ii] = -2 * xyz2[ii] + 2 * xyz0[ii] + vderiv0[ii] + vderiv2[ii];
		F[ii] = wderiv2[ii] - wderiv0[ii];
		G[ii] = vderiv3[ii] - vderiv0[ii];
		H[ii] = -2 * xyz3[ii] + 2 * xyz0[ii] + wderiv0[ii] + wderiv3[ii];
		J[ii] = uderiv3[ii] - uderiv0[ii];
		K[ii] = wderiv1[ii] - wderiv0[ii];
		// Whoever wrote the line wrap code for eclipse was smoking something they shouldn't have.
		L[ii] =
				(-6 * xyz0[ii] + 2 * xyz1[ii] + 2 * xyz2[ii] + 2 * xyz3[ii]) + (-2
						* uderiv0[ii]
																																				- uderiv1[ii]
																																				+ 0.5 * uderiv2[ii]
																																				+ 0.5 * uderiv3[ii])
				+ (-2 * vderiv0[ii] + 0.5 * vderiv1[ii] - vderiv2[ii] + 0.5
						* vderiv3[ii])
				+ (-2 * wderiv0[ii] + 0.5 * wderiv1[ii] + 0.5 * wderiv2[ii] - wderiv3[ii]);
		M[ii] = -uderiv1[ii] + 3 * xyz1[ii] - 3 * xyz0[ii] - 2 * uderiv0[ii];
		P[ii] = -vderiv2[ii] + 3 * xyz2[ii] - 3 * xyz0[ii] - 2 * vderiv0[ii];
		R[ii] = -wderiv3[ii] + 3 * xyz3[ii] - 3 * xyz0[ii] - 2 * wderiv0[ii];
		T[ii] = xyz0[ii];
		U[ii] = uderiv0[ii];
		V[ii] = vderiv0[ii];
		W[ii] = wderiv0[ii];
	}
}

void TetLengthScaleMapping::setupCoordMapping(const emInt verts[]) {
	double xyz0[3], xyz1[3], xyz2[3], xyz3[3];
	m_pMesh->getCoords(verts[0], xyz0);
	m_pMesh->getCoords(verts[1], xyz1);
	m_pMesh->getCoords(verts[2], xyz2);
	m_pMesh->getCoords(verts[3], xyz3);

	double len0 = getIsoLengthScale(verts[0]);
	double len1 = getIsoLengthScale(verts[1]);
	double len2 = getIsoLengthScale(verts[2]);
	double len3 = getIsoLengthScale(verts[3]);

	double vec01[] = DIFF(xyz1, xyz0);
	double vec02[] = DIFF(xyz2, xyz0);
	double vec03[] = DIFF(xyz3, xyz0);
	double vec12[] = DIFF(xyz2, xyz1);
	double vec13[] = DIFF(xyz3, xyz1);
	double vec23[] = DIFF(xyz3, xyz2);

	double ratio01 = sqrt(std::max(0.5, std::min(2., len0/len1)));
	double ratio02 = sqrt(std::max(0.5, std::min(2., len0/len2)));
	double ratio03 = sqrt(std::max(0.5, std::min(2., len0/len3)));
	double ratio12 = sqrt(std::max(0.5, std::min(2., len1/len2)));
	double ratio13 = sqrt(std::max(0.5, std::min(2., len1/len3)));
	double ratio23 = sqrt(std::max(0.5, std::min(2., len2/len3)));

	double uderiv0[3],
	vderiv0[3], wderiv0[3];
	double uderiv1[3], vderiv1[3], wderiv1[3];
	double uderiv2[3], vderiv2[3], wderiv2[3];
	double uderiv3[3], vderiv3[3], wderiv3[3];

	SCALE(vec01, ratio01, uderiv0);
	SCALE(vec02, ratio02, vderiv0);
	SCALE(vec03, ratio03, wderiv0);

	SCALE(vec01, 1 / ratio01, uderiv1);
	for (int ii = 0; ii < 3; ii++) {
		vderiv1[ii] = uderiv1[ii] + vec12[ii] * ratio12;
		wderiv1[ii] = uderiv1[ii] + vec13[ii] * ratio13;
	}

	SCALE(vec02, 1 / ratio02, vderiv2);
	for (int ii = 0; ii < 3; ii++) {
		uderiv2[ii] = vderiv2[ii] - vec12[ii] / ratio12;
		wderiv2[ii] = vderiv2[ii] + vec23[ii] * ratio23;
	}

	SCALE(vec03, 1 / ratio03, wderiv3);
	for (int ii = 0; ii < 3; ii++) {
		uderiv3[ii] = wderiv3[ii] - vec13[ii] / ratio13;
		vderiv3[ii] = wderiv3[ii] - vec23[ii] / ratio23;
	}

	setPolyCoeffs(xyz0, xyz1, xyz2, xyz3, uderiv0, vderiv0, wderiv0, uderiv1,
								vderiv1, wderiv1, uderiv2, vderiv2, wderiv2, uderiv3, vderiv3,
								wderiv3);
}

void TetLengthScaleMapping::computeTransformedCoords(const double uvw[3],
		double xyz[3]) const {
	const double& u = uvw[0];
	const double& v = uvw[1];
	const double& w = uvw[2];
	for (int ii = 0; ii < 3; ii++) {
		xyz[ii] = xyzOffset[ii] + u * uVec[ii] + v * vVec[ii] + w * wVec[ii];
	}
	for (int ii = 0; ii < 3; ii++) {
		xyz[ii] = u
				* (U[ii] + u * (M[ii] + u * A[ii] + v * B[ii] + w * K[ii])
						+ v * w * L[ii])
							+ v * (V[ii] + v * (P[ii] + u * C[ii] + v * E[ii] + w * F[ii]))
							+ w * (W[ii] + w * (R[ii] + u * J[ii] + v * G[ii] + w * H[ii]))
							+ T[ii];
	}
}




