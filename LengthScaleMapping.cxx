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
 * LengthScaleMapping.cxx
 *
 *  Created on: Oct. 1, 2019
 *      Author: cfog
 */

#include "ExaMesh.h"
#include <algorithm>

#include "Mapping.h"

double Mapping::getIsoLengthScale(const emInt vertInd) {
	return m_pMesh->getLengthScale(vertInd);
}

void LengthScaleTetMapping::setPolyCoeffs(const double xyz[4][3],
		double uderiv[4][3], double vderiv[4][3], double wderiv[4][3]) {
	const double *xyz0 = xyz[0];
	const double *xyz1 = xyz[1];
	const double *xyz2 = xyz[2];
	const double *xyz3 = xyz[3];

	const double *uderiv0 = uderiv[0];
	const double *uderiv1 = uderiv[1];
	const double *uderiv2 = uderiv[2];
	const double *uderiv3 = uderiv[3];

	const double *vderiv0 = vderiv[0];
	const double *vderiv1 = vderiv[1];
	const double *vderiv2 = vderiv[2];
	const double *vderiv3 = vderiv[3];

	const double *wderiv0 = wderiv[0];
	const double *wderiv1 = wderiv[1];
	const double *wderiv2 = wderiv[2];
	const double *wderiv3 = wderiv[3];

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
		L[ii] = (-6 * xyz0[ii] + 2 * xyz1[ii] + 2 * xyz2[ii] + 2 * xyz3[ii])
				+ (-2 * uderiv0[ii] - uderiv1[ii] + 0.5 * uderiv2[ii]
						+ 0.5 * uderiv3[ii])
				+ (-2 * vderiv0[ii] + 0.5 * vderiv1[ii] - vderiv2[ii]
						+ 0.5 * vderiv3[ii])
				+ (-2 * wderiv0[ii] + 0.5 * wderiv1[ii] + 0.5 * wderiv2[ii]
						- wderiv3[ii]);
		M[ii] = -uderiv1[ii] + 3 * xyz1[ii] - 3 * xyz0[ii] - 2 * uderiv0[ii];
		P[ii] = -vderiv2[ii] + 3 * xyz2[ii] - 3 * xyz0[ii] - 2 * vderiv0[ii];
		R[ii] = -wderiv3[ii] + 3 * xyz3[ii] - 3 * xyz0[ii] - 2 * wderiv0[ii];
		T[ii] = xyz0[ii];
		U[ii] = uderiv0[ii];
		V[ii] = vderiv0[ii];
		W[ii] = wderiv0[ii];
	}
}

void LengthScaleTetMapping::setupCoordMapping(const emInt verts[]) {
	double xyz[4][3];
	m_pMesh->getCoords(verts[0], xyz[0]);
	m_pMesh->getCoords(verts[1], xyz[1]);
	m_pMesh->getCoords(verts[2], xyz[2]);
	m_pMesh->getCoords(verts[3], xyz[3]);

	double len0 = getIsoLengthScale(verts[0]);
	double len1 = getIsoLengthScale(verts[1]);
	double len2 = getIsoLengthScale(verts[2]);
	double len3 = getIsoLengthScale(verts[3]);

	double vec01[] = DIFF(xyz[1], xyz[0]);
	double vec02[] = DIFF(xyz[2], xyz[0]);
	double vec03[] = DIFF(xyz[3], xyz[0]);
	double vec12[] = DIFF(xyz[2], xyz[1]);
	double vec13[] = DIFF(xyz[3], xyz[1]);
	double vec23[] = DIFF(xyz[3], xyz[2]);

	double ratio01 = sqrt(std::max(0.5, std::min(2., len0/len1)));
	double ratio02 = sqrt(std::max(0.5, std::min(2., len0/len2)));
	double ratio03 = sqrt(std::max(0.5, std::min(2., len0/len3)));
	double ratio12 = sqrt(std::max(0.5, std::min(2., len1/len2)));
	double ratio13 = sqrt(std::max(0.5, std::min(2., len1/len3)));
	double ratio23 = sqrt(std::max(0.5, std::min(2., len2/len3)));

	double uderiv[4][3], vderiv[4][3], wderiv[4][3];

	SCALE(vec01, ratio01, uderiv[0]);
	SCALE(vec02, ratio02, vderiv[0]);
	SCALE(vec03, ratio03, wderiv[0]);

	SCALE(vec01, 1 / ratio01, uderiv[1]);
	for (int ii = 0; ii < 3; ii++) {
		vderiv[1][ii] = uderiv[1][ii] + vec12[ii] * ratio12;
		wderiv[1][ii] = uderiv[1][ii] + vec13[ii] * ratio13;
	}

	SCALE(vec02, 1 / ratio02, vderiv[2]);
	for (int ii = 0; ii < 3; ii++) {
		uderiv[2][ii] = vderiv[2][ii] - vec12[ii] / ratio12;
		wderiv[2][ii] = vderiv[2][ii] + vec23[ii] * ratio23;
	}

	SCALE(vec03, 1 / ratio03, wderiv[3]);
	for (int ii = 0; ii < 3; ii++) {
		uderiv[3][ii] = wderiv[3][ii] - vec13[ii] / ratio13;
		vderiv[3][ii] = wderiv[3][ii] - vec23[ii] / ratio23;
	}

	setPolyCoeffs(xyz, uderiv, vderiv, wderiv);
}

void LengthScaleTetMapping::computeTransformedCoords(const double uvw[3],
		double xyz[3]) const {
	const double &u = uvw[0];
	const double &v = uvw[1];
	const double &w = uvw[2];
	for (int ii = 0; ii < 3; ii++) {
		xyz[ii] = u
				* (U[ii] + u * (M[ii] + u * A[ii] + v * B[ii] + w * K[ii])
						+ v * w * L[ii])
				+ v * (V[ii] + v * (P[ii] + u * C[ii] + v * E[ii] + w * F[ii]))
				+ w * (W[ii] + w * (R[ii] + u * J[ii] + v * G[ii] + w * H[ii]))
				+ T[ii];
	}
}

void LengthScalePyramidMapping::setPolyCoeffs(const double xyz[5][3],
		const double uderiv[5][3], const double vderiv[5][3],
		const double wderiv[5][3]) {
	// Define a bunch of aliases to make things slightly for readable.
	// fabc: a is the u coord, b is the v coord, c is the w coord.
	// Verts are numbered counterclockwise around the bottom of the hex (starting at
	// u = v = w = 0).
	const double *f000 = xyz[0];
	const double *f100 = xyz[1];
	const double *f110 = xyz[2];
	const double *f010 = xyz[3];
	const double *f001 = xyz[4];

	const double *fu000 = uderiv[0];
	const double *fu100 = uderiv[1];
	const double *fu110 = uderiv[2];
	const double *fu010 = uderiv[3];
	const double *fu001 = uderiv[4];

	const double *fv000 = vderiv[0];
	const double *fv100 = vderiv[1];
	const double *fv110 = vderiv[2];
	const double *fv010 = vderiv[3];
	const double *fv001 = vderiv[4];

	const double *fw000 = wderiv[0];
	const double *fw100 = wderiv[1];
	const double *fw110 = wderiv[2];
	const double *fw010 = wderiv[3];
	const double *fw001 = wderiv[4];

	for (int ii = 0; ii < 3; ii++) {
		A[ii] = f000[ii];
		B[ii] = fu000[ii];
		C[ii] = fv000[ii];
		E[ii] = -fu100[ii] + 3 * f100[ii] - 2 * fu000[ii] - 3 * f000[ii];
		F[ii] = fu010[ii] - fu000[ii] - f110[ii] + f100[ii] - f000[ii]
				+ f010[ii] + fv100[ii] - fv000[ii];
		G[ii] = -fv010[ii] + 3 * f010[ii] - 2 * fv000[ii] - 3 * f000[ii];
		H[ii] = -2 * f100[ii] + fu000[ii] + 2 * f000[ii] + fu100[ii];
		J[ii] = -fu110[ii] - 2 * fu010[ii] + 2 * fu000[ii] + 3 * f110[ii]
				- 3 * f100[ii] + 3 * f000[ii] - 3 * f010[ii] + fu100[ii];
		K[ii] = -fv110[ii] + 3 * f110[ii] - 3 * f100[ii] + 3 * f000[ii]
				- 3 * f010[ii] - 2 * fv100[ii] + 2 * fv000[ii] + fv010[ii];
		L[ii] = -2 * f010[ii] + fv000[ii] + 2 * f000[ii] + fv010[ii];
		M[ii] = fu010[ii] - fu000[ii] - 2 * f110[ii] + 2 * f100[ii]
				- 2 * f000[ii] + 2 * f010[ii] + fu110[ii] - fu100[ii];
		N[ii] = -2 * f110[ii] + 2 * f100[ii] - 2 * f000[ii] + 2 * f010[ii]
				+ fv100[ii] - fv000[ii] + fv110[ii] - fv010[ii];
		P[ii] = f001[ii];
		Q[ii] = fw000[ii];
		R[ii] = fw100[ii] - fw000[ii];
		S[ii] = fw010[ii] - fw000[ii];
		T[ii] = fw110[ii] + fw000[ii] - fw100[ii] - fw010[ii];
		U[ii] = fu001[ii];
		V[ii] = fv001[ii];
		W[ii] = fw001[ii];
	}
}

void LengthScalePyramidMapping::setupCoordMapping(const emInt verts[]) {
	double xyz[5][3], vertLen[5];
	for (int ii = 0; ii < 5; ii++) {
		m_pMesh->getCoords(verts[ii], xyz[ii]);
		vertLen[ii] = getIsoLengthScale(verts[ii]);
	}

	int edgeVerts[8][2] = { { 0, 1 }, { 1, 2 }, { 3, 2 }, { 0, 3 }, { 0, 4 }, {
			1, 4 }, { 2, 4 }, { 3, 4 } };
	double edgeVecs[8][3], ratios[8];
	for (int ii = 0; ii < 8; ii++) {
		int startVert = edgeVerts[ii][0];
		int endVert = edgeVerts[ii][1];
		edgeVecs[ii][0] = xyz[endVert][0] - xyz[startVert][0];
		edgeVecs[ii][1] = xyz[endVert][1] - xyz[startVert][1];
		edgeVecs[ii][2] = xyz[endVert][2] - xyz[startVert][2];
		ratios[ii] = sqrt(
				std::max(0.5,
						std::min(2., vertLen[startVert] / vertLen[endVert])));
	}
	double uderiv[5][3], vderiv[5][3], wderiv[5][3];

	// Set u derivatives for the base verts.
	SCALE(edgeVecs[0], ratios[0], uderiv[edgeVerts[0][0]]);
	SCALE(edgeVecs[0], 1. / ratios[0], uderiv[edgeVerts[0][1]]);
	SCALE(edgeVecs[2], ratios[2], uderiv[edgeVerts[2][0]]);
	SCALE(edgeVecs[2], 1. / ratios[2], uderiv[edgeVerts[2][1]]);

	// Set v derivatives for the base verts.
	SCALE(edgeVecs[1], ratios[1], vderiv[edgeVerts[1][0]]);
	SCALE(edgeVecs[1], 1. / ratios[1], vderiv[edgeVerts[1][1]]);
	SCALE(edgeVecs[3], ratios[3], vderiv[edgeVerts[3][0]]);
	SCALE(edgeVecs[3], 1. / ratios[3], vderiv[edgeVerts[3][1]]);

	// Set w derivatives for the base verts.
	SCALE(edgeVecs[4], ratios[4], wderiv[edgeVerts[4][0]]);
	SCALE(edgeVecs[5], ratios[5], wderiv[edgeVerts[5][0]]);
	SCALE(edgeVecs[6], ratios[6], wderiv[edgeVerts[6][0]]);
	SCALE(edgeVecs[7], ratios[7], wderiv[edgeVerts[7][0]]);

	wderiv[4][0] = wderiv[4][1] = wderiv[4][2] = 0;
	for (int ii = 4; ii <= 7; ii++) {
		double temp[3];
		SCALE(edgeVecs[ii], 1. / ratios[ii], temp);
		wderiv[4][0] += temp[0] / 4;
		wderiv[4][1] += temp[1] / 4;
		wderiv[4][2] += temp[2] / 4;
	}

	setPolyCoeffs(xyz, uderiv, vderiv, wderiv);
}

void LengthScalePyramidMapping::computeTransformedCoords(const double uvw[3],
		double xyz[3]) const {
	const double &u = uvw[0];
	const double &v = uvw[1];
	const double &w = uvw[2];
	double baseFn[3], baseDeriv[3];
	const double H1 = 2 * w * w * w - 3 * w * w + 1;
	const double H2 = 1 - H1;
	const double H3 = w * (w - 1) * (w - 1);
	const double H4 = w * w * (w - 1);
	for (int ii = 0; ii < 3; ii++) {
		baseFn[ii] = A[ii]
				+ u
						* (B[ii]
								+ u
										* (E[ii] + u * (H[ii] + v * M[ii])
												+ v * J[ii]))
				+ v
						* (C[ii] + u * F[ii]
								+ v
										* (G[ii] + u * K[ii]
												+ v * (L[ii] + u * N[ii])));
		baseDeriv[ii] = Q[ii] + u * R[ii] + v * (S[ii] + u * T[ii]);

		xyz[ii] = baseFn[ii] * H1 + P[ii] * H2 + baseDeriv[ii] * H3
				+ W[ii] * H4;
		// For now, don't code the u and v derivative terms at the top; those
		// seem hard to define properly anyway.
	}
}

void LengthScaleHexMapping::setPolyCoeffs(const double xyz[8][3],
		double uderiv[8][3], double vderiv[8][3], double wderiv[8][3]) {
	// Define a bunch of aliases to make things slightly for readable.
	// fabc: a is the u coord, b is the v coord, c is the w coord.
	// Verts are numbered counterclockwise around the bottom of the hex (starting at
	// u = v = w = 0), then CCW around the top, with 4 above 0.
	const double *f000 = xyz[0];
	const double *f100 = xyz[1];
	const double *f110 = xyz[2];
	const double *f010 = xyz[3];
	const double *f001 = xyz[4];
	const double *f101 = xyz[5];
	const double *f111 = xyz[6];
	const double *f011 = xyz[7];

	const double *fu000 = uderiv[0];
	const double *fu100 = uderiv[1];
	const double *fu110 = uderiv[2];
	const double *fu010 = uderiv[3];
	const double *fu001 = uderiv[4];
	const double *fu101 = uderiv[5];
	const double *fu111 = uderiv[6];
	const double *fu011 = uderiv[7];

	const double *fv000 = vderiv[0];
	const double *fv100 = vderiv[1];
	const double *fv110 = vderiv[2];
	const double *fv010 = vderiv[3];
	const double *fv001 = vderiv[4];
	const double *fv101 = vderiv[5];
	const double *fv111 = vderiv[6];
	const double *fv011 = vderiv[7];

	const double *fw000 = wderiv[0];
	const double *fw100 = wderiv[1];
	const double *fw110 = wderiv[2];
	const double *fw010 = wderiv[3];
	const double *fw001 = wderiv[4];
	const double *fw101 = wderiv[5];
	const double *fw111 = wderiv[6];
	const double *fw011 = wderiv[7];

	for (int ii = 0; ii < 3; ii++) {
		A[ii] = -2 * f100[ii] + 2 * f000[ii] + fu000[ii] + fu100[ii];
		B[ii] = 3 * f110[ii] - 3 * f100[ii] + 3 * f000[ii] + 2 * fu000[ii]
				- 2 * fu010[ii] - 3 * f010[ii] - fu110[ii] + fu100[ii];
		C[ii] = 3 * f110[ii] - 3 * f100[ii] + 3 * f000[ii] - 2 * fv100[ii]
				+ 2 * fv000[ii] - 3 * f010[ii] - fv110[ii] + fv010[ii];
		E[ii] = -2 * f010[ii] + 2 * f000[ii] + fv000[ii] + fv010[ii];
		F[ii] = -fv011[ii] + 3 * f011[ii] + 2 * fv000[ii] - 3 * f001[ii]
				- 3 * f010[ii] + 3 * f000[ii] - 2 * fv001[ii] + fv010[ii];
		G[ii] = -fw011[ii] - 2 * fw010[ii] + 2 * fw000[ii] + 3 * f011[ii]
				- 3 * f001[ii] - 3 * f010[ii] + 3 * f000[ii] + fw001[ii];
		H[ii] = -2 * f001[ii] + 2 * f000[ii] + fw000[ii] + fw001[ii];
		J[ii] = -fw101[ii] + 3 * f101[ii] + 2 * fw000[ii] - 3 * f001[ii]
				- 3 * f100[ii] + 3 * f000[ii] - 2 * fw100[ii] + fw001[ii];
		K[ii] = -fu101[ii] - 2 * fu001[ii] + 2 * fu000[ii] + 3 * f101[ii]
				- 3 * f001[ii] - 3 * f100[ii] + 3 * f000[ii] + fu100[ii];
		L[ii] = fw110[ii] - fw100[ii] - fv100[ii] - 2 * f100[ii] + fw000[ii]
				+ fv000[ii] + fu000[ii] + 2 * f000[ii] + 2 * f110[ii]
				- 2 * f111[ii] - fv001[ii] - fu001[ii] - 2 * f001[ii]
				- fw010[ii] - fu010[ii] - 2 * f010[ii] + 2 * f011[ii]
				+ fu011[ii] + 2 * f101[ii] + fv101[ii];
		M[ii] = -fu100[ii] + 3 * f100[ii] - 3 * f000[ii] - 2 * fu000[ii];
		N[ii] = -f110[ii] + f100[ii] - f000[ii] - fu000[ii] + fv100[ii]
				- fv000[ii] + fu010[ii] + f010[ii];
		P[ii] = -fv010[ii] + 3 * f010[ii] - 3 * f000[ii] - 2 * fv000[ii];
		Q[ii] = -f011[ii] - fv000[ii] - fw000[ii] + f001[ii] + f010[ii]
				- f000[ii] + fv001[ii] + fw010[ii];
		R[ii] = -fw001[ii] + 3 * f001[ii] - 3 * f000[ii] - 2 * fw000[ii];
		S[ii] = -f101[ii] - fu000[ii] - fw000[ii] + f001[ii] + f100[ii]
				- f000[ii] + fw100[ii] + fu001[ii];
		T[ii] = f000[ii];
		U[ii] = fu000[ii];
		V[ii] = fv000[ii];
		W[ii] = fw000[ii];

		AA[ii] = fu110[ii] - 2 * f110[ii] + 2 * f100[ii] - 2 * f000[ii]
				- fu000[ii] + fu010[ii] + 2 * f010[ii] - fu100[ii];
		BB[ii] = fv110[ii] + fv100[ii] - 2 * f110[ii] + 2 * f100[ii]
				- 2 * f000[ii] - fv000[ii] + 2 * f010[ii] - fv010[ii];
		CC[ii] = -2 * f011[ii] - fv000[ii] + 2 * f001[ii] + 2 * f010[ii]
				- 2 * f000[ii] + fv001[ii] + fv011[ii] - fv010[ii];
		DD[ii] = -2 * f011[ii] - fw000[ii] + 2 * f001[ii] + 2 * f010[ii]
				- 2 * f000[ii] + fw010[ii] + fw011[ii] - fw001[ii];
		EE[ii] = -2 * f101[ii] - fw000[ii] + 2 * f001[ii] + 2 * f100[ii]
				- 2 * f000[ii] + fw100[ii] + fw101[ii] - fw001[ii];
		FF[ii] = -2 * f101[ii] - fu000[ii] + 2 * f001[ii] + 2 * f100[ii]
				- 2 * f000[ii] + fu001[ii] + fu101[ii] - fu100[ii];
		GG[ii] = -fu100[ii] + 3 * f100[ii] - 2 * fu000[ii] - 3 * f000[ii]
				- 3 * f110[ii] + fu110[ii] + 3 * f111[ii] - fu111[ii]
				+ 2 * fu001[ii] + 3 * f001[ii] + 2 * fu010[ii] + 3 * f010[ii]
				- 3 * f011[ii] - 2 * fu011[ii] - 3 * f101[ii] + fu101[ii];
		HH[ii] = 2 * fv100[ii] + 3 * f100[ii] + fv110[ii] - 2 * fv000[ii]
				- 3 * f000[ii] - 3 * f110[ii] + 3 * f111[ii] - fv111[ii]
				+ 2 * fv001[ii] + 3 * f001[ii] - fv010[ii] + 3 * f010[ii]
				- 3 * f011[ii] + fv011[ii] - 3 * f101[ii] - 2 * fv101[ii];
		II[ii] = -2 * fw110[ii] + 2 * fw100[ii] + 3 * f100[ii] - 2 * fw000[ii]
				- 3 * f000[ii] - 3 * f110[ii] + fw011[ii] + 3 * f111[ii]
				+ fw101[ii] - fw111[ii] - fw001[ii] + 3 * f001[ii]
				+ 2 * fw010[ii] + 3 * f010[ii] - 3 * f011[ii] - 3 * f101[ii];
		JJ[ii] = fu100[ii] - 2 * f100[ii] + fu000[ii] + 2 * f000[ii]
				+ 2 * f110[ii] - fu110[ii] - 2 * f111[ii] + fu111[ii]
				- fu001[ii] - 2 * f001[ii] - fu010[ii] - 2 * f010[ii]
				+ 2 * f011[ii] + fu011[ii] + 2 * f101[ii] - fu101[ii];
		KK[ii] = -fv100[ii] - 2 * f100[ii] - fv110[ii] + fv000[ii]
				+ 2 * f000[ii] + 2 * f110[ii] - 2 * f111[ii] + fv111[ii]
				- fv001[ii] - 2 * f001[ii] + fv010[ii] - 2 * f010[ii]
				+ 2 * f011[ii] - fv011[ii] + 2 * f101[ii] + fv101[ii];
		LL[ii] = fw110[ii] - fw100[ii] - 2 * f100[ii] + fw000[ii] + 2 * f000[ii]
				+ 2 * f110[ii] - fw011[ii] - 2 * f111[ii] - fw101[ii]
				+ fw111[ii] + fw001[ii] - 2 * f001[ii] - fw010[ii]
				- 2 * f010[ii] + 2 * f011[ii] + 2 * f101[ii];
	}
}

void LengthScaleHexMapping::setupCoordMapping(const emInt verts[]) {
	double xyz[8][3], vertLen[8];
	for (int ii = 0; ii < 8; ii++) {
		m_pMesh->getCoords(verts[ii], xyz[ii]);
		vertLen[ii] = getIsoLengthScale(verts[ii]);
	}

	// Edges are ordered in this specific way to make it easier to set
	// derivatives later.
	int edgeVerts[12][2] = { { 0, 1 }, { 2, 3 }, { 4, 5 }, { 6, 7 }, { 1, 2 }, {
			3, 0 }, { 5, 6 }, { 7, 4 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
	double edgeVecs[12][3], ratios[12];
	for (int ii = 0; ii < 12; ii++) {
		int startVert = edgeVerts[ii][0];
		int endVert = edgeVerts[ii][1];
		edgeVecs[ii][0] = xyz[startVert][0] - xyz[endVert][0];
		edgeVecs[ii][1] = xyz[startVert][1] - xyz[endVert][1];
		edgeVecs[ii][2] = xyz[startVert][2] - xyz[endVert][2];
		ratios[ii] = sqrt(
				std::max(0.5,
						std::min(2., vertLen[startVert] / vertLen[endVert])));
	}
	double uderiv[8][3], vderiv[8][3], wderiv[8][3];

	for (int ii = 0; ii < 4; ii++) {
		// Set u derivatives for all eight verts.
		SCALE(edgeVecs[ii], ratios[ii], uderiv[edgeVerts[ii][0]]);
		SCALE(edgeVecs[ii], 1. / ratios[ii], uderiv[edgeVerts[ii][1]]);

		// Set v derivatives for all eight verts.
		SCALE(edgeVecs[ii + 4], ratios[ii + 4], vderiv[edgeVerts[ii + 4][0]]);
		SCALE(edgeVecs[ii + 4], 1. / ratios[ii + 4],
				vderiv[edgeVerts[ii + 4][1]]);

		// Set w derivatives for all eight verts.
		SCALE(edgeVecs[ii + 8], ratios[ii + 8], wderiv[edgeVerts[ii + 8][0]]);
		SCALE(edgeVecs[ii + 8], 1. / ratios[ii + 8],
				wderiv[edgeVerts[ii + 8][1]]);
	}

	setPolyCoeffs(xyz, uderiv, vderiv, wderiv);
}

void LengthScaleHexMapping::computeTransformedCoords(const double uvw[3],
		double xyz[3]) const {
	const double &u = uvw[0];
	const double &v = uvw[1];
	const double &w = uvw[2];
	for (int ii = 0; ii < 3; ii++) {
		xyz[ii] =
				u
						* (U[ii]
								+ u
										* (M[ii] + u * A[ii]
												+ v
														* (B[ii] + u * AA[ii]
																+ w
																		* (GG[ii]
																				+ u
																						* JJ[ii]))
												+ w * (K[ii] + u * FF[ii]))
								+ v * (N[ii] + w * L[ii]))
						+ v
								* (V[ii]
										+ v
												* (P[ii]
														+ u
																* (C[ii]
																		+ v
																				* BB[ii])
														+ v * E[ii]
														+ w
																* (F[ii]
																		+ v
																				* CC[ii]
																		+ u
																				* (HH[ii]
																						+ v
																								* KK[ii])))
										+ w * Q[ii])
						+ w
								* (W[ii] + u * S[ii]
										+ w
												* (R[ii]
														+ u
																* (J[ii]
																		+ w
																				* EE[ii]
																		+ v
																				* (II[ii]
																						+ w
																								* LL[ii]))
														+ v
																* (G[ii]
																		+ w
																				* DD[ii])
														+ w * H[ii])) + T[ii];
	}
}

void LengthScalePrismMapping::computeTransformedCoords(
		const double uvwCoords[3], double xyz[3]) const {
	const double &u = uvwCoords[0];
	const double &v = uvwCoords[1];
	const double &w = uvwCoords[2];

	const double u2 = u * u;
	const double uv = u * v;
	const double v2 = v * v;

	const double u3 = u2 * u;
	const double u2v = u2 * v;
	const double uv2 = u * v2;
	const double v3 = v2 * v;

	double Hbot = 1 + w * w * (2 * w - 3);
	double Htop = 1 - Hbot;
	double HderivBot = w * (w - 1) * (w - 1);
	double HderivTop = w * w * (w - 1);

	for (int ii = 0; ii < 3; ii++) {
		double bot = Cb[ii] + Cbu[ii] * u + Cbv[ii] * v + Cbu2[ii] * u2
				+ Cbuv[ii] * uv + Cbv2[ii] * v2 + Cbu3[ii] * u3 + Cbu2v[ii] * u2v
				+ Cbuv2[ii] * uv2 + Cbv3[ii] * v3;
		double top = Ct[ii] + Ctu[ii] * u + Ctv[ii] * v + Ctu2[ii] * u2
				+ Ctuv[ii] * uv + Ctv2[ii] * v2 + Ctu3[ii] * u3 + Ctu2v[ii] * u2v
				+ Ctuv2[ii] * uv2 + Ctv3[ii] * v3;

		double botGrad = Gb[ii] + Gbu[ii] * u + Gbv[ii] * v;
		double topGrad = Gt[ii] + Gtu[ii] * u + Gtv[ii] * v;

		xyz[ii] = bot * Hbot + top * Htop + botGrad * HderivBot
				+ topGrad * HderivTop;
	}
}

void LengthScalePrismMapping::setPolyCoeffs(const double xyz[6][3],
		double uderiv[6][3], double vderiv[6][3], double wderiv[6][3]) {
	const double *f000 = xyz[0];
	const double *f100 = xyz[1];
	const double *f010 = xyz[2];
	const double *f001 = xyz[3];
	const double *f101 = xyz[4];
	const double *f011 = xyz[5];

	const double *fu000 = uderiv[0];
	const double *fu100 = uderiv[1];
	const double *fu010 = uderiv[2];
	const double *fu001 = uderiv[3];
	const double *fu101 = uderiv[4];
	const double *fu011 = uderiv[5];

	const double *fv000 = vderiv[0];
	const double *fv100 = vderiv[1];
	const double *fv010 = vderiv[2];
	const double *fv001 = vderiv[3];
	const double *fv101 = vderiv[4];
	const double *fv011 = vderiv[5];

	const double *fw000 = wderiv[0];
	const double *fw100 = wderiv[1];
	const double *fw010 = wderiv[2];
	const double *fw001 = wderiv[3];
	const double *fw101 = wderiv[4];
	const double *fw011 = wderiv[5];

	for (int ii = 0; ii < 3; ii++) {
		Cbu3[ii] = -2 * f100[ii] + 2 * f000[ii] + fu000[ii] + fu100[ii];
		Cbu2v[ii] =  fv100[ii] - fv000[ii]; /* Wrong? */  // 6*f000[ii] + 2*fu000[ii] + fv000[ii] + fv100[ii]; //
		Cbuv2[ii] =  fu010[ii] - fu000[ii]; /* Wrong? */  // 6*f000[ii] + 2*fv000[ii] + fu000[ii] + fu010[ii]; //
		Cbv3[ii] = -2 * f010[ii] + fv000[ii] + 2 * f000[ii] + fv010[ii];

		Cbu2[ii] = -fu100[ii] + 3 * f100[ii] - 2 * fu000[ii] - 3 * f000[ii];
		Cbuv[ii] = 0; //-6*f000[ii] - 2 * fu000[ii] - 2 * fv000[ii];
		Cbv2[ii] = -fv010[ii] + 3 * f010[ii] - 2 * fv000[ii] - 3 * f000[ii];

		Cbu[ii] = fu000[ii];
		Cbv[ii] = fv000[ii];

		Cb[ii] = f000[ii];

		Ctu3[ii] = -2 * f101[ii] + 2 * f001[ii] + fu001[ii] + fu101[ii];
		Ctu2v[ii] = fv101[ii] - fv001[ii]; // 6*f001[ii] + 2*fu001[ii] + fv001[ii] + fv101[ii]; //
		Ctuv2[ii] = fu011[ii] - fu001[ii]; // 6*f001[ii] + 2*fv001[ii] + fu001[ii] + fu011[ii]; //
		Ctv3[ii] = -2 * f011[ii] + fv001[ii] + 2 * f001[ii] + fv011[ii];

		Ctu2[ii] = -fu101[ii] + 3 * f101[ii] - 2 * fu001[ii] - 3 * f001[ii];
		Ctuv[ii] = 0; // -6*f001[ii] - 2 * fu001[ii] - 2 * fv001[ii];
		Ctv2[ii] = -fv011[ii] + 3 * f011[ii] - 2 * fv001[ii] - 3 * f001[ii];

		Ctu[ii] = fu001[ii];
		Ctv[ii] = fv001[ii];

		Ct[ii] = f001[ii];

		Gb[ii] = fw000[ii];
		Gbu[ii] = fw100[ii] - fw000[ii];
		Gbv[ii] = fw010[ii] - fw000[ii];

		Gt[ii] = fw001[ii];
		Gtu[ii] = fw101[ii] - fw001[ii];
		Gtv[ii] = fw011[ii] - fw001[ii];
	}
}

void LengthScalePrismMapping::setupCoordMapping(const emInt verts[]) {
	double xyz[6][3], len[6];
	for (int ii = 0; ii < 6; ii++) {
		m_pMesh->getCoords(verts[ii], xyz[ii]);
		len[ii] = getIsoLengthScale(verts[ii]);
	}
	// All vectors are written in order of increasing u/v/w, except the
	// diagonal (12 and 45) edges.
	double vec01[] = DIFF(xyz[1], xyz[0]);
	double vec02[] = DIFF(xyz[2], xyz[0]);
	double vec12[] = DIFF(xyz[2], xyz[1]);

	double vec34[] = DIFF(xyz[4], xyz[3]);
	double vec35[] = DIFF(xyz[5], xyz[3]);
	double vec45[] = DIFF(xyz[5], xyz[4]);

	double vec03[] = DIFF(xyz[3], xyz[0]);
	double vec14[] = DIFF(xyz[4], xyz[1]);
	double vec25[] = DIFF(xyz[5], xyz[2]);

	double ratio01 = sqrt(std::max(0.5, std::min(2., len[0]/len[1])));
	double ratio12 = sqrt(std::max(0.5, std::min(2., len[1]/len[2])));
	double ratio02 = sqrt(std::max(0.5, std::min(2., len[0]/len[2])));

	double ratio34 = sqrt(std::max(0.5, std::min(2., len[3]/len[4])));
	double ratio45 = sqrt(std::max(0.5, std::min(2., len[4]/len[5])));
	double ratio35 = sqrt(std::max(0.5, std::min(2., len[3]/len[5])));

	double ratio03 = sqrt(std::max(0.5, std::min(2., len[0]/len[3])));
	double ratio14 = sqrt(std::max(0.5, std::min(2., len[1]/len[4])));
	double ratio25 = sqrt(std::max(0.5, std::min(2., len[2]/len[5])));

	double uderiv[6][3],
	vderiv[6][3], wderiv[6][3];

	// The bottom corners
	SCALE(vec01, ratio01, uderiv[0]);
	SCALE(vec02, ratio02, vderiv[0]);
	SCALE(vec03, ratio03, wderiv[0]);

	SCALE(vec01, 1 / ratio01, uderiv[1]);
	SCALE(vec14, ratio14, wderiv[1]);
	for (int ii = 0; ii < 3; ii++) {
		vderiv[1][ii] = uderiv[1][ii] + vec12[ii] * ratio12;
	}

	SCALE(vec02, 1 / ratio02, vderiv[2]);
	SCALE(vec25, ratio25, wderiv[2]);
	for (int ii = 0; ii < 3; ii++) {
		uderiv[2][ii] = vderiv[2][ii] - vec12[ii] / ratio12;
	}

	// The top corners
	SCALE(vec34, ratio34, uderiv[3]);
	SCALE(vec35, ratio35, vderiv[3]);
	SCALE(vec03, 1 / ratio03, wderiv[3]);

	SCALE(vec34, 1 / ratio34, uderiv[4]);
	SCALE(vec14, 1 / ratio14, wderiv[4]);
	for (int ii = 0; ii < 3; ii++) {
		vderiv[4][ii] = uderiv[4][ii] + vec45[ii] * ratio45;
	}

	SCALE(vec35, 1 / ratio35, vderiv[5]);
	SCALE(vec25, 1 / ratio25, wderiv[5]);
	for (int ii = 0; ii < 3; ii++) {
		uderiv[5][ii] = vderiv[5][ii] - vec45[ii] / ratio45;
	}

	setPolyCoeffs(xyz, uderiv, vderiv, wderiv);
}

