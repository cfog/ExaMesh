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
 * UniformMapping.cxx
 *
 *  Created on: Nov. 12, 2019
 *      Author: cfog
 */

#include "ExaMesh.h"
#include "Mapping.h"

void Q1TetMapping::setupCoordMapping(const emInt verts[]) {
	double coords0[3], coords1[3], coords2[3], coords3[3];
	m_pMesh->getCoords(verts[0], coords0);
	m_pMesh->getCoords(verts[1], coords1);
	m_pMesh->getCoords(verts[2], coords2);
	m_pMesh->getCoords(verts[3], coords3);
	for (int ii = 0; ii < 3; ii++) {
		A[ii] = coords0[ii];
		dU[ii] = coords1[ii] - coords0[ii];
		dV[ii] = coords2[ii] - coords0[ii];
		dW[ii] = coords3[ii] - coords0[ii];
	}
}

void Q1TetMapping::computeTransformedCoords(const double uvw[3],
		double xyz[3]) const {
	const double& u = uvw[0];
	const double& v = uvw[1];
	const double& w = uvw[2];
	for (int ii = 0; ii < 3; ii++) {
		xyz[ii] = A[ii] + u * dU[ii] + v * dV[ii] + w * dW[ii];
	}
}

void Q1PyramidMapping::setupCoordMapping(const emInt verts[]) {
	double coords0[3], coords1[3], coords2[3], coords3[3], coords4[3];
	m_pMesh->getCoords(verts[0], coords0);
	m_pMesh->getCoords(verts[1], coords1);
	m_pMesh->getCoords(verts[2], coords2);
	m_pMesh->getCoords(verts[3], coords3);
	m_pMesh->getCoords(verts[4], coords4);
	for (int ii = 0; ii < 3; ii++) {
		A[ii] = 0.25 * (coords0[ii] + coords1[ii] + coords2[ii] + coords3[ii]);
		dU[ii] = 0.25 * (-coords0[ii] + coords1[ii] + coords2[ii] - coords3[ii]);
		dV[ii] = 0.25 * (-coords0[ii] - coords1[ii] + coords2[ii] + coords3[ii]);
		dUV[ii] = 0.25 * (-coords0[ii] + coords1[ii] - coords2[ii] + coords3[ii]);
		dW[ii] = coords4[ii]
				- 0.25 * (coords0[ii] + coords1[ii] + coords2[ii] + coords3[ii]);
		Apex[ii] = coords4[ii];
	}
}

void Q1PyramidMapping::computeTransformedCoords(const double uvw[3],
		double xyz[3]) const {
	double u = uvw[0];
	double v = uvw[1];
	const double& w = uvw[2];
	if (w == 1) {
		xyz[0] = Apex[0];
		xyz[1] = Apex[1];
		xyz[2] = Apex[2];
	}
	u = 2 * u - (1 - w);
	v = 2 * v - (1 - w);
	for (int ii = 0; ii < 3; ii++) {
		xyz[ii] = (A[ii] + u * dU[ii] + v * dV[ii] + u * v * dUV[ii] / (w - 1))
				+ dW[ii] * w;
	}
}

void Q1PrismMapping::setupCoordMapping(const emInt verts[]) {
	double coords0[3], coords1[3], coords2[3], coords3[3], coords4[3], coords5[3];
	m_pMesh->getCoords(verts[0], coords0);
	m_pMesh->getCoords(verts[1], coords1);
	m_pMesh->getCoords(verts[2], coords2);
	m_pMesh->getCoords(verts[3], coords3);
	m_pMesh->getCoords(verts[4], coords4);
	m_pMesh->getCoords(verts[5], coords5);
	for (int ii = 0; ii < 3; ii++) {
		A[ii] = coords0[ii];
		dU[ii] = coords1[ii] - coords0[ii];
		dV[ii] = coords2[ii] - coords0[ii];

		dW[ii] = coords3[ii] - coords0[ii];
		dUW[ii] = coords4[ii] + coords0[ii] - coords1[ii] - coords3[ii];
		dVW[ii] = coords5[ii] + coords0[ii] - coords2[ii] - coords3[ii];
	}
}

void Q1PrismMapping::computeTransformedCoords(const double uvw[3],
		double xyz[3]) const {
	const double& u = uvw[0];
	const double& v = uvw[1];
	const double& w = uvw[2];
	for (int ii = 0; ii < 3; ii++) {
		xyz[ii] = A[ii] + u * dU[ii] + v * dV[ii]
			+ w * (dW[ii] + u * dUW[ii]	+ v * dVW[ii]);
	}
}

void Q1HexMapping::setupCoordMapping(const emInt verts[]) {
	double coords0[3], coords1[3], coords2[3], coords3[3], coords4[3], coords5[3],
			coords6[3], coords7[3];
	m_pMesh->getCoords(verts[0], coords0);
	m_pMesh->getCoords(verts[1], coords1);
	m_pMesh->getCoords(verts[2], coords2);
	m_pMesh->getCoords(verts[3], coords3);
	m_pMesh->getCoords(verts[4], coords4);
	m_pMesh->getCoords(verts[5], coords5);
	m_pMesh->getCoords(verts[6], coords6);
	m_pMesh->getCoords(verts[7], coords7);
	for (int ii = 0; ii < 3; ii++) {
		// Representing these coefficients as a matrix, wtih the rows (in order)
		// being A, dU, dV, dW, dUV, dUW, dVW, dUVW; and the columns numbered from 0 to 7
		// (left to right), we've got:
		// /   1   									  \.
		// |  -1  1                   |
		// |  -1        1             |
		// |  -1           1          |
		// |   1 -1  1 -1             |
		// |   1 -1       -1  1       |
		// |   1       -1 -1        1 |
		// \  -1  1 -1  1  1 -1  1 -1 /
		//
		// Those appear to be correct when you have more than one of u, v, w set to
		// one: you get the coords for the correct corner.
		A[ii] = coords0[ii];
		dU[ii] = coords1[ii] - coords0[ii];
		dV[ii] = coords3[ii] - coords0[ii];
		dW[ii] = coords4[ii] - coords0[ii];

		dUV[ii] = coords2[ii] + coords0[ii] - coords1[ii] - coords3[ii];
		dUW[ii] = coords5[ii] + coords0[ii] - coords1[ii] - coords4[ii];
		dVW[ii] = coords7[ii] + coords0[ii] - coords3[ii] - coords4[ii];

		dUVW[ii] = coords6[ii] - coords7[ii] - coords5[ii] - coords0[ii]
				+ coords4[ii] + coords3[ii]
								- coords2[ii]
								+ coords1[ii];
	}
}

void Q1HexMapping::computeTransformedCoords(const double uvw[3],
		double xyz[3]) const {
	const double& u = uvw[0];
	const double& v = uvw[1];
	const double& w = uvw[2];
	for (int ii = 0; ii < 3; ii++) {
		// Original, naive implementation of polynomial evaluation
		xyz[ii] = A[ii] + u * dU[ii] + v * dV[ii] + w * dW[ii] + u * v * dUV[ii]
							+ u * w * dUW[ii] + v * w * dVW[ii] + u * v * w * dUVW[ii];
		// Faster version that has only one multiply-add per term.
//		xyz[ii] = A[ii] + u * (dU[ii] + v * dUV[ii] + w * (dUW[ii] + v * dUVW[ii]))
//								 + v * dV[ii] + w * (dW[ii] + v * dVW[ii]);
	}
}


