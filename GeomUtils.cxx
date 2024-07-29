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
 * GeomUtils.cxx
 *
 *  Created on: Jul. 3, 2019
 *      Author: cfog
 */

#include "exa-defs.h"

#define MAX(a, b) \
	((a) > (b) ? (a) : (b))

#define MINMAX4(da, db, dc, dd, dMin, dMax)                                    \
  {                                                                            \
    double dMax1, dMax2, dMin1, dMin2;                                         \
    if (da > db) {                                                             \
      dMax1 = da;                                                              \
      dMin1 = db;                                                              \
    } else {                                                                   \
      dMax1 = db;                                                              \
      dMin1 = da;                                                              \
    }                                                                          \
    if (dc > dd) {                                                             \
      dMax2 = dc;                                                              \
      dMin2 = dd;                                                              \
    } else {                                                                   \
      dMax2 = dd;                                                              \
      dMin2 = dc;                                                              \
    }                                                                          \
    if (dMax1 > dMax2)                                                         \
      dMax = dMax1;                                                            \
    else                                                                       \
      dMax = dMax2;                                                            \
    if (dMin1 < dMin2)                                                         \
      dMin = dMin1;                                                            \
    else                                                                       \
      dMin = dMin2;                                                            \
  }

int checkOrient3D(const double adA[3], const double adB[3], const double adC[3],
		const double adD[3])
// Computes the orientation of four verts in 3D, using as nearly
// exact arithmetic as required.
		{

	double dXa = adA[0];
	double dYa = adA[1];
	double dZa = adA[2];

	double dXb = adB[0];
	double dYb = adB[1];
	double dZb = adB[2];

	double dXc = adC[0];
	double dYc = adC[1];
	double dZc = adC[2];

	double dXd = adD[0];
	double dYd = adD[1];
	double dZd = adD[2];

	double dMaxX, dMinX, dMaxY, dMinY, dMaxZ, dMinZ;

	double dDX2 = dXb - dXa;
	double dDX3 = dXc - dXa;
	double dDX4 = dXd - dXa;
	MINMAX4(dXa, dXb, dXc, dXd, dMinX, dMaxX);

	double dDY2 = dYb - dYa;
	double dDY3 = dYc - dYa;
	double dDY4 = dYd - dYa;
	MINMAX4(dYa, dYb, dYc, dYd, dMinY, dMaxY);

	double dDZ2 = dZb - dZa;
	double dDZ3 = dZc - dZa;
	double dDZ4 = dZd - dZa;
	MINMAX4(dZa, dZb, dZc, dZd, dMinZ, dMaxZ);
	double dMax = MAX(MAX(dMaxX - dMinX, dMaxY - dMinY), dMaxZ - dMinZ);

	// dDet is proportional to the cell volume
	double dDet = (dDX2 * (dDY3 * dDZ4 - dDY4 * dDZ3)
			+ dDX3 * (dDY4 * dDZ2 - dDY2 * dDZ4)
			+ dDX4 * (dDY2 * dDZ3 - dDY3 * dDZ2));

	//   // Compute a length scale based on edge lengths.
	//   double dScale = ( dDIST3D(adA, adB) + dDIST3D(adA, adC) +
	// 		    dDIST3D(adA, adD) + dDIST3D(adB, adC) +
	// 		    dDIST3D(adB, adD) + dDIST3D(adC, adD) ) / 6.;

	//   dDet /= (dScale*dScale*dScale);

	//   double dError = 1.e-13;

	// Compute an upper bound on the error bound.

	const double dMachEps = 2.22044605e-13; // about 2^-52 * 1000;

	double dError = dMachEps * dMax * dMax * dMax;

	if (dDet > dError)
		return (1);
	else if (dDet < -dError)
		return (-1);

	//     If neither of those two worked, compute a more accurate error bound.
	//     The stuff in parentheses is the result of perturbing each term in
	//     the determinant by tweaking each multiplicand by the same amount.
	//     That amount is a multiple of machine zero.
	//       dError = dMachEps *
	//         (fabs(dDX2*dDY3) + fabs(dDX2*dDZ4) + fabs(dDY3*dDZ4) +
	//          fabs(dDX3*dDY4) + fabs(dDX3*dDZ2) + fabs(dDY4*dDZ2) +
	//          fabs(dDX4*dDY2) + fabs(dDX4*dDZ3) + fabs(dDY2*dDZ3) -
	//          fabs(dDZ2*dDY3) + fabs(dDZ2*dDX4) + fabs(dDY3*dDX4) +
	//          fabs(dDZ3*dDY4) + fabs(dDZ3*dDX2) + fabs(dDY4*dDX2) +
	//          fabs(dDZ4*dDY2) + fabs(dDZ4*dDX3) + fabs(dDY2*dDX3));
	//       if (dDet > dError)
	//         return (1);
	//       else if (dDet < -dError)
	//         return (-1);
	//       else
	return (0);
}

double tetVolume(const double coords0[], const double coords1[],
		const double coords2[], const double coords3[]) {
	double edge01[] = DIFF(coords1, coords0);
	double edge02[] = DIFF(coords2, coords0);
	double edge03[] = DIFF(coords3, coords0);
	double normal[3];
	CROSS(edge01, edge02, normal);
	double retVal = DOT(normal,edge03) / 6;
//	assert(retVal > 0);
	return retVal;
}

double pyrVolume(const double coords0[], const double coords1[],
		const double coords2[], const double coords3[], double coords4[]) {
	// Point 4 is the apex.
	double vecB[3], vecC[3], vecE[3];
	for (int ii = 0; ii < 3; ii++) {
		vecB[ii] = 0.25
				* (coords0[ii] + coords3[ii] - coords1[ii] - coords2[ii]);
		vecC[ii] = 0.25
				* (coords0[ii] + coords1[ii] - coords3[ii] - coords2[ii]);
		vecE[ii] =
				coords4[ii]
						- 0.25
								* (coords0[ii] + coords1[ii] + coords2[ii]
										+ coords3[ii]);
	}
	double normal[3];
	CROSS(vecB, vecC, normal);
	double retVal = DOT(normal, vecE) / 0.75;
//	assert(retVal > 0);
	return retVal;
}

