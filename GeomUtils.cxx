/*
 * GeomUtils.cxx
 *
 *  Created on: Jul. 3, 2019
 *      Author: cfog
 */

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
			+ dDX3 * (dDY4 * dDZ2 - dDY2 * dDZ4) + dDX4 * (dDY2 * dDZ3 - dDY3 * dDZ2));

	//   // Compute a length scale based on edge lengths.
	//   double dScale = ( dDIST3D(adA, adB) + dDIST3D(adA, adC) +
	// 		    dDIST3D(adA, adD) + dDIST3D(adB, adC) +
	// 		    dDIST3D(adB, adD) + dDIST3D(adC, adD) ) / 6.;

	//   dDet /= (dScale*dScale*dScale);

	//   double dError = 1.e-13;

	// Compute an upper bound on the error bound.

	const double dMachEps = 2.22044605e-13; // about 2^-52 * 1000;

	double dError = dMachEps * dMax * dMax * dMax;

	if (dDet > dError) return (1);
	else if (dDet < -dError) return (-1);

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




