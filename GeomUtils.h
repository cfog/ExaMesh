/*
 * GeomUtils.h
 *
 *  Created on: Jul. 3, 2019
 *      Author: cfog
 */

#include <cmath>

bool checkOrient3D(const double coordsA[3], const double coordsB[3],
		const double coordsC[3], const double coordsD[3]);

#define dDISTSQ3D(a, b)                                                          \
  ((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +         \
      (a[2] - b[2]) * (a[2] - b[2]))

