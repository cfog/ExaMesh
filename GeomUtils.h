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
 * GeomUtils.h
 *
 *  Created on: Jul. 3, 2019
 *      Author: cfog
 */

#include <cmath>

double tetVolume(const double coords0[], const double coords1[],
		const double coords2[], const double coords3[]);

double pyrVolume(const double coords0[], const double coords1[],
		const double coords2[], const double coords3[], double coords4[]);

bool checkOrient3D(const double coordsA[3], const double coordsB[3],
		const double coordsC[3], const double coordsD[3]);

#define dDISTSQ3D(a, b)                                                          \
  ((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +         \
      (a[2] - b[2]) * (a[2] - b[2]))

