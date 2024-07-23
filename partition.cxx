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
 * partition.cxx
 *
 *  Created on: Oct. 21, 2019
 *      Author: cfog
 */

#include <deque>
#include <vector>

#include "ExaMesh.h"
#include "Part.h"

void ExaMesh::findCentroidOfVerts(const emInt* verts,
		emInt nPts, double& x,
		double& y, double& z) const {
	x = y = z = 0;
	for (emInt jj = 0; jj < nPts; jj++) {
		x += getX(verts[jj]);
		y += getY(verts[jj]);
		z += getZ(verts[jj]);
	}
	x /= nPts;
	y /= nPts;
	z /= nPts;
}

void extentBoundingBox(double x, double y, double z, double& xmin, double& ymin,
		double& zmin, double& xmax, double& ymax, double& zmax) {
	xmin = std::min(xmin, x);
	ymin = std::min(ymin, y);
	zmin = std::min(zmin, z);
	xmax = std::max(xmax, x);
	ymax = std::max(ymax, y);
	zmax = std::max(zmax, z);
}

void ExaMesh::addCellToPartitionData(const emInt* verts,
		emInt nPts, emInt ii, int type, std::vector<CellPartData>& vecCPD,
		double& xmin, double& ymin, double& zmin, double& xmax, double& ymax,
		double& zmax) const {
	{
		double x(0), y(0), z(0);
		findCentroidOfVerts(verts, nPts, x, y, z);
		extentBoundingBox(x, y, z, xmin, ymin, zmin, xmax, ymax, zmax);
		CellPartData CPD(ii, type, x, y, z);
		vecCPD.push_back(CPD);
	}
}

bool partitionCells(const ExaMesh* const pEM, const emInt nPartsToMake,
		std::vector<Part>& parts, std::vector<CellPartData>& vecCPD) {
	// Create collection of all cell (and bdry face) data, including info about
	// which entity it is.  Along the way, find the global bounding box.
	double xmin, xmax, ymin, ymax, zmin, zmax;
	xmin = ymin = zmin = DBL_MAX;
	xmax = ymax = zmax = -DBL_MAX;

	// Partitioning only cells, not bdry faces.  Also, currently no
	// cost differential for different cell types.
	pEM->setupCellDataForPartitioning(vecCPD, xmin, ymin, zmin, xmax, ymax, zmax);
	// Create a single part that contains all the cells, and put it in a deque.
	std::deque<Part> partsToSplit;

	Part P(0, vecCPD.size(), nPartsToMake, xmin, xmax, ymin, ymax, zmin, zmax);
	partsToSplit.push_back(P);

	// Grab a part from the deque to split.
	while (!partsToSplit.empty()) {
		P = partsToSplit.front();
		partsToSplit.pop_front();

		Part P1, P2;
		P.split(vecCPD, P1, P2);
		if (P1.numParts() > 1) {
			partsToSplit.push_back(P1);
		}
		else {
			parts.push_back(P1);
		}
		if (P2.numParts() > 1) {
			partsToSplit.push_back(P2);
		}
		else {
			parts.push_back(P2);
		}
	}
	return true;
}


