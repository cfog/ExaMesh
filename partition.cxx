/*
 * partition.cxx
 *
 *  Created on: Oct. 21, 2019
 *      Author: cfog
 */

#include <deque>
#include <vector>

#include "cgnslib.h"

#include "ExaMesh.h"
#include "Part.h"

void findCentroidOfVerts(const ExaMesh* const pEM, const emInt* verts,
		emInt nPts, double& x, double& y, double& z) {
	x = y = z = 0;
	for (emInt jj = 0; jj < nPts; jj++) {
		x += pEM->getX(verts[jj]);
		y += pEM->getY(verts[jj]);
		z += pEM->getZ(verts[jj]);
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

void addCellToPartitionData(const ExaMesh* const pEM, const emInt* verts,
		emInt nPts, emInt ii, int type, std::vector<CellPartData>& vecCPD,
		double& xmin, double& ymin, double& zmin, double& xmax, double& ymax,
		double& zmax) {
	{
		double x(0), y(0), z(0);
		findCentroidOfVerts(pEM, verts, nPts, x, y, z);
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
	for (emInt ii = 0; ii < pEM->numTets(); ii++) {
		const emInt* verts = pEM->getTetConn(ii);
		addCellToPartitionData(pEM, verts, 4, ii, TETRA_4, vecCPD, xmin, ymin, zmin,
														xmax, ymax, zmax);
	}
	for (emInt ii = 0; ii < pEM->numPyramids(); ii++) {
		const emInt* verts = pEM->getPyrConn(ii);
		addCellToPartitionData(pEM, verts, 5, ii, PYRA_5, vecCPD, xmin, ymin, zmin,
														xmax, ymax, zmax);
	}
	for (emInt ii = 0; ii < pEM->numPrisms(); ii++) {
		const emInt* verts = pEM->getPrismConn(ii);
		addCellToPartitionData(pEM, verts, 6, ii, PENTA_6, vecCPD, xmin, ymin, zmin,
														xmax, ymax, zmax);
	}
	for (emInt ii = 0; ii < pEM->numHexes(); ii++) {
		const emInt* verts = pEM->getHexConn(ii);
		addCellToPartitionData(pEM, verts, 8, ii, HEXA_8, vecCPD, xmin, ymin, zmin,
														xmax, ymax, zmax);
	}

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


