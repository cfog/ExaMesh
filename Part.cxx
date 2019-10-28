/*
 * Part.cxx
 *
 *  Created on: Oct. 21, 2019
 *      Author: cfog
 */

#include <algorithm>

#include <assert.h>

#include "exa-defs.h"
#include "Part.h"

class CellPartDataComparator {
	int whichVar;
public:
	CellPartDataComparator(const int which) :
			whichVar(which) {
	}
	bool operator()(const CellPartData& CPD1, const CellPartData& CPD2) const {
		return (CPD1.getCoord(whichVar) < CPD2.getCoord(whichVar));
	}
};

void Part::split(std::vector<CellPartData>& vCPD, Part& P1, Part& P2) const {
	assert(m_nParts > 1);
	// Find lengths in each coord direction; we're going to split the longest.
	double extents[] = { m_xmax - m_xmin, m_ymax - m_ymin, m_zmax - m_zmin };
	double mins[] = { m_xmin, m_ymin, m_zmin };
	double maxes[] = { m_xmax, m_ymax, m_zmax };
	int whichVar = -1;

	if (extents[0] > extents[1] && extents[0] > extents[2]) {
		// Sort in order of increasing x
		whichVar = 0;
	}
	else if (extents[1] > extents[2]) {
		// Sort in order of increasing y
		whichVar = 1;
	}
	else {
		// Sort in order of increasing z
		whichVar = 2;
	}

	CellPartDataComparator CPDC(whichVar);
	std::sort(vCPD.begin() + m_first, vCPD.begin() + m_last, CPDC);

	// Identify split point.  If there are going to be N parts made from this
	// one, then check at every 1/N of the number of cells, seeking the value
	// that is closest to bisecting cells in the direction we just sorted.
	// We could identify the point by binary search, but this should be so
	// much faster than the sort that it can't possibly matter at all.
	emInt numCells = m_last - m_first + 1;
	emInt divider = numCells / m_nParts + m_first;
	double divCoord = vCPD[divider].getCoord(whichVar);
	double bestFraction = (divCoord - mins[whichVar]) / extents[whichVar];
	emInt bestNParts = 1;
	for (emInt ii = 2; ii < m_nParts; ii++) {
		emInt candDivider = m_first + size_t(numCells) * ii / m_nParts;
		double candDivCoord = vCPD[candDivider].getCoord(whichVar);
		double thisFrac = (candDivCoord - mins[whichVar]) / extents[whichVar];
		if (fabs(thisFrac - 0.5) < fabs(bestFraction - 0.5)) {
			bestFraction = thisFrac;
			bestNParts = ii;
			divider = candDivider;
			divCoord = candDivCoord;
		}
		else {
			// Once we get past halfway, it'll never get any better again.
			break;
		}
	}

	// Now set up the new parts.
	double newMin1[] = { mins[0], mins[1], mins[2] };
	double newMax1[] = { maxes[0], maxes[1], maxes[2] };

	double newMin2[] = { mins[0], mins[1], mins[2] };
	double newMax2[] = { maxes[0], maxes[1], maxes[2] };

	newMax1[whichVar] = newMin2[whichVar] = divCoord;

	P1.setData(m_first, divider, bestNParts, newMin1, newMax1);
	P2.setData(divider, m_last, m_nParts - bestNParts, newMin2, newMax2);
}


