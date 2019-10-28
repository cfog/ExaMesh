/*
 * Part.h
 *
 *  Created on: Oct. 21, 2019
 *      Author: cfog
 */

#ifndef SRC_PART_H_
#define SRC_PART_H_

#include <vector>

#include <values.h>

class CellPartData {
	emInt m_index, m_cellType;
	double m_coords[3];
public:
	CellPartData(const emInt ind, const emInt type, const double x,
			const double y, const double z) :
			m_index(ind), m_cellType(type) {
		m_coords[0] = x;
		m_coords[1] = y;
		m_coords[2] = z;
	}
	double getCoord(const int which) const {
		assert(which >= 0 && which < 3);
		return m_coords[which];
	}

	emInt getCellType() const {
		return m_cellType;
	}

	emInt getIndex() const {
		return m_index;
	}
};

class Part {
	double m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax;
	emInt m_first, m_last, m_nParts;
public:
	Part() :
			m_xmin(DBL_MAX), m_xmax(-DBL_MAX), m_ymin(DBL_MAX), m_ymax(-DBL_MAX),
					m_zmin(DBL_MAX), m_zmax(-DBL_MAX), m_first(0), m_last(0), m_nParts(0) {
	}
	Part(const emInt first, const emInt last, const emInt partsToMake,
			const double xmin, const double xmax, const double ymin,
			const double ymax, const double zmin, const double zmax) :
			m_xmin(xmin), m_xmax(xmax), m_ymin(ymin), m_ymax(ymax), m_zmin(zmin),
					m_zmax(zmax), m_first(first), m_last(last), m_nParts(partsToMake) {
	}
	void setData(const emInt first, const emInt last, const emInt partsToMake,
			const double mins[3], const double maxes[3]) {
		m_xmin = mins[0];
		m_xmax = maxes[0];
		m_ymin = mins[1];
		m_ymax = maxes[1];
		m_zmin = mins[2];
		m_zmax = maxes[2];
		m_first = first;
		m_last = last;
		m_nParts = partsToMake;
	}
	emInt numParts() const {
		return m_nParts;
	}
	void split(std::vector<CellPartData>& vCPD, Part& P1, Part& P2) const;

	emInt getFirst() const {
		return m_first;
	}

	emInt getLast() const {
		return m_last;
	}

	double getXmax() const {
		return m_xmax;
	}

	double getXmin() const {
		return m_xmin;
	}

	double getYmax() const {
		return m_ymax;
	}

	double getYmin() const {
		return m_ymin;
	}

	double getZmax() const {
		return m_zmax;
	}

	double getZmin() const {
		return m_zmin;
	}
};


#endif /* SRC_PART_H_ */
