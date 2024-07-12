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
 * Part.h
 *
 *  Created on: Oct. 21, 2019
 *      Author: cfog
 */

#ifndef SRC_PART_H_
#define SRC_PART_H_

#include <vector>

#include <values.h>
#include "mpi.h"
#include <boost/serialization/access.hpp>
#include <boost/mpi/datatype.hpp>
class CellPartData {
private: 
	friend class boost::serialization::access; 
	template <class Archive> 
	void serialize(Archive &ar, const unsigned int /*version*/)	{
		ar &m_index; 
		ar &m_cellType; 
		ar &m_coords; 
	}
	emInt m_index, m_cellType;
	double m_coords[3];
public:
	// for registering MPI-type-CAUTION 
	CellPartData(){}; 
	CellPartData(const emInt ind, const emInt type, const double x,
			const double y, const double z) :
			m_index(ind), m_cellType(type) {
		m_coords[0] = x;
		m_coords[1] = y;
		m_coords[2] = z;
	}
	CellPartData (const emInt ind, const emInt type): m_index(ind),m_cellType(type){}; 
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
	friend MPI_Datatype register_mpi_type(CellPartData const&);
	friend bool operator==(const CellPartData& a, const CellPartData& b); 

};

class Part {
private: 
	friend class boost::serialization::access; 
	template <class Archive> 
	void serialize(Archive &ar, const unsigned int /*version*/)	{
		ar &m_xmin; 
		ar &m_xmax; 
		ar &m_ymin; 
		ar &m_ymax; 
		ar &m_zmin; 
		ar &m_zmax; 
		ar &m_first; 
		ar &m_last; 
		ar &m_nParts; 
	}

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
	
	friend MPI_Datatype register_mpi_type(Part const&);
	friend bool operator==(const Part& a, const Part& b); 

};
namespace boost { namespace mpi {
	template <>
struct is_mpi_datatype<Part> : mpl::true_ { };
} }
namespace boost { namespace mpi {
	template <>
struct is_mpi_datatype<CellPartData> : mpl::true_ { };
} }
#endif /* SRC_PART_H_ */
