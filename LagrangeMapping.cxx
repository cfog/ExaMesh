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
 * LagrangeMapping.cxx
 *
 *  Created on: Oct. 1, 2019
 *      Author: cfog
 */

#include <assert.h>
#include "ExaMesh.h"

#include "Mapping.h"

LagrangeMapping::LagrangeMapping(const ExaMesh *const EM, const int nVals) :
		Mapping(EM), m_numValues(nVals) {
	m_nodalValues = new double[nVals][3];
}

LagrangeMapping::~LagrangeMapping() {
	delete[] m_nodalValues;
}

void LagrangeMapping::setupCoordMapping(const emInt verts[]) {
	double coords[m_numValues][3];
	for (int ii = 0; ii < m_numValues; ii++) {
		coords[ii][0] = m_pMesh->getX(verts[ii]);
		coords[ii][1] = m_pMesh->getY(verts[ii]);
		coords[ii][2] = m_pMesh->getZ(verts[ii]);
	}
	setNodalValues(coords);
	setModalValues();
}

void LagrangeMapping::setNodalValues(double inputValues[][3]) {
	for (int ii = 0; ii < m_numValues; ii++) {
		m_nodalValues[ii][0] = inputValues[ii][0];
		m_nodalValues[ii][1] = inputValues[ii][1];
		m_nodalValues[ii][2] = inputValues[ii][2];
	}
}

