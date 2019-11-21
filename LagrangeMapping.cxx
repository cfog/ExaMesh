/*
 * Mapping.cxx
 *
 *  Created on: Oct. 1, 2019
 *      Author: cfog
 */

#include <assert.h>
#include <ExaMesh.h>

#include "Mapping.h"

LagrangeMapping::LagrangeMapping(const ExaMesh* const EM, const int nVals) :
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

