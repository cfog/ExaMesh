/*
 * BdryQuadDivider.cxx
 *
 *  Created on: July 8, 2019
 *      Author: cfog
 */

#include "BdryQuadDivider.h"

void BdryQuadDivider::divideInterior(const emInt[]) {
}

void BdryQuadDivider::createNewCells() {
	// Okay, sure, these aren't actually cells in the usual sense, but so what?
	for (int jj = 0; jj <= nDivs - 1; jj++) {
		int ii = -1;
		for (ii = 0; ii <= nDivs - 1; ii++) {
			emInt vertsNew1[] = { localVerts[ii][jj][0], localVerts[ii + 1][jj][0],
														localVerts[ii + 1][jj + 1][0],
														localVerts[ii][jj + 1][0] };
			emInt bdryQuad = m_pMesh->addBdryQuad(vertsNew1);
			assert(
					bdryQuad < m_pMesh->numBdryQuads() && bdryQuad
							< m_pMesh->maxNBdryQuads());
		} // Done with all quads for this row.
	} // Done with this row (constant j)
}

