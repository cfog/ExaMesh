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
 * BdryTriDivider.cxx
 *
 *  Created on: July 8, 2019
 *      Author: cfog
 */

#include "BdryTriDivider.h"

void BdryTriDivider::divideInterior() {
}

void BdryTriDivider::createNewCells() {
	// Okay, sure, these aren't actually cells in the usual sense, but so what?
	// Create topologically up-pointing triangles.
	for (int jj = 0; jj <= nDivs - 1; jj++) {
		int ii = -1;
		for (ii = 0; ii <= nDivs - jj - 2; ii++) {
			emInt vertsNew1[] = { localVerts[ii][jj][0], localVerts[ii + 1][jj][0],
														localVerts[ii][jj + 1][0] };
			m_pMesh->addBdryTri(vertsNew1);
			assert(m_pMesh->numBdryTris() <= m_pMesh->maxNBdryTris());

			// And now the other in that pair:
			emInt vertsNew2[] = { vertsNew1[1], localVerts[ii + 1][jj + 1][0],
														vertsNew1[2] };
			m_pMesh->addBdryTri(vertsNew2);
			assert(m_pMesh->numBdryTris() < m_pMesh->maxNBdryTris());
		} // Done with all prism pairs for this row.
		// Now one more at the end.
		ii = nDivs - jj - 1;
		emInt vertsNewLast[] = { localVerts[ii][jj][0], localVerts[ii + 1][jj][0],
															localVerts[ii][jj + 1][0] };

		m_pMesh->addBdryTri(vertsNewLast);
		assert(m_pMesh->numBdryTris() <= m_pMesh->maxNBdryTris());
	} // Done with this row (constant j)
}

void BdryTriDivider::setRefinedVerts(TriFaceVerts &TF){
	for (int ii = 0; ii <= nDivs ; ii++) {
	 	for (int jj = 0; jj <= nDivs-ii ; jj++) {
			int trueI;
			int trueJ;
			emInt vert= localVerts[ii][jj][0]; 
			TF.setIntVertInd(ii,jj,vert); 
		}
	}	
}