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
			emInt bdryTri = m_pMesh->addBdryTri(vertsNew1);
			assert(
					bdryTri < m_pMesh->numBdryTris() && bdryTri
							< m_pMesh->maxNBdryTris());

			// And now the other in that pair:
			emInt vertsNew2[] = { vertsNew1[1], localVerts[ii + 1][jj + 1][0],
														vertsNew1[2] };
			bdryTri = m_pMesh->addBdryTri(vertsNew2);
			assert(
					bdryTri < m_pMesh->numBdryTris() && bdryTri
							< m_pMesh->maxNBdryTris());
		} // Done with all prism pairs for this row.
		// Now one more at the end.
		ii = nDivs - jj - 1;
		emInt vertsNewLast[] = { localVerts[ii][jj][0], localVerts[ii + 1][jj][0],
															localVerts[ii][jj + 1][0] };
		emInt bdryTri = m_pMesh->addBdryTri(vertsNewLast);
		assert(
				bdryTri < m_pMesh->numBdryTris() && bdryTri < m_pMesh->maxNBdryTris());
	} // Done with this row (constant j)
}

