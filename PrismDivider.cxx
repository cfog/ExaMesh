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
 * PrismDivider.cxx
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#include "PrismDivider.h"

void PrismDivider::setupCoordMapping(const emInt verts[]) {
	for (int ii = 0; ii < 6; ii++) {
		cellVerts[ii] = verts[ii];
	}
	m_Map->setupCoordMapping(verts);
}

void PrismDivider::getPhysCoordsFromParamCoords(const double uvw[3],
		double xyz[3]) {
	m_Map->computeTransformedCoords(uvw, xyz);
}

//
//void PrismDivider::setupCoordMapping(const emInt verts[]) {
//	for (int ii = 0; ii < 6; ii++) {
//		cellVerts[ii] = verts[ii];
//	}
//
//	double coords0[3], coords1[3], coords2[3], coords3[3], coords4[3], coords5[3];
//
//	m_pMesh->getCoords(verts[0], coords0);
//	m_pMesh->getCoords(verts[1], coords1);
//	m_pMesh->getCoords(verts[2], coords2);
//	m_pMesh->getCoords(verts[3], coords3);
//	m_pMesh->getCoords(verts[4], coords4);
//	m_pMesh->getCoords(verts[5], coords5);
//	for (int ii = 0; ii < 3; ii++) {
//		xyzOffsetBot[ii] = coords0[ii];
//		uVecBot[ii] = coords1[ii] - coords0[ii];
//		vVecBot[ii] = coords2[ii] - coords0[ii];
//		xyzOffsetTop[ii] = coords3[ii];
//		uVecTop[ii] = coords4[ii] - coords3[ii];
//		vVecTop[ii] = coords5[ii] - coords3[ii];
//	}
//}
//
//void PrismDivider::getPhysCoordsFromParamCoords(const double uvw[3],
//		double xyz[3]) {
//	const double& u = uvw[0];
//	const double& v = uvw[1];
//	const double& w = uvw[2];
//	double coordsBot[3], coordsTop[3];
//	for (int ii = 0; ii < 3; ii++) {
//		coordsBot[ii] = xyzOffsetBot[ii] + u * uVecBot[ii] + v * vVecBot[ii];
//		coordsTop[ii] = xyzOffsetTop[ii] + u * uVecTop[ii] + v * vVecTop[ii];
//		xyz[ii] = coordsBot[ii] * (1 - w) + coordsTop[ii] * w;
//	}
//}

void PrismDivider::divideInterior() {
  // Number of verts added:
  //    Tets:      (nD-1)(nD-2)(nD-3)/6
	if (nDivs < 3) return;
	double uvw[3];
	double &u = uvw[0];
	double& v = uvw[1];
	double& w = uvw[2];
	for (int kk = 1; kk < nDivs; kk++) {
		w = (1 - double(kk) / nDivs);
		for (int jj = 1; jj <= nDivs - 2; jj++) {
			v = double(jj) / nDivs;
			for (int ii = 1; ii <= nDivs - 1 - jj; ii++) {
				u = double(ii) / nDivs;
				double coordsNew[3];
				getPhysCoordsFromParamCoords(uvw, coordsNew);
				emInt vNew = m_pMesh->addVert(coordsNew);
				localVerts[ii][jj][kk] = vNew;
				m_uvw[ii][jj][nDivs - kk][0] = u;
				m_uvw[ii][jj][nDivs - kk][1] = v;
				m_uvw[ii][jj][nDivs - kk][2] = w;
      }
    } // Done looping over all interior verts for the triangle.
  }   // Done looping over all levels for the prism.
}

void PrismDivider::createNewCells() {
	// Output info about the points for this Prism, layer by layer
//	double newVol = 0;
//	for (int level = 0; level <= nDivs; level++) {
//		fprintf(stderr, "Level = %d\n", level);
//		for (int jj = 0; jj <= nDivs; jj++) {
//			for (int ii = 0; ii <= nDivs - jj; ii++) {
//				const double* coords = m_pMesh->getCoords(localVerts[ii][jj][level]);
//				fprintf(stderr, "  i = %2d, j = %2d:  (%10.6f, %10.6f, %10.6f)\n", ii,
//								jj, coords[0], coords[1], coords[2]);
//			}
//		}
//	}

	for (int level = 1; level <= nDivs; level++) {
//		fprintf(stderr, "Level: %d\n", level);
		// Create up-pointing Prisms.
		for (int jj = 0; jj <= nDivs - 1; jj++) {
      int ii = -1;
			for (ii = 0; ii <= nDivs - jj - 2; ii++) {
				emInt vertsNew1[] = { localVerts[ii][jj][level],
															localVerts[ii + 1][jj][level], localVerts[ii][jj
																	+ 1][level],
															localVerts[ii][jj][level - 1],
															localVerts[ii + 1][jj][level - 1],
															localVerts[ii][jj + 1][level - 1] };
				m_pMesh->addPrism(vertsNew1);
				assert(m_pMesh->numPrisms() < m_pMesh->maxNPrisms());

				// And now the other in that pair:
				emInt vertsNew2[] = { vertsNew1[1], localVerts[ii + 1][jj + 1][level],
															vertsNew1[2], vertsNew1[4], localVerts[ii + 1][jj
																	+ 1][level - 1],
															vertsNew1[5] };
				m_pMesh->addPrism(vertsNew2);
				assert(m_pMesh->numPrisms() < m_pMesh->maxNPrisms());
      } // Done with all prism pairs for this row.
      // Now one more at the end.
      ii = nDivs - jj - 1;
			emInt vertsNewLast[] = { localVerts[ii][jj][level],
																localVerts[ii + 1][jj][level], localVerts[ii][jj
																		+ 1][level],
																localVerts[ii][jj][level - 1],
																localVerts[ii + 1][jj][level - 1],
																localVerts[ii][jj + 1][level - 1] };
			m_pMesh->addPrism(vertsNewLast);
			assert(m_pMesh->numPrisms() <= m_pMesh->maxNPrisms());
    } // Done with this row (constant j)
  }   // Done with this level
//	logMessage(MSG_MANAGER, "  final volume: %G\n", newVol);
}
