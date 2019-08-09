/*
 * PrismDivider.cxx
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#include "PrismDivider.h"

void PrismDivider::setupCoordMapping(const emInt verts[]) {
	const double *coords0 = m_pMesh->getCoords(verts[0]);
	const double *coords1 = m_pMesh->getCoords(verts[1]);
	const double *coords2 = m_pMesh->getCoords(verts[2]);
	const double *coords3 = m_pMesh->getCoords(verts[3]);
	const double *coords4 = m_pMesh->getCoords(verts[4]);
	const double *coords5 = m_pMesh->getCoords(verts[5]);
	for (int ii = 0; ii < 3; ii++) {
		xyzOffsetBot[ii] = coords0[ii];
		uVecBot[ii] = coords1[ii] - coords0[ii];
		vVecBot[ii] = coords2[ii] - coords0[ii];
		xyzOffsetTop[ii] = coords3[ii];
		uVecTop[ii] = coords4[ii] - coords3[ii];
		vVecTop[ii] = coords5[ii] - coords3[ii];
	}
}

void PrismDivider::getPhysCoordsFromParamCoords(const double uvw[3],
		double xyz[3]) {
	const double& u = uvw[0];
	const double& v = uvw[1];
	const double& w = uvw[2];
	double coordsBot[3], coordsTop[3];
	for (int ii = 0; ii < 3; ii++) {
		coordsBot[ii] = xyzOffsetBot[ii] + u * uVecBot[ii] + v * vVecBot[ii];
		coordsTop[ii] = xyzOffsetTop[ii] + u * uVecTop[ii] + v * vVecTop[ii];
		xyz[ii] = coordsBot[ii] * (1 - w) + coordsTop[ii] * w;
	}
}

void PrismDivider::divideInterior(const emInt verts[]) {
  // Number of verts added:
  //    Tets:      (nD-1)(nD-2)(nD-3)/6
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
				emInt prism = m_pMesh->addPrism(vertsNew1);
				assert(prism < m_pMesh->numPrisms() && prism < m_pMesh->maxNPrisms());

				// And now the other in that pair:
				emInt vertsNew2[] = { vertsNew1[1], localVerts[ii + 1][jj + 1][level],
															vertsNew1[2], vertsNew1[4], localVerts[ii + 1][jj
																	+ 1][level - 1],
															vertsNew1[5] };
				prism = m_pMesh->addPrism(vertsNew2);
				assert(prism < m_pMesh->numPrisms() && prism < m_pMesh->maxNPrisms());
      } // Done with all prism pairs for this row.
      // Now one more at the end.
      ii = nDivs - jj - 1;
			emInt vertsNewLast[] = { localVerts[ii][jj][level],
																localVerts[ii + 1][jj][level], localVerts[ii][jj
																		+ 1][level],
																localVerts[ii][jj][level - 1],
																localVerts[ii + 1][jj][level - 1],
																localVerts[ii][jj + 1][level - 1] };
			emInt prism = m_pMesh->addPrism(vertsNewLast);
			assert(prism < m_pMesh->numPrisms() && prism < m_pMesh->maxNPrisms());
    } // Done with this row (constant j)
  }   // Done with this level
//	logMessage(MSG_MANAGER, "  final volume: %G\n", newVol);
}
