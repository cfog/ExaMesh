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
 * PyrDivider.cxx
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#include "GeomUtils.h"
#include "PyrDivider.h"

void PyrDivider::setupCoordMapping(const emInt verts[]) {
	for (int ii = 0; ii < 5; ii++) {
		cellVerts[ii] = verts[ii];
	}
	m_Map->setupCoordMapping(verts);
}

void PyrDivider::getPhysCoordsFromParamCoords(const double uvw[3],
		double xyz[3]) {
	m_Map->computeTransformedCoords(uvw, xyz);
}

//void PyrDivider::setupCoordMapping(const emInt verts[]) {
//	for (int ii = 0; ii < 5; ii++) {
//		cellVerts[ii] = verts[ii];
//	}
//
//	double coords0[3], coords1[3], coords2[3], coords3[3], coords4[3];
//
//	m_pMesh->getCoords(verts[0], coords0);
//	m_pMesh->getCoords(verts[1], coords1);
//	m_pMesh->getCoords(verts[2], coords2);
//	m_pMesh->getCoords(verts[3], coords3);
//	m_pMesh->getCoords(verts[4], coords4);
//
//	for (int ii = 0; ii < 3; ii++) {
//		xyzOffset[ii] = coords0[ii];
//		uVec[ii] = coords1[ii] - coords0[ii];
//		vVec[ii] = coords3[ii] - coords0[ii];
//		uvVec[ii] = coords2[ii] + coords0[ii] - coords1[ii] - coords3[ii];
//		xyzApex[ii] = coords4[ii];
//	}
//}

//void PyrDivider::getPhysCoordsFromParamCoords(const double uvw[3],
//		double xyz[3]) {
//	double u = uvw[0];
//	double v = uvw[1];
//	double w = uvw[2];
//	if (u != 0 && u != 1 && w != 0) u /= (1 - w);
//	if (v != 0 && v != 1 && w != 0) v /= (1 - w);
//	double coordsBase[3];
//	for (int ii = 0; ii < 3; ii++) {
//		coordsBase[ii] = xyzOffset[ii] + u * uVec[ii] + v * vVec[ii]
//											+ u * v * uvVec[ii];
//		xyz[ii] = coordsBase[ii] * (1 - w) + xyzApex[ii] * w;
//	}
//}


//void PyrDivider::divideInterior() {
//  // Number of verts added:
//	//    Pyrs:      (nD-1)(nD-2)(2 nD-3)/6
//	if (nDivs < 3) return;
//
//	for (int kk = 1; kk < nDivs; kk++) {
//		int jMax = maxJ(1, kk);
//		for (int jj = 1; jj < jMax; jj++) {
//			int iMax = maxI(1, kk);
//			for (int ii = 1; ii < iMax; ii++) {
//				double uvw[] = { double(ii) / nDivs,
//						double(jj) / nDivs,
//						1 - double(kk) / nDivs };
//				double &u = uvw[0];
//				double &v = uvw[1];
//				double &w = uvw[2];
//				double coords[3];
//				m_Map->computeTransformedCoords(uvw, coords);
//				emInt vNew = m_pMesh->addVert(coords);
//				localVerts[ii][jj][kk] = vNew;
//				m_uvw[ii][jj][kk][0] = u;
//				m_uvw[ii][jj][kk][1] = v;
//				m_uvw[ii][jj][kk][2] = w;
////				printf("%3d %3d %3d %5f %5f %5f\n", ii, jj,
////						kk, u, v, w);
//			}
//		}
//	} // Done looping to create all verts inside the tet.
////	printf("\n");
////	for (int kk = 2; kk <= nDivs - 1; kk++) {
////		for (int jj = 1; jj <= kk - 1; jj++) {
////			for (int ii = 1; ii <= kk - 1; ii++) {
////				double uvw[] = { double(ii) / nDivs,
////						double(jj) / nDivs,
////						1 - double(kk) / nDivs };
////				double &u = uvw[0];
////				double &v = uvw[1];
////				double &w = uvw[2];
////				double coords[3];
////				m_Map->computeTransformedCoords(uvw, coords);
////				emInt vNew = m_pMesh->addVert(coords);
////				localVerts[ii][jj][kk] = vNew;
////				m_uvw[ii][jj][kk][0] = u;
////				m_uvw[ii][jj][kk][1] = v;
////				m_uvw[ii][jj][kk][2] = w;
////				printf("%3d %3d %3d %5f %5f %5f\n", ii, jj,
////						kk, u, v, w);
////			}
////		}
////	} // Done looping to create all verts inside the tet.
//}

void PyrDivider::createNewCells() {
  // Output info about the points for this pyramid, layer by layer

//	for (int level = 0; level <= nDivs; level++) {
//		fprintf(stderr, "Level = %d\n", level);
//    for (int jj = 0; jj <= level; jj++) {
//      for (int ii = 0; ii <= level; ii++) {
//				emInt vert = localVerts[ii][jj][level];
//				fprintf(stderr, "  i = %2d, j = %2d:  (%10.6f, %10.6f, %10.6f)\n", ii,
//								jj, m_pMesh->getCoords(vert)[0], m_pMesh->getCoords(vert)[1],
//								m_pMesh->getCoords(vert)[2]);
//      }
//    }
//  }

	for (int level = 1; level <= nDivs; level++) {
//		fprintf(stderr, "Level: %d\n", level);
    // Create up-pointing pyrs.  For a given level, there are
    // level^2 of these.
    for (int jj = 0; jj < level; jj++) {
      for (int ii = 0; ii < level; ii++) {
				emInt vertsNew[] = { localVerts[ii][jj][level],
															localVerts[ii + 1][jj][level],
															localVerts[ii + 1][jj + 1][level],
															localVerts[ii][jj + 1][level],
															localVerts[ii][jj][level - 1] };
				m_pMesh->addPyramid(vertsNew);
				assert(m_pMesh->numPyramids() < m_pMesh->maxNPyrs());
      }
    }

    // Down-pointing pyramids, hanging from quads on the previous level.  There
    // are (level-1)^2 of these.
    for (int jj = 0; jj <= level - 2; jj++) {
      for (int ii = 0; ii <= level - 2; ii++) {
				emInt vertsNew[] = { localVerts[ii][jj][level - 1], localVerts[ii][jj
						+ 1][level - 1],
															localVerts[ii + 1][jj + 1][level - 1],
															localVerts[ii + 1][jj][level - 1], localVerts[ii
																	+ 1][jj + 1][level] };
				m_pMesh->addPyramid(vertsNew);
				assert(m_pMesh->numPyramids() <= m_pMesh->maxNPyrs());
      }
    }

    // Now there are tets in the gaps between these pyramids.  There are
    // 2 (level-1) (level-2) of these.

    // The set on lines of constant j on level l.
    for (int jj = 1; jj <= level - 1; jj++) {
      for (int ii = 0; ii <= level - 1; ii++) {
				emInt vertsNew[] = { localVerts[ii][jj][level],
															localVerts[ii + 1][jj][level],
															localVerts[ii][jj][level - 1], localVerts[ii][jj
																	- 1][level - 1] };
				m_pMesh->addTet(vertsNew);
				assert(m_pMesh->numTets() < m_pMesh->maxNTets());
				assert(checkOrient3D(vertsNew) == 1);
      }
    }

    // The set on lines of constant i on level l.
    for (int jj = 0; jj <= level - 1; jj++) {
      for (int ii = 1; ii <= level - 1; ii++) {
				emInt vertsNew[] = { localVerts[ii][jj][level],
															localVerts[ii][jj + 1][level],
															localVerts[ii - 1][jj][level - 1],
															localVerts[ii][jj][level - 1] };
				m_pMesh->addTet(vertsNew);
				assert(m_pMesh->numTets() <= m_pMesh->maxNTets());
				assert(checkOrient3D(vertsNew) == 1);
      }
    }

  } // Done with this level
}
