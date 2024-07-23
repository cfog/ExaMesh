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
 * TetDivider.cxx
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#include "GeomUtils.h"
#include "TetDivider.h"

void TetDivider::setupCoordMapping(const emInt verts[]) {
	for (int ii = 0; ii < 4; ii++) {
		cellVerts[ii] = verts[ii];
	}
	m_Map->setupCoordMapping(verts);
}

void TetDivider::getPhysCoordsFromParamCoords(const double uvw[3],
		double xyz[3]) {
	m_Map->computeTransformedCoords(uvw, xyz);
}

//void TetDivider::divideInterior() {
//	// Number of verts added:
//	//    Tets:      (nD-1)(nD-2)(nD-3)/6
//	if (nDivs <= 3) return;
//	// No need to include bdry points in the iteration
//	for (int kk = 1; kk < nDivs; kk++) {
//		int jMax = maxJ(1, kk);
//		for (int jj = 1; jj < jMax; jj++) {
//			int iMax = maxI(jj, kk);
//			for (int ii = 1; ii < iMax; ii++) {
//				assert(ii + jj < kk);
//
//				double uvw[] = {double(ii) / nDivs,
//						double(jj) / nDivs,
//						1 - double(kk) / nDivs
//				};
//				double coords[3];
//
//				m_Map->computeTransformedCoords(uvw, coords);
//				emInt vNew = m_pMesh->addVert(coords);
//
//				localVerts[ii][jj][kk] = vNew;
//				m_uvw[ii][jj][kk][0] = uvw[0];
//				m_uvw[ii][jj][kk][1] = uvw[1];
//				m_uvw[ii][jj][kk][2] = uvw[2];
////				printf("%3d %3d %3d %8f %8f %8f\n", ii, jj, kk,
////						u, v, w);
//			}
//		}
//	}
//}

void TetDivider::stuffTetsIntoOctahedron(emInt vertsNew[][4]) {
	m_pMesh->addTet(vertsNew[0]);
	m_pMesh->addTet(vertsNew[1]);
	m_pMesh->addTet(vertsNew[2]);
	m_pMesh->addTet(vertsNew[3]);
	assert(m_pMesh->numTets() <= m_pMesh->maxNTets());
#ifndef NDEBUG
	assert(checkOrient3D(vertsNew[0]) != -1);
	assert(checkOrient3D(vertsNew[1]) != -1);
	assert(checkOrient3D(vertsNew[2]) != -1);
	assert(checkOrient3D(vertsNew[3]) != -1);
#endif
}

void TetDivider::createNewCells() {
//		// Output info about the points for this tet, layer by layer
//
//	for (int level = 0; level <= nDivs; level++) {
//		fprintf(stderr, "Level = %d\n", level);
//		for (int jj = 0; jj <= level; jj++) {
//			for (int ii = 0; ii <= level - jj; ii++) {
//				fprintf(stderr, "  i = %2d, j = %2d, v = %4u:  ", ii, jj,
//								localVerts[ii][jj][level]);
//				const double* coords = m_pMesh->getCoords(localVerts[ii][jj][level]);
//				fprintf(stderr, "(%10.6f, %10.6f, %10.6f)\n", coords[0], coords[1],
//								coords[2]);
//			}
//		}
//	}

	// Level 1 is the tip of the original tet; level nDivs
	// is at the bottom (level with most entities).
	for (int level = 1; level <= int(nDivs); level++) {
//		fprintf(stderr, "Level: %d\n", level);
		// Create up-pointing tets.  For a given level, there are
		// (level+1)(level)/2 of these.
		for (int jj = 0; jj < level; jj++) {
			for (int ii = 0; ii < level - jj; ii++) {
				int kk = nDivs - level;
				//					logMessage(
				//							MSG_DEBUG,
				//							"Tet: (%d, %d, %d), (%d, %d,
				//%d), (%d, %d, %d), (%d, %d, %d)\n", 							ii, jj, level, ii + 1, jj, level,
				//ii, jj + 1, level, ii, jj, 							level - 1);
				emInt vert0 = localVerts[ii][jj][kk];
				emInt vert1 = localVerts[ii + 1][jj][kk];
				emInt vert2 = localVerts[ii][jj + 1][kk];
				emInt vert3 = localVerts[ii][jj][kk + 1];
				emInt verts[] = { vert0, vert1, vert2, vert3 };

				m_pMesh->addTet(verts);
				assert(m_pMesh->numTets() <= m_pMesh->maxNTets());
				assert(checkOrient3D(verts) == 1);
			}
		}

		// Each down-point triangle on the previous level has a tet
		// that extends down to a point on this level. There are
		// (levels-1)*(levels-2)/2 of these.
		for (int jj = 0; jj <= level - 3; jj++) {
			for (int ii = 1; ii <= level - jj - 2; ii++) {
				int kk = nDivs - level;
				emInt vert0 = localVerts[ii][jj][kk + 1];
				emInt vert1 = localVerts[ii - 1][jj + 1][kk + 1];
				emInt vert2 = localVerts[ii][jj + 1][kk + 1];
				emInt vert3 = localVerts[ii][jj + 1][kk];
				emInt verts[] = { vert0, vert1, vert2, vert3 };
				//					logMessage(
				//							MSG_DEBUG,
				//							"Tet: (%d, %d, %d), (%d, %d,
				//%d), (%d, %d, %d), (%d, %d, %d)\n", 							ii, jj, level - 1, ii - 1, jj + 1,
				//level - 1, ii, jj + 1, 							level - 1, ii, jj + 1, level);
				m_pMesh->addTet(verts);
				assert(m_pMesh->numTets() <= m_pMesh->maxNTets());
				assert(checkOrient3D(verts) == 1);
			}
		}

		// The rest of the tris (down-pointing) in this level have an octahedron
		// on them, connecting them to an (up-pointing) tri the level above.
		// There are (level)(level-1)/2 octahedra, each of which will be split
		// into four tetrahedra.

		bool useDiagAF = true, useDiagBD = false, useDiagCE = false;
		for (int jj = 0; jj <= level - 2; jj++) {
			for (int ii = 1; ii <= level - jj - 1; ii++) {
				int kk = nDivs - level;
				emInt vertA = localVerts[ii][jj][kk];
				emInt vertB = localVerts[ii][jj + 1][kk];
				emInt vertC = localVerts[ii - 1][jj + 1][kk];
				emInt vertD = localVerts[ii - 1][jj][kk + 1];
				emInt vertE = localVerts[ii][jj][kk + 1];
				emInt vertF = localVerts[ii - 1][jj + 1][kk + 1];

				double coordsA[3], coordsB[3], coordsC[3], coordsD[3], coordsE[3],
						coordsF[3];
				m_pMesh->getCoords(vertA, coordsA);
				m_pMesh->getCoords(vertB, coordsB);
				m_pMesh->getCoords(vertC, coordsC);
				m_pMesh->getCoords(vertD, coordsD);
				m_pMesh->getCoords(vertE, coordsE);
				m_pMesh->getCoords(vertF, coordsF);

				useDiagAF = useDiagBD = useDiagCE = false;
				double distsqAF = dDISTSQ3D(coordsA, coordsF);
				double distsqBD = dDISTSQ3D(coordsB, coordsD);
				double distsqCE = dDISTSQ3D(coordsC, coordsE);
				if (distsqAF <= distsqBD && distsqAF <= distsqCE) {
					useDiagAF = true;
				}
				else if (distsqBD <= distsqCE) {
					useDiagBD = true;
				}
				else {
					useDiagCE = true;
				}

				if (useDiagAF) {
					assert(!useDiagBD && !useDiagCE);
					emInt vertsNew[][4] = { { vertB, vertC, vertA, vertF },
							{ vertC, vertD, vertA, vertF },
							{ vertD, vertE, vertA, vertF },
							{ vertE, vertB, vertA, vertF } };
					stuffTetsIntoOctahedron(vertsNew);
				}
				else if (useDiagBD) {
					assert(!useDiagAF && !useDiagCE);
					emInt vertsNew[][4] = { { vertC, vertA, vertB, vertD },
							{ vertA, vertE, vertB, vertD },
							{ vertE, vertF, vertB, vertD },
							{ vertF, vertC, vertB, vertD } };
					stuffTetsIntoOctahedron(vertsNew);
				}
				else {
					assert(!useDiagAF && !useDiagBD);
					emInt vertsNew[][4] = { { vertA, vertB, vertC, vertE },
							{ vertB, vertF, vertC, vertE },
							{ vertF, vertD, vertC, vertE },
							{ vertD, vertA, vertC, vertE } };
					stuffTetsIntoOctahedron(vertsNew);
				}
			}
		} // Done with octahedra
	}   // Done with this level
}
