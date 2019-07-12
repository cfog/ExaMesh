/*
 * TetDivider.cxx
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#include "GeomUtils.h"
#include "TetDivider.h"

void TetDivider::divideInterior(const emInt verts[]) {
  // Number of verts added:
  //    Tets:      (nD-1)(nD-2)(nD-3)/6
	const double* coords0 = m_pMesh->getCoords(verts[0]);
	const double* coords1 = m_pMesh->getCoords(verts[1]);
	const double* coords2 = m_pMesh->getCoords(verts[2]);
	const double* coords3 = m_pMesh->getCoords(verts[3]);

	double deltaInI[] = { (coords1[0] - coords0[0]) / nDivs, (coords1[1]
			- coords0[1])
																														/ nDivs,
												(coords1[2] - coords0[2]) / nDivs };
	double deltaInJ[] = { (coords2[0] - coords0[0]) / nDivs, (coords2[1]
			- coords0[1])
																														/ nDivs,
												(coords2[2] - coords0[2]) / nDivs };
	double deltaInK[] = { (coords3[0] - coords0[0]) / nDivs, (coords3[1]
			- coords0[1])
																														/ nDivs,
												(coords3[2] - coords0[2]) / nDivs };
	for (int kk = 0; kk <= nDivs - 4; kk++) {
		for (int jj = 0; jj <= nDivs - 4 - kk; jj++) {
			for (int ii = 0; ii <= nDivs - 4 - kk - jj; ii++) {
				assert(ii + jj + kk <= nDivs - 4);
				double coords[] = { coords0[0] + deltaInI[0] * (ii + 1)
														+
                               deltaInJ[0] * (jj + 1) + deltaInK[0] * (kk + 1),
														coords0[1] + deltaInI[1] * (ii + 1)
														+
                               deltaInJ[1] * (jj + 1) + deltaInK[1] * (kk + 1),
														coords0[2] + deltaInI[2] * (ii + 1)
														+
                               deltaInJ[2] * (jj + 1) + deltaInK[2] * (kk + 1)};
				emInt vNew = m_pMesh->addVert(coords);
				localVerts[ii + 1][jj + 1][nDivs - (kk + 1)] = vNew;
      }
    }
  } // Done looping to create all verts inside the tet.
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

	for (int level = 1; level <= nDivs; level++) {
//		fprintf(stderr, "Level: %d\n", level);
    // Create up-pointing tets.  For a given level, there are
    // (level+1)(level)/2 of these.
    for (int jj = 0; jj < level; jj++) {
      for (int ii = 0; ii < level - jj; ii++) {
        //					logMessage(
        //							MSG_DEBUG,
        //							"Tet: (%d, %d, %d), (%d, %d,
        //%d), (%d, %d, %d), (%d, %d, %d)\n", 							ii, jj, level, ii + 1, jj, level,
        //ii, jj + 1, level, ii, jj, 							level - 1);
				emInt vert0 = localVerts[ii][jj][level];
				emInt vert1 = localVerts[ii + 1][jj][level];
				emInt vert2 = localVerts[ii][jj + 1][level];
				emInt vert3 = localVerts[ii][jj][level - 1];
				emInt verts[] = { vert0, vert1, vert2, vert3 };

				emInt tet = m_pMesh->addTet(verts);
				assert(tet < m_pMesh->maxNTets());
				assert(
						checkOrient3D(m_pMesh->getCoords(vert0), m_pMesh->getCoords(vert1),
													m_pMesh->getCoords(vert2),
													m_pMesh->getCoords(vert3)));
      }
    }

		// Each down-point triangle on the previous level has a tet
    // that extends down to a point on this level. There are
    // (levels-1)*(levels-2)/2 of thes.
    for (int jj = 0; jj <= level - 3; jj++) {
      for (int ii = 1; ii <= level - jj - 2; ii++) {
				emInt vert0 = localVerts[ii][jj][level - 1];
				emInt vert1 = localVerts[ii - 1][jj + 1][level - 1];
				emInt vert2 = localVerts[ii][jj + 1][level - 1];
				emInt vert3 = localVerts[ii][jj + 1][level];
				emInt verts[] = { vert0, vert1, vert2, vert3 };
        //					logMessage(
        //							MSG_DEBUG,
        //							"Tet: (%d, %d, %d), (%d, %d,
        //%d), (%d, %d, %d), (%d, %d, %d)\n", 							ii, jj, level - 1, ii - 1, jj + 1,
        //level - 1, ii, jj + 1, 							level - 1, ii, jj + 1, level);
				emInt tet = m_pMesh->addTet(verts);
				assert(tet < m_pMesh->maxNTets());
				assert(
						checkOrient3D(m_pMesh->getCoords(vert0), m_pMesh->getCoords(vert1),
													m_pMesh->getCoords(vert2),
													m_pMesh->getCoords(vert3)));
      }
    }

    // The rest of the tris (down-pointing) in this level have an octahedron
    // on them, connecting them to an (up-pointing) tri the level above.
    // There are (level)(level-1)/2 octahedra, each of which will be split
    // into four tetrahedra.

    // TODO: Decide which octahedron template to use based on mesh data
    bool useDiagAF = true, useDiagBD = false, useDiagCE = false;
		for (int jj = 0; jj <= level - 2; jj++) {
			for (int ii = 1; ii <= level - jj - 1; ii++) {
				emInt vertA = localVerts[ii][jj][level];
				emInt vertB = localVerts[ii][jj + 1][level];
				emInt vertC = localVerts[ii - 1][jj + 1][level];
				emInt vertD = localVerts[ii - 1][jj][level - 1];
				emInt vertE = localVerts[ii][jj][level - 1];
				emInt vertF = localVerts[ii - 1][jj + 1][level - 1];

				useDiagAF = useDiagBD = useDiagCE = false;
				double distAF = dDIST3D(m_pMesh->getCoords(vertA),
																m_pMesh->getCoords(vertF));
				double distBD = dDIST3D(m_pMesh->getCoords(vertB),
																m_pMesh->getCoords(vertD));
				double distCE = dDIST3D(m_pMesh->getCoords(vertC),
																m_pMesh->getCoords(vertE));
				if (distAF <= distBD && distAF <= distCE) {
					useDiagAF = true;
				}
				else if (distBD <= distCE) {
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
					for (int ii = 0; ii < 4; ii++) {

						emInt tet = m_pMesh->addTet(vertsNew[ii]);
						assert(tet < m_pMesh->maxNTets());
						assert(
								checkOrient3D(m_pMesh->getCoords(vertsNew[ii][0]),
															m_pMesh->getCoords(vertsNew[ii][1]),
															m_pMesh->getCoords(vertsNew[ii][2]),
															m_pMesh->getCoords(vertsNew[ii][3])));
					}
        } else if (useDiagBD) {
          assert(!useDiagAF && !useDiagCE);
					emInt vertsNew[][4] = { { vertC, vertA, vertB, vertD },
																	{ vertA, vertE, vertB, vertD },
																	{ vertE, vertF, vertB, vertD },
																	{ vertF, vertC, vertB, vertD } };
					for (int ii = 0; ii < 4; ii++) {

						emInt tet = m_pMesh->addTet(vertsNew[ii]);
						assert(tet < m_pMesh->maxNTets());
						assert(
								checkOrient3D(m_pMesh->getCoords(vertsNew[ii][0]),
															m_pMesh->getCoords(vertsNew[ii][1]),
															m_pMesh->getCoords(vertsNew[ii][2]),
															m_pMesh->getCoords(vertsNew[ii][3])));
					}
        } else {
          assert(!useDiagAF && !useDiagBD);
					emInt vertsNew[][4] = { { vertA, vertB, vertC, vertE },
																	{ vertB, vertF, vertC, vertE },
																	{ vertF, vertD, vertC, vertE },
																	{ vertD, vertA, vertC, vertE } };
					for (int ii = 0; ii < 4; ii++) {

						emInt tet = m_pMesh->addTet(vertsNew[ii]);
						assert(tet < m_pMesh->maxNTets());
						assert(
								checkOrient3D(m_pMesh->getCoords(vertsNew[ii][0]),
															m_pMesh->getCoords(vertsNew[ii][1]),
															m_pMesh->getCoords(vertsNew[ii][2]),
															m_pMesh->getCoords(vertsNew[ii][3])));
					}
        }
      }
    } // Done with octahedra
  }   // Done with this level
}