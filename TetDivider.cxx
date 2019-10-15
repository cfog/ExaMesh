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

void TetDivider::divideInterior() {
  // Number of verts added:
  //    Tets:      (nD-1)(nD-2)(nD-3)/6
	double uvw[3];
	double& u = uvw[0];
	double& v = uvw[1];
	double& w = uvw[2];
	for (int kk = 0; kk <= nDivs - 4; kk++) {
		w = double(kk + 1) / nDivs;
		for (int jj = 0; jj <= nDivs - 4 - kk; jj++) {
			v = double(jj + 1) / nDivs;
			for (int ii = 0; ii <= nDivs - 4 - kk - jj; ii++) {
				u = double(ii + 1) / nDivs;
				assert(ii + jj + kk <= nDivs - 4);
				double coords[3];
				m_Map->computeTransformedCoords(uvw, coords);
				emInt vNew = m_pMesh->addVert(coords);
				localVerts[ii + 1][jj + 1][nDivs - (kk + 1)] = vNew;
      }
    }
  } // Done looping to create all verts inside the tet.
}

void TetDivider::stuffTetsIntoOctahedron(emInt vertsNew[][4]) {
	for (int ii = 0; ii < 4; ii++) {
		emInt tet = m_pMesh->addTet(vertsNew[ii]);
		assert(tet < m_pMesh->maxNTets());
#ifndef NDEBUG
		// TODO Re-enable this after you have a better criterion for
		// octa splitting.
//		assert(checkOrient3D(vertsNew[ii]) == 1);
#endif
	}
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
				assert(checkOrient3D(verts) == 1);
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
				emInt vertA = localVerts[ii][jj][level];
				emInt vertB = localVerts[ii][jj + 1][level];
				emInt vertC = localVerts[ii - 1][jj + 1][level];
				emInt vertD = localVerts[ii - 1][jj][level - 1];
				emInt vertE = localVerts[ii][jj][level - 1];
				emInt vertF = localVerts[ii - 1][jj + 1][level - 1];

				double coordsA[3], coordsB[3], coordsC[3], coordsD[3], coordsE[3],
						coordsF[3];
				m_pMesh->getCoords(vertA, coordsA);
				m_pMesh->getCoords(vertB, coordsB);
				m_pMesh->getCoords(vertC, coordsC);
				m_pMesh->getCoords(vertD, coordsD);
				m_pMesh->getCoords(vertE, coordsE);
				m_pMesh->getCoords(vertF, coordsF);

				useDiagAF = useDiagBD = useDiagCE = false;
				double distAF = dDIST3D(coordsA, coordsF);
				double distBD = dDIST3D(coordsB, coordsD);
				double distCE = dDIST3D(coordsC, coordsE);
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
					stuffTetsIntoOctahedron(vertsNew);
        } else if (useDiagBD) {
          assert(!useDiagAF && !useDiagCE);
					emInt vertsNew[][4] = { { vertC, vertA, vertB, vertD },
																	{ vertA, vertE, vertB, vertD },
																	{ vertE, vertF, vertB, vertD },
																	{ vertF, vertC, vertB, vertD } };
					stuffTetsIntoOctahedron(vertsNew);
        } else {
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
