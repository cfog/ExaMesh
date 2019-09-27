/*
 * TetDivider.cxx
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#include "GeomUtils.h"
#include "TetDivider.h"

void TetDivider::setPolyCoeffs(const double* xyz0, const double* xyz1,
		const double* xyz2, const double* xyz3, double uderiv0[3],
		double vderiv0[3], double wderiv0[3], double uderiv1[3], double vderiv1[3],
		double wderiv1[3], double uderiv2[3], double vderiv2[3], double wderiv2[3],
		double uderiv3[3], double vderiv3[3], double wderiv3[3]) {
	for (int ii = 0; ii < 3; ii++) {
		//		xyzOffset[ii] = xyz0[ii];
		//		uVec[ii] = xyz1[ii] - xyz0[ii];
		//		vVec[ii] = xyz2[ii] - xyz0[ii];
		//		wVec[ii] = xyz3[ii] - xyz0[ii];
		//
		A[ii] = -2 * xyz1[ii] + 2 * xyz0[ii] + uderiv0[ii] + uderiv1[ii];
		B[ii] = vderiv1[ii] - vderiv0[ii];
		C[ii] = uderiv2[ii] - uderiv0[ii];
		E[ii] = -2 * xyz2[ii] + 2 * xyz0[ii] + vderiv0[ii] + vderiv2[ii];
		F[ii] = wderiv2[ii] - wderiv0[ii];
		G[ii] = vderiv3[ii] - vderiv0[ii];
		H[ii] = -2 * xyz3[ii] + 2 * xyz0[ii] + wderiv0[ii] + wderiv3[ii];
		J[ii] = uderiv3[ii] - uderiv0[ii];
		K[ii] = wderiv1[ii] - wderiv0[ii];
		// Whoever wrote the line wrap code for eclipse was smoking something they shouldn't have.
		L[ii] = (-6 * xyz0[ii] + 2 * xyz1[ii] + 2 * xyz2[ii] + 2 * xyz3[ii])
				+ (-2
						* uderiv0[ii]
																																				- uderiv1[ii]
																																				+ 0.5 * uderiv2[ii]
																																				+ 0.5 * uderiv3[ii])
				+ (-2 * vderiv0[ii] + 0.5 * vderiv1[ii] - vderiv2[ii] + 0.5
						* vderiv3[ii])
				+ (-2 * wderiv0[ii] + 0.5 * wderiv1[ii] + 0.5 * wderiv2[ii] - wderiv3[ii]);
		M[ii] = -uderiv1[ii] + 3 * xyz1[ii] - 3 * xyz0[ii] - 2 * uderiv0[ii];
		P[ii] = -vderiv2[ii] + 3 * xyz2[ii] - 3 * xyz0[ii] - 2 * vderiv0[ii];
		R[ii] = -wderiv3[ii] + 3 * xyz3[ii] - 3 * xyz0[ii] - 2 * wderiv0[ii];
		T[ii] = xyz0[ii];
		U[ii] = uderiv0[ii];
		V[ii] = vderiv0[ii];
		W[ii] = wderiv0[ii];
	}
}

void TetDivider::setupCoordMapping(const emInt verts[]) {
	for (int ii = 0; ii < 4; ii++) {
		cellVerts[ii] = verts[ii];
	}
	const double *xyz0 = m_pMesh->getCoords(verts[0]);
	const double *xyz1 = m_pMesh->getCoords(verts[1]);
	const double *xyz2 = m_pMesh->getCoords(verts[2]);
	const double *xyz3 = m_pMesh->getCoords(verts[3]);

	double len0 = getIsoLengthScale(verts[0]);
	double len1 = getIsoLengthScale(verts[1]);
	double len2 = getIsoLengthScale(verts[2]);
	double len3 = getIsoLengthScale(verts[3]);

	double vec01[] = DIFF(xyz1, xyz0);
	double vec02[] = DIFF(xyz2, xyz0);
	double vec03[] = DIFF(xyz3, xyz0);
	double vec12[] = DIFF(xyz2, xyz1);
	double vec13[] = DIFF(xyz3, xyz1);
	double vec23[] = DIFF(xyz3, xyz2);

	double ratio01 = sqrt(std::max(0.5, std::min(2., len0/len1)));
	double ratio02 = sqrt(std::max(0.5, std::min(2., len0/len2)));
	double ratio03 = sqrt(std::max(0.5, std::min(2., len0/len3)));
	double ratio12 = sqrt(std::max(0.5, std::min(2., len1/len2)));
	double ratio13 = sqrt(std::max(0.5, std::min(2., len1/len3)));
	double ratio23 = sqrt(std::max(0.5, std::min(2., len2/len3)));

	double uderiv0[3], vderiv0[3], wderiv0[3];
	double uderiv1[3], vderiv1[3], wderiv1[3];
	double uderiv2[3], vderiv2[3], wderiv2[3];
	double uderiv3[3], vderiv3[3], wderiv3[3];

	SCALE(vec01, ratio01, uderiv0);
	SCALE(vec02, ratio02, vderiv0);
	SCALE(vec03, ratio03, wderiv0);

	SCALE(vec01, 1 / ratio01, uderiv1);
	for (int ii = 0; ii < 3; ii++) {
		vderiv1[ii] = uderiv1[ii] + vec12[ii] * ratio12;
		wderiv1[ii] = uderiv1[ii] + vec13[ii] * ratio13;
	}

	SCALE(vec02, 1 / ratio02, vderiv2);
	for (int ii = 0; ii < 3; ii++) {
		uderiv2[ii] = vderiv2[ii] - vec12[ii] / ratio12;
		wderiv2[ii] = vderiv2[ii] + vec23[ii] * ratio23;
	}

	SCALE(vec03, 1 / ratio03, wderiv3);
	for (int ii = 0; ii < 3; ii++) {
		uderiv3[ii] = wderiv3[ii] - vec13[ii] / ratio13;
		vderiv3[ii] = wderiv3[ii] - vec23[ii] / ratio23;
	}

	setPolyCoeffs(xyz0, xyz1, xyz2, xyz3, uderiv0, vderiv0, wderiv0, uderiv1,
								vderiv1, wderiv1, uderiv2, vderiv2, wderiv2, uderiv3, vderiv3,
								wderiv3);
}

void TetDivider::getPhysCoordsFromParamCoords(const double uvw[3],
		double xyz[3]) {
	const double& u = uvw[0];
	const double& v = uvw[1];
	const double& w = uvw[2];
	for (int ii = 0; ii < 3; ii++) {
		xyz[ii] = xyzOffset[ii] + u * uVec[ii] + v * vVec[ii] + w * wVec[ii];
	}
	for (int ii = 0; ii < 3; ii++) {
		xyz[ii] = u
				* (U[ii] + u * (M[ii] + u * A[ii] + v * B[ii] + w * K[ii])
						+ v * w * L[ii])
							+ v * (V[ii] + v * (P[ii] + u * C[ii] + v * E[ii] + w * F[ii]))
							+ w * (W[ii] + w * (R[ii] + u * J[ii] + v * G[ii] + w * H[ii]))
							+ T[ii];
	}
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
				getPhysCoordsFromParamCoords(uvw, coords);
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
