/*
 * PyrDivider.cxx
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#include "GeomUtils.h"
#include "PyrDivider.h"

void PyrDivider::divideInterior(const emInt verts[]) {
  // Number of verts added:
  //    Tets:      (nD-1)(nD-2)(nD-3)/6
	const double* coords0 = m_pMesh->getCoords(verts[0]);
	const double* coords1 = m_pMesh->getCoords(verts[1]);
	const double* coords2 = m_pMesh->getCoords(verts[2]);
	const double* coords3 = m_pMesh->getCoords(verts[3]);
	const double* coords4 = m_pMesh->getCoords(verts[4]);

	double deltaInI[] = { (coords1[0] - coords0[0]), (coords1[1] - coords0[1]),
												(coords1[2] - coords0[2]) };
	double deltaInJ[] = { (coords3[0] - coords0[0]), (coords3[1] - coords0[1]),
												(coords3[2] - coords0[2]) };
	double crossDelta[] = {
			(coords2[0] + coords0[0] - coords1[0] - coords3[0]),
													(coords2[1] + coords0[1] - coords1[1] - coords3[1]),
													(coords2[2] + coords0[2] - coords1[2] - coords3[2]) };

	for (int kk = 2; kk <= nDivs - 1; kk++) {
    double zeta = 1 - double(kk) / nDivs;
    for (int jj = 1; jj <= kk - 1; jj++) {
      double eta = double(jj) / kk;
      for (int ii = 1; ii <= kk - 1; ii++) {
        double xi = double(ii) / kk;
				double coordsBase[] = { coords0[0] + deltaInI[0] * xi
																+
                                   deltaInJ[0] * eta + crossDelta[0] * xi * eta,
																coords0[1] + deltaInI[1] * xi
																+
                                   deltaInJ[1] * eta + crossDelta[1] * xi * eta,
																coords0[2] + deltaInI[2] * xi
																+
                                   deltaInJ[2] * eta +
                                   crossDelta[2] * xi * eta};

				double coords[] = { coordsBase[0] * (1 - zeta) + coords4[0] * zeta,
														coordsBase[1] * (1 - zeta) + coords4[1] * zeta,
														coordsBase[2] * (1 - zeta) + coords4[2] * zeta };
				emInt vNew = m_pMesh->addVert(coords);
				localVerts[ii][jj][kk] = vNew;
      }
    }
  } // Done looping to create all verts inside the tet.
}

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
				emInt pyr = m_pMesh->addPyramid(vertsNew);
				assert(pyr < m_pMesh->numPyrs() && pyr < m_pMesh->maxNPyrs());
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
				emInt pyr = m_pMesh->addPyramid(vertsNew);
				assert(pyr < m_pMesh->numPyrs() && pyr < m_pMesh->maxNPyrs());
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
				emInt tet = m_pMesh->addTet(vertsNew);
				assert(tet < m_pMesh->numTets() && tet < m_pMesh->maxNTets());
				assert(
						checkOrient3D(m_pMesh->getCoords(vertsNew[0]),
													m_pMesh->getCoords(vertsNew[1]),
													m_pMesh->getCoords(vertsNew[2]),
													m_pMesh->getCoords(vertsNew[3])));
      }
    }

    // The set on lines of constant i on level l.
    for (int jj = 0; jj <= level - 1; jj++) {
      for (int ii = 1; ii <= level - 1; ii++) {
				emInt vertsNew[] = { localVerts[ii][jj][level],
															localVerts[ii][jj + 1][level],
															localVerts[ii - 1][jj][level - 1],
															localVerts[ii][jj][level - 1] };
				emInt tet = m_pMesh->addTet(vertsNew);
				assert(tet < m_pMesh->numTets() && tet < m_pMesh->maxNTets());
				assert(
						checkOrient3D(m_pMesh->getCoords(vertsNew[0]),
													m_pMesh->getCoords(vertsNew[1]),
													m_pMesh->getCoords(vertsNew[2]),
													m_pMesh->getCoords(vertsNew[3])));
      }
    }

  } // Done with this level
}
