/*
 * PrismDivider.cxx
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#include "HexDivider.h"

void HexDivider::divideInterior(const emInt verts[]) {
  // Number of verts added:
  //    Tets:      (nD-1)(nD-2)(nD-3)/6
	const double* coords0 = m_pMesh->getCoords(verts[0]);
	const double* coords1 = m_pMesh->getCoords(verts[1]);
	const double* coords2 = m_pMesh->getCoords(verts[2]);
	const double* coords3 = m_pMesh->getCoords(verts[3]);
	const double* coords4 = m_pMesh->getCoords(verts[4]);
	const double* coords5 = m_pMesh->getCoords(verts[5]);
	const double* coords6 = m_pMesh->getCoords(verts[6]);
	const double* coords7 = m_pMesh->getCoords(verts[7]);

	double deltaInXiBot[] = { coords1[0] - coords0[0], coords1[1] - coords0[1],
														coords1[2] - coords0[2] };
	double deltaInEtaBot[] = { coords3[0] - coords0[0], coords3[1] - coords0[1],
															coords3[2] - coords0[2] };
	double deltaCrossBot[] =
			{ coords2[0] + coords0[0] - coords1[0] - coords3[0], coords2[1]
					+ coords0[1]
																														- coords1[1]
																														- coords3[1],
				coords2[2] + coords0[2] - coords1[2] - coords3[2] };

	double deltaInXiTop[] = { coords5[0] - coords4[0], coords5[1] - coords4[1],
														coords5[2] - coords4[2] };
	double deltaInEtaTop[] = { coords7[0] - coords4[0], coords7[1] - coords4[1],
															coords7[2] - coords4[2] };
	double deltaCrossTop[] =
			{ coords6[0] + coords4[0] - coords5[0] - coords7[0], coords6[1]
					+ coords4[1]
																														- coords5[1]
																														- coords7[1],
				coords6[2] + coords4[2] - coords5[2] - coords7[2] };

	for (int kk = 1; kk <= nDivs - 1; kk++) {
    double zeta = (1 - double(kk) / nDivs);
		for (int jj = 1; jj <= nDivs - 1; jj++) {
      double eta = double(jj) / nDivs;
			for (int ii = 1; ii <= nDivs - 1; ii++) {
        double xi = double(ii) / nDivs;
        double coordsBot[] = {
            coords0[0] + deltaInXiBot[0] * xi
																+ deltaInEtaBot[0] * eta
																+
                deltaCrossBot[0] * xi * eta,
																coords0[1] + deltaInXiBot[1] * xi
																+ deltaInEtaBot[1] * eta
																+
                deltaCrossBot[1] * xi * eta,
																coords0[2] + deltaInXiBot[2] * xi
																+ deltaInEtaBot[2] * eta
																+
                deltaCrossBot[2] * xi * eta};
        double coordsTop[] = {
            coords4[0] + deltaInXiTop[0] * xi
																+ deltaInEtaTop[0] * eta
																+
                deltaCrossTop[0] * xi * eta,
																coords4[1] + deltaInXiTop[1] * xi
																+ deltaInEtaTop[1] * eta
																+
                deltaCrossTop[1] * xi * eta,
																coords4[2] + deltaInXiTop[2] * xi
																+ deltaInEtaTop[2] * eta
																+
                deltaCrossTop[2] * xi * eta};
				double coordsNew[] =
						{ coordsBot[0] * (1 - zeta) + coordsTop[0] * zeta,
							coordsBot[1] * (1 - zeta) + coordsTop[1] * zeta,
							coordsBot[2] * (1 - zeta) + coordsTop[2] * zeta };
				emInt vNew = m_pMesh->addVert(coordsNew);
				localVerts[ii][jj][kk] = vNew;
      }
    } // Done looping over all interior verts for the triangle.
  }   // Done looping over all levels for the prism.
}

void HexDivider::createNewCells() {
	// Output info about the points for this hex, layer by layer

//	for (int level = 0; level <= nDivs; level++) {
//		fprintf(stderr, "Level = %d\n", level);
//		for (int jj = 0; jj <= level; jj++) {
//			for (int ii = 0; ii <= level; ii++) {
//				const double* coords = m_pMesh->getCoords(localVerts[ii][jj][level]);
//				fprintf(stderr, "  i = %2d, j = %2d:  (%10.6f, %10.6f, %10.6f)\n", ii,
//								jj, coords[0], coords[1], coords[2]);
//      }
//    }
//  }

	for (int level = 1; level <= nDivs; level++) {
//		fprintf(stderr, "Level: %d\n", level);
		// Create new hexes.  Always (nDivs-1)^2 for each level.
		for (int jj = 0; jj <= nDivs - 1; jj++) {
			for (int ii = 0; ii <= nDivs - 1; ii++) {
				emInt vertsNew[] = { localVerts[ii][jj][level],
															localVerts[ii + 1][jj][level],
															localVerts[ii + 1][jj + 1][level],
															localVerts[ii][jj + 1][level],
															localVerts[ii][jj][level - 1],
															localVerts[ii + 1][jj][level - 1], localVerts[ii
																	+ 1][jj + 1][level - 1],
															localVerts[ii][jj + 1][level - 1] };
				emInt hex = m_pMesh->addHex(vertsNew);
				assert(hex < m_pMesh->maxNHexes() && hex < m_pMesh->numHexes());
      }
    } // Done with this row (constant j)
  }   // Done with this level
}
