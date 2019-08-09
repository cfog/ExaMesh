//////////////////////////////////////////////////////////////////////////
//
// Refine a mesh (a coarse mesh on a single part in a parallel mesh)
// by smooth subdivision, and write that mesh to a file.
//
//////////////////////////////////////////////////////////////////////////

#include <string.h>
#include <locale.h>
#include <unistd.h>
#include <time.h>

#include <set>

#include "HexDivider.h"
#include "PrismDivider.h"
#include "PyrDivider.h"
#include "TetDivider.h"
#include "BdryTriDivider.h"
#include "BdryQuadDivider.h"
#include "examesh.h"

emInt subdividePartMesh(const UMesh * const pVM_input,
		UMesh * const pVM_output, const int nDivs) {
  assert(nDivs >= 2);
  // Assumption:  the mesh is already ordered in a way that seems sensible
  // to the caller, both cells and vertices.  As a result, we can create new
  // verts and cells on the fly, with the expectation that the new ones will
  // remain about as well-ordered as the old ones.

  // Base of the tet is tri 012.
  //   Pt 0 is at indices (i,j,k) = (0,0,nDivs).
  //   Pt 1 is at indices (i,j,k) = (nDivs,0,nDivs).
  //   Pt 2 is at indices (i,j,k) = (0,nDivs,nDivs).
  //   Pt 3 is at indices (i,j,k) = (0,0,0).
  // Other created points are laid into this framework.
  // The i,j,k system isn't topologically right-handed.  It's set up so that k
  // corresponds to layers in the divided tet, with each layer having
  // progressively more points / tris.  Nevertheless, the tets produces should
  // be geometrically right-handed.

  // TODO: Potentially, identify in advance how many times each edge is used,
  // so that when all of them have appeared, the edge can be removed from the
  // map.
	std::map<Edge, EdgeVerts> vertsOnEdges;
	std::set<TriFaceVerts> vertsOnTris;
	std::set<QuadFaceVerts> vertsOnQuads;

	// Copy vertex data into the new mesh.
	for (emInt iV = 0; iV < pVM_input->numVerts(); iV++) {
		pVM_output->addVert(pVM_input->getCoords(iV));
	}
	assert(pVM_input->numVerts() == pVM_output->numVerts());

	TetDivider TD(pVM_output, nDivs);
  // TODO:  Factor edge handling, face handling, and overall tet handling
	for (emInt iT = 0; iT < pVM_input->numTets(); iT++) {
    // Divide all the edges, including storing info about which new verts
    // are on which edges
		const emInt* const thisTet = pVM_input->getTetConn(iT);
		TD.setupCoordMapping(thisTet);
		TD.divideEdges(vertsOnEdges, thisTet);

    // Divide all the faces, including storing info about which new verts
    // are on which faces
		TD.divideFaces(vertsOnTris, vertsOnQuads, thisTet);

    // Divide the cell
    if (nDivs > 3) {
			TD.divideInterior();
    } // Done with internal division

		// And now the moment of truth:  create a flock of new tets.
    TD.createNewCells();
  } // Done looping over all tets

	PyrDivider PD(pVM_output, nDivs);
  // TODO:  Factor edge handling, face handling, and overall tet handling
	for (emInt iP = 0; iP < pVM_input->numPyrs(); iP++) {
    // Divide all the edges, including storing info about which new verts
    // are on which edges
		const emInt* const thisPyr = pVM_input->getPyrConn(iP);
		PD.setupCoordMapping(thisPyr);

		PD.divideEdges(vertsOnEdges, thisPyr);

    // Divide all the faces, including storing info about which new verts
    // are on which faces
		PD.divideFaces(vertsOnTris, vertsOnQuads, thisPyr);

    // Divide the cell
    if (nDivs >= 3) {
			PD.divideInterior();
    } // Done with internal division

		// And now the moment of truth:  create a flock of new pyramids.
    PD.createNewCells();
  } // Done looping over all pyramids

	PrismDivider PrismD(pVM_output, nDivs);
  // TODO:  Factor edge handling, face handling, and overall tet handling
	for (emInt iP = 0; iP < pVM_input->numPrisms(); iP++) {
    // Divide all the edges, including storing info about which new verts
    // are on which edges
		const emInt* const thisPrism = pVM_input->getPrismConn(iP);
		PrismD.setupCoordMapping(thisPrism);

		PrismD.divideEdges(vertsOnEdges, thisPrism);

    // Divide all the faces, including storing info about which new verts
    // are on which faces
		PrismD.divideFaces(vertsOnTris, vertsOnQuads, thisPrism);

    // Divide the cell
    if (nDivs >= 3) {
			PrismD.divideInterior();
    } // Done with internal division

		// And now the moment of truth:  create a flock of new prisms.
    PrismD.createNewCells();
  } // Done looping over all tets

	HexDivider HD(pVM_output, nDivs);
  // TODO:  Factor edge handling, face handling, and overall tet handling
	for (emInt iH = 0; iH < pVM_input->numHexes(); iH++) {
    // Divide all the edges, including storing info about which new verts
    // are on which edges
		const emInt* const thisHex = pVM_input->getHexConn(iH);
		HD.setupCoordMapping(thisHex);

		HD.divideEdges(vertsOnEdges, thisHex);

    // Divide all the faces, including storing info about which new verts
    // are on which faces
		HD.divideFaces(vertsOnTris, vertsOnQuads, thisHex);

    // Divide the cell
    if (nDivs >= 2) {
			HD.divideInterior();
    } // Done with internal division

		// And now the moment of truth:  create a flock of new hexes.
    HD.createNewCells();
  } // Done looping over all pyramids

	BdryTriDivider BTD(pVM_output, nDivs);
	for (emInt iBT = 0; iBT < pVM_input->numBdryTris(); iBT++) {
		const emInt* const thisBdryTri = pVM_input->getBdryTriConn(iBT);

		// Shouldn't need to divide anything at all here, but these function
		// copy the vertices into the CellDivider internal data structure.
		BTD.divideEdges(vertsOnEdges, thisBdryTri);
		BTD.divideFaces(vertsOnTris, vertsOnQuads, thisBdryTri);

		BTD.createNewCells();
	}

	BdryQuadDivider BQD(pVM_output, nDivs);
	for (emInt iBT = 0; iBT < pVM_input->numBdryQuads(); iBT++) {
		const emInt* const thisBdryQuad = pVM_input->getBdryQuadConn(iBT);

		// Shouldn't need to divide anything at all here, but this function
		// copies the triangle vertices into the CellDivider internal data structure.
		BQD.divideEdges(vertsOnEdges, thisBdryQuad);
		BQD.divideFaces(vertsOnTris, vertsOnQuads, thisBdryQuad);

		BQD.createNewCells();
	}

//	assert(vertsOnTris.empty());
//	assert(vertsOnQuads.empty());
//
	fprintf(stderr, "Final size of edge list: %'lu\n",
							vertsOnEdges.size());

	return pVM_output->numCells();
}

bool computeMeshSize(const struct MeshSize& MSIn, const emInt nDivs,
		struct MeshSize& MSOut) {
	// It's relatively easy to compute some of these quantities:
	const emInt surfFactor = nDivs * nDivs;
	const emInt volFactor = surfFactor * nDivs;

	const emInt maxFaces = EMINT_MAX / surfFactor;
	const emInt maxCells = EMINT_MAX / volFactor;

	if (MSIn.nBdryTris > maxFaces || MSIn.nBdryQuads > maxFaces
			|| MSIn.nTets > maxCells || MSIn.nPyrs > maxCells
			|| MSIn.nPrisms > maxCells || MSIn.nHexes > maxCells) {
		fprintf(stderr, "Output mesh will exceed max index size!\n");
		return false;
	}

	MSOut.nBdryTris = MSIn.nBdryTris * surfFactor;
	MSOut.nBdryQuads = MSIn.nBdryQuads * surfFactor;

	MSOut.nTets = MSIn.nTets * volFactor
			+ MSIn.nPyrs * (volFactor - nDivs) * 2 / 3;
	MSOut.nPyrs = MSIn.nPyrs * (2 * volFactor + nDivs) / 3;
	MSOut.nPrisms = MSIn.nPrisms * volFactor;
	MSOut.nHexes = MSIn.nHexes * volFactor;

	// Use signed 64-bit ints for these calculations.  It's possible someone will ask for
	// something that blows out 32-bit unsigned ints, and will need to be stopped.
	ssize_t inputTriCount = (MSIn.nBdryTris + MSIn.nTets * 4 + MSIn.nPyrs * 4
														+ MSIn.nPrisms * 2)
													/ 2;
	ssize_t inputQuadCount = (MSIn.nBdryQuads + MSIn.nPyrs + MSIn.nPrisms * 3
														+ MSIn.nHexes * 6)
														/ 2;
	ssize_t inputFaceCount = inputTriCount + inputQuadCount;

	ssize_t inputCellCount = MSIn.nTets + MSIn.nPyrs + MSIn.nPrisms + MSIn.nHexes;
	ssize_t inputBdryEdgeCount = (MSIn.nBdryTris * 3 + MSIn.nBdryQuads * 4) / 2;
	// Upcast the first arg explicitly, and the rest should follow.
	int inputGenus = (ssize_t(MSIn.nBdryVerts) - inputBdryEdgeCount
			+ MSIn.nBdryTris
										+ MSIn.nBdryQuads
										- 2)
										/ 2;

	ssize_t inputEdges = (ssize_t(MSIn.nVerts) + inputFaceCount - inputCellCount
												- 1 - inputGenus);

	ssize_t outputFaceVerts = inputTriCount * (nDivs - 2) * (nDivs - 1) / 2
			+ inputQuadCount * (nDivs - 1) * (nDivs - 1);
	ssize_t outputCellVerts = MSIn.nTets * (nDivs - 3) * (nDivs - 2) * (nDivs - 1)
			/ 6
														+ MSIn.nPyrs * (2 * nDivs - 3) * (nDivs - 2)
															* (nDivs - 1)
															/ 6
														+ MSIn.nPrisms * (nDivs - 1) * (nDivs - 2)
															* (nDivs - 1)
															/ 2
														+ MSIn.nHexes * (nDivs - 1) * (nDivs - 1)
															* (nDivs - 1);
	ssize_t outputEdgeVerts = inputEdges * (nDivs - 1);
	ssize_t outputVerts = outputFaceVerts + outputEdgeVerts + outputCellVerts
												+ MSIn.nVerts;
//	ssize_t outputEdges = outputVerts + inputFaceCount * surfFactor
//												- inputCellCount * volFactor - 1 - inputGenus;
	MSOut.nVerts = outputVerts;

	return true;
}

