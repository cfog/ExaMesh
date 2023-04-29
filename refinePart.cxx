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

#include "ExaMesh.h"
#include "HexDivider.h"
#include "PrismDivider.h"
#include "PyrDivider.h"
#include "TetDivider.h"
#include "BdryTriDivider.h"
#include "BdryQuadDivider.h"
#include <stdio.h>
void printLocalVerts(const exa_set<TriFaceVerts> tris, const emInt nDivs){
	for(auto itr=tris.begin(); itr!=tris.end(); itr++){
		std::cout<<"For this tri: "<< itr->getCorner(0)<<
		" "<<itr->getCorner(1)<<" "<<itr->getCorner(2)<<std::endl; 
		for (int jj = 1; jj <= nDivs - 2; jj++){
			for (int ii = 1; ii <= nDivs - 1 - jj; ii++){
				std::cout<<itr->getIntVertInd(ii,jj)<<" "; 
			} 
		}
		std::cout<<std::endl; 
	}
}
emInt subdividePartMesh(const ExaMesh *const pVM_input, 
	UMesh *const pVM_output,
	const int nDivs, const emInt partID) {
	assert(nDivs >= 1);
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
	exa_map<Edge, EdgeVerts> vertsOnEdges;
	exa_set<TriFaceVerts> vertsOnTris;
	exa_set<QuadFaceVerts> vertsOnQuads;

	// Copy vertex data into the new mesh.
	for (emInt iV = 0; iV < pVM_input->numVertsToCopy(); iV++) {
		double coords[3];
		pVM_input->getCoords(iV, coords);
		pVM_output->addVert(coords);
//		double len = pVM_input->getLengthScale(iV);
//		printf("Vert: %4d  Coords: (%8f %8f %8f)  Len: %8f\n", iV, coords[0],
//				coords[1], coords[2], len);
	}
	assert(pVM_input->numVertsToCopy() == pVM_output->numVerts());

	// Need to explicitly specify the type of mapping here.
	TetDivider TD(pVM_output, pVM_input, nDivs);
	for (emInt iT = 0; iT < pVM_input->numTets(); iT++) {
		// Divide all the edges, including storing info about which new verts
		// are on which edges
		const emInt *const thisTet = pVM_input->getTetConn(iT);
		TD.setupCoordMapping(thisTet);
		TD.createDivisionVerts(vertsOnEdges, vertsOnTris,
				vertsOnQuads);

		// And now the moment of truth:  create a flock of new tets.
		TD.createNewCells();
		if ((iT + 1) % 100000 == 0)
			fprintf(
			stderr, "Refined %'12d tets.  Tree sizes: %'12lu %'12lu %'12lu\r",
					iT + 1, vertsOnEdges.size(), vertsOnTris.size(),
					vertsOnQuads.size());
	} // Done looping over all tets
#ifndef NDEBUG
	fprintf(stderr, "\nDone with tets\n");
#endif

	PyrDivider PD(pVM_output, pVM_input, nDivs);
	for (emInt iP = 0; iP < pVM_input->numPyramids(); iP++) {
		// Divide all the edges, including storing info about which new verts
		// are on which edges
		const emInt *const thisPyr = pVM_input->getPyrConn(iP);
		PD.setupCoordMapping(thisPyr);

		PD.createDivisionVerts(vertsOnEdges, vertsOnTris,
				vertsOnQuads);

		// And now the moment of truth:  create a flock of new pyramids.
		PD.createNewCells();
		if ((iP + 1) % 100000 == 0)
			fprintf(
			stderr, "Refined %'12d pyrs.  Tree sizes: %'12lu %'12lu %'12lu\r",
					iP + 1, vertsOnEdges.size(), vertsOnTris.size(),
					vertsOnQuads.size());
	} // Done looping over all pyramids
#ifndef NDEBUG
	fprintf(stderr, "\nDone with pyramids\n");
#endif

	PrismDivider PrismD(pVM_output, pVM_input, nDivs);
	for (emInt iP = 0; iP < pVM_input->numPrisms(); iP++) {
		// Divide all the edges, including storing info about which new verts
		// are on which edges
		const emInt *const thisPrism = pVM_input->getPrismConn(iP);
		PrismD.setupCoordMapping(thisPrism);
		PrismD.createDivisionVerts(vertsOnEdges, vertsOnTris,
				vertsOnQuads);

		// And now the moment of truth:  create a flock of new prisms.
		PrismD.createNewCells();
		if ((iP + 1) % 100000 == 0)
			fprintf(
			stderr, "Refined %'12d prisms.  Tree sizes: %'12lu %'12lu %'12lu\r",
					iP + 1, vertsOnEdges.size(), vertsOnTris.size(),
					vertsOnQuads.size());
	} // Done looping over all prisms
#ifndef NDEBUG
	fprintf(stderr, "\nDone with prisms\n");
#endif

	HexDivider HD(pVM_output, pVM_input, nDivs);
	for (emInt iH = 0; iH < pVM_input->numHexes(); iH++) {
		// Divide all the edges, including storing info about which new verts
		// are on which edges
		const emInt *const thisHex = pVM_input->getHexConn(iH);
		HD.setupCoordMapping(thisHex);

		HD.createDivisionVerts(vertsOnEdges, vertsOnTris,
				vertsOnQuads);

		// And now the moment of truth:  create a flock of new hexes.
		HD.createNewCells();
		if ((iH + 1) % 100000 == 0)
			fprintf(
			stderr, "Refined %'12d hexes.  Tree sizes: %'12lu %'12lu %'12lu\r",
					iH + 1, vertsOnEdges.size(), vertsOnTris.size(),
					vertsOnQuads.size());
	} // Done looping over all hexes
#ifndef NDEBUG
	fprintf(stderr, "\nDone with hexes\n");
#endif
	exa_set<TriFaceVerts> tris=pVM_input->getTriPart(); 
	// std::cout<<"verts on Tri size: "<<vertsOnTris.size()<<
	// std::endl; 
	BdryTriDivider BTD(pVM_output, nDivs);
	for (emInt iBT = 0; iBT < pVM_input->numBdryTris(); iBT++) {
		const emInt *const thisBdryTri = pVM_input->getBdryTriConn(iBT);
		TriFaceVerts TFV (nDivs,thisBdryTri[0],thisBdryTri[1],thisBdryTri[2]); 
		
		BTD.setupCoordMapping(thisBdryTri);
		// Shouldn't need to divide anything at all here, but these function
		// copy the vertices into the CellDivider internal data structure.
		BTD.createDivisionVerts(vertsOnEdges, vertsOnTris, vertsOnQuads);

		BTD.createNewCells();
		if ((iBT + 1) % 100000 == 0)
			fprintf(
			stderr,
					"Refined %'12d bdry tris.  Tree sizes: %'12lu %'12lu %'12lu\r",
					iBT + 1, vertsOnEdges.size(), vertsOnTris.size(),
					vertsOnQuads.size());
		//BTD.getRefinedVerts(pVM_input);	
		auto it=tris.find(TFV); 
		if(it!=tris.end()){
			TFV.setPartID(it->getPartid()); 
			TFV.setRemotePartID(it->getRemotePartid()); 
			emInt remoteIndices [3]={
				it->getRemoteIndices(0),
				it->getRemoteIndices(1),
				it->getRemoteIndices(2)
			}; 
			emInt global[3]={
				it->getGlobalCorner(0),
				it->getGlobalCorner(1),
				it->getGlobalCorner(2)

			}; 
			TFV.setGlobalCorners(global[0],global[1],global[2]); 
			TFV.setRemoteIndices(remoteIndices); 
			BTD.setRefinedVerts(TFV);
			pVM_output->updateRefinedPartTris(TFV);
			// for (int ii = 0; ii <= nDivs ; ii++) {
	 		// 	for (int jj = 0; jj <= nDivs-ii ; jj++) {
			// 		std::cout<<"ii: "<<
			// 		ii<<" jj: "<<jj<<" "<<TFV.getIntVertInd(ii,jj)<<std::endl;
			// 	}
			// }	

			// std::cout<<std::endl; 
		}
			
	

	}
#ifndef NDEBUG
	fprintf(stderr, "\nDone with bdry tris\n");
#endif

	exa_set<QuadFaceVerts> quads= pVM_input->getQuadPart(); 
	BdryQuadDivider BQD(pVM_output, nDivs);
	for (emInt iBQ = 0; iBQ < pVM_input->numBdryQuads(); iBQ++) {
		const emInt *const thisBdryQuad = pVM_input->getBdryQuadConn(iBQ);
		QuadFaceVerts QFV(nDivs,thisBdryQuad[0],thisBdryQuad[1],
		thisBdryQuad[2],thisBdryQuad[3]); 
		BQD.setupCoordMapping(thisBdryQuad);

		// Shouldn't need to divide anything at all here, but this function
		// copies the triangle vertices into the CellDivider internal data structure.
		BQD.createDivisionVerts(vertsOnEdges, vertsOnTris, vertsOnQuads);

		BQD.createNewCells();
		if ((iBQ + 1) % 100000 == 0)
			fprintf(
			stderr,
					"Refined %'12d bdry quads.  Tree sizes: %'12lu %'12lu %'12lu\r",
					iBQ + 1, vertsOnEdges.size(), vertsOnTris.size(),
					vertsOnQuads.size());

		auto it=quads.find(QFV); 
		if(it!=quads.end()){
			QFV.setPartID(it->getPartid()); 
			QFV.setRemotePartID(it->getRemotePartid()); 
			emInt remoteIndices [4]={
				it->getRemoteIndices(0),
				it->getRemoteIndices(1),
				it->getRemoteIndices(2),
				it->getRemoteIndices(3)
			}; 
			emInt global[4]={
				it->getGlobalCorner(0),
				it->getGlobalCorner(1),
				it->getGlobalCorner(2),
				it->getGlobalCorner(3)

			}; 
			QFV.setGlobalCorners(global[0],global[1],global[2],global[3]); 
			QFV.setRemoteIndices(remoteIndices); 
			BQD.setRefinedVerts(QFV);
			pVM_output->updateRefinedPartQuads(QFV);

		}
	}
#ifndef NDEBUG
	fprintf(stderr, "\nDone with bdry quads\n");
#endif

//	assert(vertsOnTris.empty());
//	assert(vertsOnQuads.empty());
//
#ifndef NDEBUG
	fprintf(stderr, "Final size of edge list: %'lu\n", vertsOnEdges.size());
	fprintf(stderr, "Final size of tri list: %'lu\n", vertsOnTris.size());
	fprintf(stderr, "Final size of quad list: %'lu\n", vertsOnQuads.size());
#endif

	return pVM_output->numCells();
}

bool computeMeshSize(const struct MeshSize &MSIn, const emInt nDivs,
		struct MeshSize &MSOut) {
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
			+ MSIn.nPrisms * 2) / 2;
	ssize_t inputQuadCount = (MSIn.nBdryQuads + MSIn.nPyrs + MSIn.nPrisms * 3
			+ MSIn.nHexes * 6) / 2;
	ssize_t inputFaceCount = inputTriCount + inputQuadCount;

	ssize_t inputCellCount = MSIn.nTets + MSIn.nPyrs + MSIn.nPrisms
			+ MSIn.nHexes;
	ssize_t inputBdryEdgeCount = (MSIn.nBdryTris * 3 + MSIn.nBdryQuads * 4) / 2;
	// Upcast the first arg explicitly, and the rest should follow.
	int inputGenus = (ssize_t(MSIn.nBdryVerts) - inputBdryEdgeCount
			+ MSIn.nBdryTris + MSIn.nBdryQuads - 2) / 2;

	ssize_t inputEdges = (ssize_t(MSIn.nVerts) + inputFaceCount - inputCellCount
			- 1 - inputGenus);

	ssize_t outputFaceVerts = inputTriCount * (nDivs - 2) * (nDivs - 1) / 2
			+ inputQuadCount * (nDivs - 1) * (nDivs - 1);
	ssize_t outputCellVerts = MSIn.nTets * (nDivs - 3) * (nDivs - 2)
			* (nDivs - 1) / 6
			+ MSIn.nPyrs * (2 * nDivs - 3) * (nDivs - 2) * (nDivs - 1) / 6
			+ MSIn.nPrisms * (nDivs - 1) * (nDivs - 2) * (nDivs - 1) / 2
			+ MSIn.nHexes * (nDivs - 1) * (nDivs - 1) * (nDivs - 1);
	ssize_t outputEdgeVerts = inputEdges * (nDivs - 1);
	ssize_t outputVerts = outputFaceVerts + outputEdgeVerts + outputCellVerts
			+ MSIn.nVerts;
//	ssize_t outputEdges = outputVerts + inputFaceCount * surfFactor
//												- inputCellCount * volFactor - 1 - inputGenus;
	MSOut.nVerts = outputVerts;

	return true;
}

