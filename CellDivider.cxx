/*
 * CellDivider.cxx
 *
 *  Created on: Jun. 2, 2019
 *      Author: cfog
 */

#include "examesh.h"

#include <map>
#include <set>

#include "CellDivider.h"

void getEdgeVerts(UMesh *pVM, std::map<Edge, EdgeVerts> &vertsOnEdges,
		const int nDivs, const emInt v0, const emInt v1, EdgeVerts &EV) {
	typename std::map<Edge, EdgeVerts>::iterator iterEdges;
	Edge E(v0, v1);
	iterEdges = vertsOnEdges.find(E);

	if (iterEdges == vertsOnEdges.end()) {
		// Doesn't exist yet, so create it.
		EV.verts[0] = E.getV0();
		EV.verts[nDivs] = E.getV1();
		const double *startLoc = pVM->getCoords(EV.verts[0]);
		const double *endLoc = pVM->getCoords(EV.verts[nDivs]);
		// TODO:  Coords shouldn't be equally spaced, but smoothly placed
		// based on length scale at the ends.
		double delta[] = { (endLoc[0] - startLoc[0]) / nDivs, (endLoc[1]
				- startLoc[1])
																													/ nDivs,
												(endLoc[2] - startLoc[2]) / nDivs };
		for (int ii = 1; ii < nDivs; ii++) {
			double newCoords[] = { startLoc[0] + ii * delta[0], startLoc[1]
					+ ii * delta[1],
															startLoc[2] + ii * delta[2] };
			EV.verts[ii] = pVM->addVert(newCoords);
		}
		vertsOnEdges.insert(std::make_pair(E, EV));
	}
	else {
		EV = iterEdges->second;
	}
}

TriFaceVerts::TriFaceVerts(const emInt v0, const emInt v1, const emInt v2) {
	corners[0] = v0;
	corners[1] = v1;
	corners[2] = v2;
}

QuadFaceVerts::QuadFaceVerts(const emInt v0, const emInt v1, const emInt v2,
		const emInt v3) {
	corners[0] = v0;
	corners[1] = v1;
	corners[2] = v2;
	corners[3] = v3;
}

void sortVerts3(const emInt input[3], emInt output[3]) {
	// This is insertion sort, specialized for three inputs.
	if (input[1] < input[0]) {
		output[0] = input[1];
		output[1] = input[0];
	}
	else {
		output[0] = input[0];
		output[1] = input[1];
	}
	if (input[2] < output[1]) {
		output[2] = output[1];
		if (input[2] < output[0]) {
			output[1] = output[0];
			output[0] = input[2];
		}
		else {
			output[1] = input[2];
		}
	}
	else {
		output[2] = input[2];
	}
}

void sortVerts4(const emInt input[4], emInt output[4]) {
	// This is insertion sort, specialized for four inputs.
	if (input[1] < input[0]) {
		output[0] = input[1];
		output[1] = input[0];
	}
	else {
		output[0] = input[0];
		output[1] = input[1];
	}

	if (input[2] < output[1]) {
		output[2] = output[1];
		if (input[2] < output[0]) {
			output[1] = output[0];
			output[0] = input[2];
		}
		else {
			output[1] = input[2];
		}
	}
	else {
		output[2] = input[2];
	}

	if (input[3] < output[2]) {
		output[3] = output[2];
		if (input[3] < output[1]) {
			output[2] = output[1];
			if (input[3] < output[0]) {
				output[1] = output[0];
				output[0] = input[3];
			}
			else {
				output[1] = input[3];
			}
		}
		else {
			output[2] = input[3];
		}
	}
	else {
		output[3] = input[3];
	}
}

bool operator<(const TriFaceVerts& a, const TriFaceVerts& b) {
	emInt aSorted[3], bSorted[3];
	sortVerts3(a.corners, aSorted);
	sortVerts3(b.corners, bSorted);
	return (aSorted[0] < bSorted[0]
			|| (aSorted[0] == bSorted[0] && aSorted[1] < bSorted[1])
			|| (aSorted[0] == bSorted[0] && aSorted[1] == bSorted[1]
					&& aSorted[2] < bSorted[2]));
}

bool operator<(const QuadFaceVerts& a, const QuadFaceVerts& b) {
	emInt aSorted[4], bSorted[4];
	sortVerts4(a.corners, aSorted);
	sortVerts4(b.corners, bSorted);
	return (aSorted[0] < bSorted[0]
			|| (aSorted[0] == bSorted[0] && aSorted[1] < bSorted[1])
			|| (aSorted[0] == bSorted[0] && aSorted[1] == bSorted[1]
					&& aSorted[2] < bSorted[2])
			|| (aSorted[0] == bSorted[0] && aSorted[1] == bSorted[1]
					&& aSorted[2] == bSorted[2] && aSorted[3] < bSorted[3]));
}

void getTriVerts(UMesh *pVM, std::set<TriFaceVerts> &vertsOnTris,
		const int nDivs, const emInt vert0, const emInt vert1, const emInt vert2,
		TriFaceVerts &TFV) {
	typename std::set<TriFaceVerts>::iterator iterTris;
	TriFaceVerts TFVTemp(vert0, vert1, vert2);

	iterTris = vertsOnTris.find(TFVTemp);
	if (iterTris == vertsOnTris.end()) {
		const double *coords0 = pVM->getCoords(vert0);
		const double *coords1 = pVM->getCoords(vert1);
		const double *coords2 = pVM->getCoords(vert2);

		double deltaInI[] = { (coords1[0] - coords0[0]) / nDivs, (coords1[1]
				- coords0[1])
																															/ nDivs,
													(coords1[2] - coords0[2]) / nDivs };
		double deltaInJ[] = { (coords2[0] - coords0[0]) / nDivs, (coords2[1]
				- coords0[1])
																															/ nDivs,
													(coords2[2] - coords0[2]) / nDivs };
		TFV.corners[0] = vert0;
		TFV.corners[1] = vert1;
		TFV.corners[2] = vert2;

		for (int jj = 0; jj < nDivs - 2; jj++) {
			for (int ii = 0; ii < nDivs - 2 - jj; ii++) {
				double newCoords[] = { coords0[0] + deltaInI[0] * (ii + 1)
																+ deltaInJ[0] * (jj + 1),
																coords0[1] + deltaInI[1] * (ii + 1)
																+ deltaInJ[1] * (jj + 1),
																coords0[2] + deltaInI[2] * (ii + 1)
																+ deltaInJ[2] * (jj + 1) };
				emInt vNew = pVM->addVert(newCoords);
				TFV.intVerts[ii][jj] = vNew;
			}
		} // Done looping over all interior verts for the triangle.
		vertsOnTris.insert(TFV);
	}
	else {
		TFV = *iterTris;
		vertsOnTris.erase(iterTris); // Will never need this again.
	}
}

void getQuadVerts(UMesh *pVM, std::set<QuadFaceVerts> &vertsOnQuads,
		const int nDivs, const emInt vert0, const emInt vert1, const emInt vert2,
		const emInt vert3, QuadFaceVerts &QFV) {
	typename std::set<QuadFaceVerts>::iterator iterQuads;
	QuadFaceVerts QFVTemp(vert0, vert1, vert2, vert3);

	iterQuads = vertsOnQuads.find(QFVTemp);
	if (iterQuads == vertsOnQuads.end()) {
		const double *coords0 = pVM->getCoords(vert0);
		const double *coords1 = pVM->getCoords(vert1);
		const double *coords2 = pVM->getCoords(vert2);
		const double *coords3 = pVM->getCoords(vert3);

		double deltaInI[] = { (coords1[0] - coords0[0]) / nDivs, (coords1[1]
				- coords0[1])
																															/ nDivs,
													(coords1[2] - coords0[2]) / nDivs };
		double deltaInJ[] = { (coords3[0] - coords0[0]) / nDivs, (coords3[1]
				- coords0[1])
																															/ nDivs,
													(coords3[2] - coords0[2]) / nDivs };
		double crossDelta[] = {
				(coords2[0] + coords0[0] - coords1[0] - coords3[0]) / (nDivs * nDivs),
				(coords2[1] + coords0[1] - coords1[1] - coords3[1]) / (nDivs * nDivs),
				(coords2[2] + coords0[2] - coords1[2] - coords3[2]) / (nDivs * nDivs) };
		QFV.corners[0] = vert0;
		QFV.corners[1] = vert1;
		QFV.corners[2] = vert2;
		QFV.corners[3] = vert3;

		for (int jj = 1; jj <= nDivs - 1; jj++) {
			for (int ii = 1; ii <= nDivs - 1; ii++) {
				double newCoords[] = { coords0[0] + deltaInI[0] * ii + deltaInJ[0] * jj
																+ crossDelta[0] * ii * jj,
																coords0[1] + deltaInI[1] * ii + deltaInJ[1] * jj
																+ crossDelta[1] * ii * jj,
																coords0[2] + deltaInI[2] * ii + deltaInJ[2] * jj
																+ crossDelta[2] * ii * jj };
				emInt vNew = pVM->addVert(newCoords);
				QFV.intVerts[ii - 1][jj - 1] = vNew;
			}
		} // Done looping over all interior verts for the triangle.
		vertsOnQuads.insert(QFV);
	}
	else {
	QFV = *iterQuads;
		vertsOnQuads.erase(iterQuads); // Will never need this again.
	}
}

void CellDivider::divideEdges(std::map<Edge, EdgeVerts> &vertsOnEdges,
		const emInt verts[]) {
	// Divide all the edges, including storing info about which new verts
	// are on which edges
	for (int iE = 0; iE < numEdges; iE++) {

		EdgeVerts EV;
		getEdgeVerts(m_pMesh, vertsOnEdges, nDivs, verts[edgeVertIndices[iE][0]],
									verts[edgeVertIndices[iE][1]], EV);

		// Now transcribe these into the master table for this cell.
		emInt startIndex = 1000, endIndex = 1000;
		if (EV.verts[0] == verts[edgeVertIndices[iE][0]]) {
			// Transcribe this edge forward.
			startIndex = edgeVertIndices[iE][0];
			endIndex = edgeVertIndices[iE][1];
		}
		else {
			startIndex = edgeVertIndices[iE][1];
			endIndex = edgeVertIndices[iE][0];
		}
		int startI = vertIJK[startIndex][0];
		int startJ = vertIJK[startIndex][1];
		int startK = vertIJK[startIndex][2];
		int incrI = (vertIJK[endIndex][0] - startI) / nDivs;
		int incrJ = (vertIJK[endIndex][1] - startJ) / nDivs;
		int incrK = (vertIJK[endIndex][2] - startK) / nDivs;

		for (int ii = 0; ii <= nDivs; ii++) {
			int II = startI + ii * incrI;
			int JJ = startJ + ii * incrJ;
			int KK = startK + ii * incrK;
			assert(II >= 0 && II <= nDivs);
			assert(JJ >= 0 && JJ <= nDivs);
			assert(KK >= 0 && KK <= nDivs);
			localVerts[II][JJ][KK] = EV.verts[ii];
		}
	}
}

void CellDivider::divideFaces(std::set<TriFaceVerts> &vertsOnTris,
		std::set<QuadFaceVerts> &vertsOnQuads, const emInt verts[]) {
	// Divide all the faces, including storing info about which new verts
	// are on which faces

	// The quad faces are first.
	for (int iF = 0; iF < numQuadFaces; iF++) {
		QuadFaceVerts QFV;
		const emInt vert0 = verts[faceVertIndices[iF][0]];
		const emInt vert1 = verts[faceVertIndices[iF][1]];
		const emInt vert2 = verts[faceVertIndices[iF][2]];
		const emInt vert3 = verts[faceVertIndices[iF][3]];
		getQuadVerts(m_pMesh, vertsOnQuads, nDivs, vert0, vert1, vert2, vert3, QFV);
		// Now extract info from the TFV and stuff it into the cell's point
		// array.

		// Critical first step: identify which vert is which.
		emInt corner[] = { 1000, 1000, 1000, 1000 };

		for (int iC = 0; iC < 4; iC++) {
			const emInt corn = QFV.corners[iC];
			for (int iV = 0; iV < numVerts; iV++) {
				const emInt cand = verts[iV];
				if (corn == cand) {
					corner[iC] = iV;
					break;
				}
			}
		}

		int startI = vertIJK[corner[0]][0];
		int startJ = vertIJK[corner[0]][1];
		int startK = vertIJK[corner[0]][2];
		int incrIi = (vertIJK[corner[1]][0] - startI) / nDivs;
		int incrJi = (vertIJK[corner[1]][1] - startJ) / nDivs;
		int incrKi = (vertIJK[corner[1]][2] - startK) / nDivs;
		int incrIj = (vertIJK[corner[3]][0] - startI) / nDivs;
		int incrJj = (vertIJK[corner[3]][1] - startJ) / nDivs;
		int incrKj = (vertIJK[corner[3]][2] - startK) / nDivs;

		for (int jj = 1; jj <= nDivs - 1; jj++) {
			for (int ii = 1; ii <= nDivs - 1; ii++) {
				int II = startI + incrIi * ii + incrIj * jj;
				int JJ = startJ + incrJi * ii + incrJj * jj;
				int KK = startK + incrKi * ii + incrKj * jj;
				assert(II >= 0 && II <= nDivs);
				assert(JJ >= 0 && JJ <= nDivs);
				assert(KK >= 0 && KK <= nDivs);

				localVerts[II][JJ][KK] = QFV.intVerts[ii - 1][jj - 1];
			}
		}
	}

	for (int iF = numQuadFaces; iF < numQuadFaces + numTriFaces; iF++) {
		TriFaceVerts TFV;
		const emInt vert0 = verts[faceVertIndices[iF][0]];
		const emInt vert1 = verts[faceVertIndices[iF][1]];
		const emInt vert2 = verts[faceVertIndices[iF][2]];
		getTriVerts(m_pMesh, vertsOnTris, nDivs, vert0, vert1, vert2, TFV);
		// Now extract info from the TFV and stuff it into the Prismamid's point
		// array.

		// 1000 is way more points than cells have.
		emInt corner[] = { 1000, 1000, 1000 };
		// Critical first step: identify which vert is which.
		for (int iC = 0; iC < 3; iC++) {
			const emInt corn = TFV.corners[iC];
			for (int iV = 0; iV < numVerts; iV++) {
				const emInt cand = verts[iV];
				if (corn == cand) {
					corner[iC] = iV;
					break;
				}
			}
		}

		int startI = vertIJK[corner[0]][0];
		int startJ = vertIJK[corner[0]][1];
		int startK = vertIJK[corner[0]][2];
		int incrIi = (vertIJK[corner[1]][0] - startI) / nDivs;
		int incrJi = (vertIJK[corner[1]][1] - startJ) / nDivs;
		int incrKi = (vertIJK[corner[1]][2] - startK) / nDivs;
		int incrIj = (vertIJK[corner[2]][0] - startI) / nDivs;
		int incrJj = (vertIJK[corner[2]][1] - startJ) / nDivs;
		int incrKj = (vertIJK[corner[2]][2] - startK) / nDivs;

		for (int jj = 0; jj < nDivs - 2; jj++) {
			for (int ii = 0; ii < nDivs - 2 - jj; ii++) {
				int II = startI + incrIi * (ii + 1) + incrIj * (jj + 1);
				int JJ = startJ + incrJi * (ii + 1) + incrJj * (jj + 1);
				int KK = startK + incrKi * (ii + 1) + incrKj * (jj + 1);
				assert(II >= 0 && II <= nDivs);
				assert(JJ >= 0 && JJ <= nDivs);
				assert(KK >= 0 && KK <= nDivs);

				localVerts[II][JJ][KK] = TFV.intVerts[ii][jj];
			}
		}
	}
}





