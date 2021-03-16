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
 * CellDivider.cxx
 *
 *  Created on: Jun. 2, 2019
 *      Author: cfog
 */

#include "ExaMesh.h"
#include "GeomUtils.h"
#include "CellDivider.h"

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

TriFaceVerts::TriFaceVerts(const emInt v0, const emInt v1, const emInt v2,
		const emInt type, const emInt elemInd) {
	corners[0] = v0;
	corners[1] = v1;
	corners[2] = v2;
	volElement = elemInd;
	volElementType = type;
	intVerts = nullptr;
	setupSorted();
}

void TriFaceVerts::setupSorted() {
	sortVerts3(corners, sorted);
}

bool operator==(const TriFaceVerts& a, const TriFaceVerts& b) {
	return (a.sorted[0] == b.sorted[0] && a.sorted[1] == b.sorted[1]
					&& a.sorted[2] == b.sorted[2]);
}

bool operator<(const TriFaceVerts& a, const TriFaceVerts& b) {
	return (a.sorted[0] < b.sorted[0]
			|| (a.sorted[0] == b.sorted[0] && a.sorted[1] < b.sorted[1])
			|| (a.sorted[0] == b.sorted[0] && a.sorted[1] == b.sorted[1]
					&& a.sorted[2] < b.sorted[2]));
}

QuadFaceVerts::QuadFaceVerts(const emInt v0, const emInt v1, const emInt v2,
		const emInt v3, const emInt type, const emInt elemInd) {
	corners[0] = v0;
	corners[1] = v1;
	corners[2] = v2;
	corners[3] = v3;
	volElement = elemInd;
	volElementType = type;
	setupSorted();
}

void QuadFaceVerts::setupSorted() {
	sortVerts4(corners, sorted);
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

bool operator==(const QuadFaceVerts& a, const QuadFaceVerts& b) {
	return (a.sorted[0] == b.sorted[0] && a.sorted[1] == b.sorted[1]
					&& a.sorted[2] == b.sorted[2] && a.sorted[3] == b.sorted[3]);
}

bool operator<(const QuadFaceVerts& a, const QuadFaceVerts& b) {
	return (a.sorted[0] < b.sorted[0]
			|| (a.sorted[0] == b.sorted[0] && a.sorted[1] < b.sorted[1])
			|| (a.sorted[0] == b.sorted[0] && a.sorted[1] == b.sorted[1]
					&& a.sorted[2] < b.sorted[2])
			|| (a.sorted[0] == b.sorted[0] && a.sorted[1] == b.sorted[1]
					&& a.sorted[2] == b.sorted[2] && a.sorted[3] < b.sorted[3]));
}

int CellDivider::checkOrient3D(const emInt verts[4]) const {
	double coords0[3], coords1[3], coords2[3], coords3[3];
	m_pMesh->getCoords(verts[0], coords0);
	m_pMesh->getCoords(verts[1], coords1);
	m_pMesh->getCoords(verts[2], coords2);
	m_pMesh->getCoords(verts[3], coords3);
	return ::checkOrient3D(coords0, coords1, coords2, coords3);
}

void CellDivider::getEdgeVerts(exa_map<Edge, EdgeVerts> &vertsOnEdges,
		const int edge, const double dihedral, EdgeVerts &EV) {
	int ind0 = edgeVertIndices[edge][0];
	int ind1 = edgeVertIndices[edge][1];

	emInt vert0 = cellVerts[ind0];
	emInt vert1 = cellVerts[ind1];

	Edge E(vert0, vert1);
	auto iterEdges = vertsOnEdges.find(E);

	if (iterEdges == vertsOnEdges.end()) {
		// Doesn't exist yet, so create it.
		EV.verts[0] = E.getV0();
		EV.verts[nDivs] = E.getV1();
		EV.m_totalDihed = dihedral;

		bool forward = true;
		if (EV.verts[0] != vert0) {
			forward = false;
		}

		double uvwStart[3], uvwEnd[3];
		if (forward) {
			uvwStart[0] = uvwIJK[ind0][0];
			uvwStart[1] = uvwIJK[ind0][1];
			uvwStart[2] = uvwIJK[ind0][2];

			uvwEnd[0] = uvwIJK[ind1][0];
			uvwEnd[1] = uvwIJK[ind1][1];
			uvwEnd[2] = uvwIJK[ind1][2];
		}
		else {
			uvwStart[0] = uvwIJK[ind1][0];
			uvwStart[1] = uvwIJK[ind1][1];
			uvwStart[2] = uvwIJK[ind1][2];

			uvwEnd[0] = uvwIJK[ind0][0];
			uvwEnd[1] = uvwIJK[ind0][1];
			uvwEnd[2] = uvwIJK[ind0][2];
		}
		double delta[] = { (uvwEnd[0] - uvwStart[0]) / nDivs, (uvwEnd[1]
				- uvwStart[1])
																													/ nDivs,
												(uvwEnd[2] - uvwStart[2]) / nDivs };
		for (int ii = 1; ii < nDivs; ii++) {
			double uvw[] = { uvwStart[0] + ii * delta[0], uvwStart[1] + ii * delta[1],
												uvwStart[2] + ii * delta[2] };
			double newCoords[3];
			getPhysCoordsFromParamCoords(uvw, newCoords);
			EV.verts[ii] = m_pMesh->addVert(newCoords);
		}
		vertsOnEdges.insert(std::make_pair(E, EV));
	}
	else {
		iterEdges->second.m_totalDihed += dihedral;
		EV = iterEdges->second;
		if (EV.m_totalDihed > (2 - 1.e-8) * M_PI) {
			vertsOnEdges.erase(iterEdges);
		}
	}
}


typename exa_set<TriFaceVerts>::iterator CellDivider::getTriVerts(
		exa_set<TriFaceVerts> &vertsOnTris,
		const int face,
		bool& shouldErase) {
	int ind0 = faceVertIndices[face][0];
	int ind1 = faceVertIndices[face][1];
	int ind2 = faceVertIndices[face][2];

	emInt vert0 = cellVerts[ind0];
	emInt vert1 = cellVerts[ind1];
	emInt vert2 = cellVerts[ind2];
	TriFaceVerts TFVTemp(vert0, vert1, vert2);
	auto iterTris = vertsOnTris.find(TFVTemp);
	TFVTemp.freeVertMemory();
	if (iterTris == vertsOnTris.end()) {
		const double inv_nDivs = 1. / (nDivs);
		const double uvw0[] = { uvwIJK[ind0][0], uvwIJK[ind0][1], uvwIJK[ind0][2] };
		const double uvw1[] = { uvwIJK[ind1][0], uvwIJK[ind1][1], uvwIJK[ind1][2] };
		const double uvw2[] = { uvwIJK[ind2][0], uvwIJK[ind2][1], uvwIJK[ind2][2] };

		double deltaUVWInI[] = { (uvw1[0] - uvw0[0]) * inv_nDivs,
															(uvw1[1] - uvw0[1]) * inv_nDivs, (uvw1[2]
																	- uvw0[2])
																																* inv_nDivs };
		double deltaUVWInJ[] = { (uvw2[0] - uvw0[0]) * inv_nDivs,
															(uvw2[1] - uvw0[1]) * inv_nDivs, (uvw2[2]
																	- uvw0[2])
																																* inv_nDivs };
		TriFaceVerts TFV(vert0, vert1, vert2);
		TFV.allocVertMemory();

		for (int jj = 0; jj < nDivs - 2; jj++) {
			for (int ii = 0; ii < nDivs - 2 - jj; ii++) {
				double uvw[] = { uvw0[0] + deltaUVWInI[0] * (ii + 1)
													+ deltaUVWInJ[0] * (jj + 1),
													uvw0[1] + deltaUVWInI[1] * (ii + 1)
													+ deltaUVWInJ[1] * (jj + 1),
													uvw0[2] + deltaUVWInI[2] * (ii + 1)
													+ deltaUVWInJ[2] * (jj + 1) };
				double newCoords[3];
				getPhysCoordsFromParamCoords(uvw, newCoords);
				emInt vNew = m_pMesh->addVert(newCoords);
				TFV.intVerts[ii][jj] = vNew;
			}
		} // Done looping over all interior verts for the triangle.
		iterTris = vertsOnTris.insert(TFV).first;
		shouldErase = false;
	}
	else {
		shouldErase = true;
//		TFV = *iterTris;
//		vertsOnTris.erase(iterTris); // Will never need this again.
	}
	return iterTris;
}

void CellDivider::getQuadVerts(exa_set<QuadFaceVerts> &vertsOnQuads,
		const int face, QuadFaceVerts &QFV) {
	int ind0 = faceVertIndices[face][0];
	int ind1 = faceVertIndices[face][1];
	int ind2 = faceVertIndices[face][2];
	int ind3 = faceVertIndices[face][3];

	emInt vert0 = cellVerts[ind0];
	emInt vert1 = cellVerts[ind1];
	emInt vert2 = cellVerts[ind2];
	emInt vert3 = cellVerts[ind3];

	QuadFaceVerts QFVTemp(vert0, vert1, vert2, vert3);

	auto iterQuads = vertsOnQuads.find(QFVTemp);
	if (iterQuads == vertsOnQuads.end()) {
		const double inv_nDivs = 1. / (nDivs);
		const double uvw0[] = { uvwIJK[ind0][0], uvwIJK[ind0][1], uvwIJK[ind0][2] };
		const double uvw1[] = { uvwIJK[ind1][0], uvwIJK[ind1][1], uvwIJK[ind1][2] };
		const double uvw2[] = { uvwIJK[ind2][0], uvwIJK[ind2][1], uvwIJK[ind2][2] };
		const double uvw3[] = { uvwIJK[ind3][0], uvwIJK[ind3][1], uvwIJK[ind3][2] };

		double deltaInI[] = { (uvw1[0] - uvw0[0]) * inv_nDivs, (uvw1[1] - uvw0[1])
				* inv_nDivs,
													(uvw1[2] - uvw0[2]) * inv_nDivs };
		double deltaInJ[] = { (uvw3[0] - uvw0[0]) * inv_nDivs, (uvw3[1] - uvw0[1])
				* inv_nDivs,
													(uvw3[2] - uvw0[2]) * inv_nDivs };

		double crossDelta[] = { (uvw2[0] + uvw0[0] - uvw1[0] - uvw3[0])
				* (inv_nDivs * inv_nDivs),
														(uvw2[1] + uvw0[1] - uvw1[1] - uvw3[1]) * (inv_nDivs
																* inv_nDivs),
														(uvw2[2] + uvw0[2] - uvw1[2] - uvw3[2]) * (inv_nDivs
																* inv_nDivs) };
		QFV.corners[0] = vert0;
		QFV.corners[1] = vert1;
		QFV.corners[2] = vert2;
		QFV.corners[3] = vert3;
		QFV.setupSorted();


		for (int jj = 1; jj <= nDivs - 1; jj++) {
			for (int ii = 1; ii <= nDivs - 1; ii++) {
				double newCoords[3];
				double uvw[] = { uvw0[0] + deltaInI[0] * ii + deltaInJ[0] * jj
																+ crossDelta[0] * ii * jj,
																uvw0[1] + deltaInI[1] * ii + deltaInJ[1] * jj
																+ crossDelta[1] * ii * jj,
																uvw0[2] + deltaInI[2] * ii + deltaInJ[2] * jj
																+ crossDelta[2] * ii * jj };
				getPhysCoordsFromParamCoords(uvw, newCoords);
				emInt vNew = m_pMesh->addVert(newCoords);
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

void CellDivider::divideEdges(exa_map<Edge, EdgeVerts> &vertsOnEdges) {
	// Divide all the edges, including storing info about which new verts
	// are on which edges
	for (int iE = 0; iE < numEdges; iE++) {

		EdgeVerts EV;
		double dihedral = 0;
		getEdgeVerts(vertsOnEdges, iE, dihedral, EV);

		// Now transcribe these into the master table for this cell.
		emInt startIndex = 1000, endIndex = 1000;
		if (EV.verts[0] == cellVerts[edgeVertIndices[iE][0]]) {
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

void CellDivider::divideFaces(exa_set<TriFaceVerts> &vertsOnTris,
exa_set<QuadFaceVerts> &vertsOnQuads) {
	// Divide all the faces, including storing info about which new verts
	// are on which faces

	// The quad faces are first.
	for (int iF = 0; iF < numQuadFaces; iF++) {
		QuadFaceVerts QFV;
		getQuadVerts(vertsOnQuads, iF, QFV);
		// Now extract info from the QFV and stuff it into the cell's point
		// array.

		// Critical first step: identify which vert is which.
		emInt corner[] = { 1000, 1000, 1000, 1000 };

		for (int iC = 0; iC < 4; iC++) {
			const emInt corn = QFV.corners[iC];
			for (int iV = 0; iV < numVerts; iV++) {
				const emInt cand = cellVerts[iV];
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
		bool shouldErase = false;
		auto iterTris = getTriVerts(vertsOnTris, iF, shouldErase);
		// Now extract info from the TFV and stuff it into the Prismamid's point
		// array.

		// 1000 is way more points than cells have.
		emInt corner[] = { 1000, 1000, 1000 };
		// Critical first step: identify which vert is which.
		for (int iC = 0; iC < 3; iC++) {
			const emInt corn = iterTris->corners[iC];
//			const emInt corn = TFV.corners[iC];

			for (int iV = 0; iV < numVerts; iV++) {
				const emInt cand = cellVerts[iV];
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

				localVerts[II][JJ][KK] = iterTris->intVerts[ii][jj];
//				localVerts[II][JJ][KK] = TFV.intVerts[ii][jj];
			}
		}
		if (shouldErase) {
			iterTris->freeVertMemory();
			vertsOnTris.erase(iterTris);
		}
	}
}

