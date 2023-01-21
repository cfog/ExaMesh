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

#include <math.h>

#include "ExaMesh.h"
#include "GeomUtils.h"
#include "CellDivider.h"
#include <stdio.h>
void sortVerts3(const emInt input[3], emInt output[3]) {
	// This is insertion sort, specialized for three inputs.
	if (input[1] < input[0]) {
		output[0] = input[1];
		output[1] = input[0];
	} else {
		output[0] = input[0];
		output[1] = input[1];
	}
	if (input[2] < output[1]) {
		output[2] = output[1];
		if (input[2] < output[0]) {
			output[1] = output[0];
			output[0] = input[2];
		} else {
			output[1] = input[2];
		}
	} else {
		output[2] = input[2];
	}
}

TriFaceVerts::TriFaceVerts(const int nDivs, const emInt v0, const emInt v1,
		const emInt v2, const emInt type, const emInt elemInd) :
		FaceVerts(nDivs, 3) {
	m_volElem = elemInd;
	m_volElemType = type;
	setCorners(v0, v1, v2);
}

void TriFaceVerts::setupSorted() {
	sortVerts3(m_corners, m_sorted);
}

void TriFaceVerts::getVertAndST(const int ii, const int jj, emInt& vert,
		double st[2], const int rotCase) const {
	assert(isValidIJ(ii, jj));
	int trueI = -1, trueJ = -1;
	switch (rotCase) {
	case 1:
		// Just reversed orientation, no rotation
		trueI = ii;
		trueJ = jj;
		break;
	case 2:
		trueI = jj;
		trueJ = m_nDivs - ii - jj;
		break;
	case 3:
		trueI = m_nDivs - jj - ii;
		trueJ = ii;
		break;
	case -1:
		// Just reversed orientation, no rotation
		trueI = jj;
		trueJ = ii;
		break;
	case -2:
		trueI = m_nDivs - ii - jj;
		trueJ = jj;
		break;
	case -3:
		trueI = ii;
		trueJ = m_nDivs - jj - ii;
		break;
	default:
		assert(0);
		break;
	}
	double sTmp = m_param_st[trueI][trueJ][0];
	double tTmp = m_param_st[trueI][trueJ][1];
	switch (rotCase) {
	case 1:
		st[0] = sTmp;
		st[1] = tTmp;
		break;
	case 2:
		st[0] = tTmp;
		st[1] = 1 - sTmp - tTmp;
		break;
	case 3:
		st[0] = 1 - sTmp - tTmp;
		st[1] = sTmp;
		break;
	case -1:
		st[0] = tTmp;
		st[1] = sTmp;
		break;
	case -2:
		st[0] = 1 - sTmp - tTmp;
		st[1] = tTmp;
		break;
	case -3:
		st[0] = sTmp;
		st[1] = 1 - tTmp - sTmp;
		break;
	default:
		assert(0);
		break;
	}
	assert(isValidParam(st[0]));
	assert(isValidParam(st[1]));
	vert = getIntVertInd(trueI, trueJ);
}

bool operator==(const TriFaceVerts &a, const TriFaceVerts &b) {
	return (a.m_sorted[0] == b.m_sorted[0] && a.m_sorted[1] == b.m_sorted[1]
			&& a.m_sorted[2] == b.m_sorted[2]);
}

bool operator<(const TriFaceVerts &a, const TriFaceVerts &b) {
	return (a.m_sorted[0] < b.m_sorted[0]
			|| (a.m_sorted[0] == b.m_sorted[0] && a.m_sorted[1] < b.m_sorted[1])
			|| (a.m_sorted[0] == b.m_sorted[0] && a.m_sorted[1] == b.m_sorted[1]
					&& a.m_sorted[2] < b.m_sorted[2]));
}

void TriFaceVerts::computeParaCoords(const int ii, const int jj,
		double st[2]) const {
	assert(ii >= 1 && ii <= m_nDivs - jj - 1);
	assert(jj >= 1 && jj <= m_nDivs - 1);

	// This is sufficient only for the very simplest of tests.
	int iLeft = 0;
	int jLeft = jj;

	int iRight = m_nDivs - jj;
	int jRight = jj;

	int iBot = ii;
	int jBot = 0;

	int iTop = ii;
	int jTop = m_nDivs - ii;

	double stLeft[] =
			{ m_param_st[iLeft][jLeft][0], m_param_st[iLeft][jLeft][1] };
	double stRight[] = { m_param_st[iRight][jRight][0],
			m_param_st[iRight][jRight][1] };
	double stBot[] = { m_param_st[iBot][jBot][0], m_param_st[iBot][jBot][1] };
	double stTop[] = { m_param_st[iTop][jTop][0], m_param_st[iTop][jTop][1] };
	getFaceParametricIntersectionPoint(stLeft, stRight, stBot, stTop, st);
	assert(isValidParam(st[0]));
	assert(isValidParam(st[1]));
}

QuadFaceVerts::QuadFaceVerts(const int nDivs, const emInt v0, const emInt v1,
		const emInt v2, const emInt v3, const emInt type, const emInt elemInd) :
		FaceVerts(nDivs, 4) {
	m_volElem = elemInd;
	m_volElemType = type;
	setCorners(v0, v1, v2, v3);
}

void QuadFaceVerts::setupSorted() {
	sortVerts4(m_corners, m_sorted);
}

void sortVerts4(const emInt input[4], emInt output[4]) {
	// This is insertion sort, specialized for four inputs.
	if (input[1] < input[0]) {
		output[0] = input[1];
		output[1] = input[0];
	} else {
		output[0] = input[0];
		output[1] = input[1];
	}

	if (input[2] < output[1]) {
		output[2] = output[1];
		if (input[2] < output[0]) {
			output[1] = output[0];
			output[0] = input[2];
		} else {
			output[1] = input[2];
		}
	} else {
		output[2] = input[2];
	}

	if (input[3] < output[2]) {
		output[3] = output[2];
		if (input[3] < output[1]) {
			output[2] = output[1];
			if (input[3] < output[0]) {
				output[1] = output[0];
				output[0] = input[3];
			} else {
				output[1] = input[3];
			}
		} else {
			output[2] = input[3];
		}
	} else {
		output[3] = input[3];
	}
}

bool operator==(const QuadFaceVerts &a, const QuadFaceVerts &b) {
	return (a.m_sorted[0] == b.m_sorted[0] && a.m_sorted[1] == b.m_sorted[1]
			&& a.m_sorted[2] == b.m_sorted[2] && a.m_sorted[3] == b.m_sorted[3]);
}

bool operator<(const QuadFaceVerts &a, const QuadFaceVerts &b) {
	return (a.m_sorted[0] < b.m_sorted[0]
			|| (a.m_sorted[0] == b.m_sorted[0] && a.m_sorted[1] < b.m_sorted[1])
			|| (a.m_sorted[0] == b.m_sorted[0] && a.m_sorted[1] == b.m_sorted[1]
					&& a.m_sorted[2] < b.m_sorted[2])
			|| (a.m_sorted[0] == b.m_sorted[0] && a.m_sorted[1] == b.m_sorted[1]
					&& a.m_sorted[2] == b.m_sorted[2]
					&& a.m_sorted[3] < b.m_sorted[3]));
}

void QuadFaceVerts::computeParaCoords(const int ii, const int jj,
		double st[2]) const {
	assert(ii >= 1 && ii <= m_nDivs - 1);
	assert(jj >= 1 && jj <= m_nDivs - 1);

	int iLeft = 0;
	int jLeft = jj;

	int iRight = m_nDivs;
	int jRight = jj;

	int iBot = ii;
	int jBot = 0;

	int iTop = ii;
	int jTop = m_nDivs;

	double stLeft[] =
			{ m_param_st[iLeft][jLeft][0], m_param_st[iLeft][jLeft][1] };
	double stRight[] = { m_param_st[iRight][jRight][0],
			m_param_st[iRight][jRight][1] };
	double stBot[] = { m_param_st[iBot][jBot][0], m_param_st[iBot][jBot][1] };
	double stTop[] = { m_param_st[iTop][jTop][0], m_param_st[iTop][jTop][1] };
	getFaceParametricIntersectionPoint(stLeft, stRight, stBot, stTop, st);
	assert(isValidParam(st[0]));
	assert(isValidParam(st[1]));
}

void QuadFaceVerts::getVertAndST(const int ii, const int jj, emInt& vert,
		double st[2], const int rotCase) const {
	assert(isValidIJ(ii, jj));
	int trueI = -1, trueJ = -1;
	switch (rotCase) {
	case 1:
		// No change from the original
		trueI = ii;
		trueJ = jj;
		break;
	case 2:
		trueI = jj;
		trueJ = m_nDivs - ii;
		break;
	case 3:
		trueI = m_nDivs - ii;
		trueJ = m_nDivs - jj;
		break;
	case 4:
		trueI = m_nDivs - jj;
		trueJ = ii;
		break;
	case -1:
		// Just reversed orientation, no rotation
		trueI = jj;
		trueJ = ii;
		break;
	case -2:
		trueI = m_nDivs - ii;
		trueJ = jj;
		break;
	case -3:
		trueI = m_nDivs - jj;
		trueJ = m_nDivs - ii;
		break;
	case -4:
		trueI = ii;
		trueJ = m_nDivs - jj;
		break;
	default:
		assert(0);
		break;
	}
	double sTmp = m_param_st[trueI][trueJ][0];
	double tTmp = m_param_st[trueI][trueJ][1];
	switch (rotCase) {
	case 1:
		st[0] = sTmp;
		st[1] = tTmp;
		break;
	case 2:
		st[0] = tTmp;
		st[1] = 1 - sTmp;
		break;
	case 3:
		st[0] = 1 - sTmp;
		st[1] = 1 - tTmp;
		break;
	case 4:
		st[0] = 1 - tTmp;
		st[1] = sTmp;
		break;
	case -1:
		st[0] = tTmp;
		st[1] = sTmp;
		break;
	case -2:
		st[0] = 1 - sTmp;
		st[1] = tTmp;
		break;
	case -3:
		st[0] = 1 - tTmp;
		st[1] = 1 - sTmp;
		break;
	case -4:
		st[0] = sTmp;
		st[1] = 1 - tTmp;
		break;
	default:
		assert(0);
		break;
	}
	assert(isValidParam(st[0]));
	assert(isValidParam(st[1]));
	vert = getIntVertInd(trueI, trueJ);
}

int CellDivider::checkOrient3D(const emInt verts[4]) const {
	double coords0[3], coords1[3], coords2[3], coords3[3];
	m_pMesh->getCoords(verts[0], coords0);
	m_pMesh->getCoords(verts[1], coords1);
	m_pMesh->getCoords(verts[2], coords2);
	m_pMesh->getCoords(verts[3], coords3);
	return ::checkOrient3D(coords0, coords1, coords2, coords3);
}

void CellDivider::getEdgeParametricDivision(EdgeVerts &EV) const {
	emInt startInd = EV.m_verts[0];
	emInt endInd = EV.m_verts[nDivs];
	double startLenOrig = m_Map->getIsoLengthScale(startInd);
	double endLenOrig = m_Map->getIsoLengthScale(endInd);

	// Working out parametric coordinate along the edge, so (see Carl's
	// notes for June 16, 2020, or the journal article) x_A = 0, dx = 1,
	// and:
	//
	// u = startLen * xi + (3 - 2*startLen - endLen) * xi^2
	//     + (startLen + endLen - 2) * xi^3
	assert(isfinite(startLenOrig) && startLenOrig > 0);
	assert(isfinite(endLenOrig) && endLenOrig > 0);

	double startLen = sqrt(startLenOrig / endLenOrig);
	double endLen = 1. / startLen;
	assert(isfinite(startLen));
	assert(isfinite(endLen));

	for (int ii = 1; ii < nDivs; ii++) {
		double xi = ((double) ii) / nDivs;
		EV.m_param_t[ii] = startLen * xi + (3 - 2 * startLen - endLen) * xi * xi
				+ (startLen + endLen - 2) * xi * xi * xi;
	}
	EV.m_param_t[0] = 0;
	EV.m_param_t[nDivs] = 1;
}

void CellDivider::getEdgeVerts(exa_map<Edge, EdgeVerts> &vertsOnEdges,
		const int edge, const double dihedral, EdgeVerts &EV) {
	int ind0 = edgeVertIndices[edge][0];
	int ind1 = edgeVertIndices[edge][1];

	emInt vert0 = cellVerts[ind0];
	emInt vert1 = cellVerts[ind1];
//	printf("Edge: %5d %5d ", vert0, vert1);

	Edge E(vert0, vert1);
	auto iterEdges = vertsOnEdges.find(E);

	if (iterEdges == vertsOnEdges.end()) {
//		printf("new\n");
		// Doesn't exist yet, so create it.
		EV.m_verts[0] = E.getV0();
		EV.m_verts[nDivs] = E.getV1();
		EV.m_totalDihed = dihedral;

		bool forward = true;
		if (EV.m_verts[0] != vert0) {
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
		} else {
			uvwStart[0] = uvwIJK[ind1][0];
			uvwStart[1] = uvwIJK[ind1][1];
			uvwStart[2] = uvwIJK[ind1][2];

			uvwEnd[0] = uvwIJK[ind0][0];
			uvwEnd[1] = uvwIJK[ind0][1];
			uvwEnd[2] = uvwIJK[ind0][2];
		}
		double delta[] = { (uvwEnd[0] - uvwStart[0]), (uvwEnd[1] - uvwStart[1]),
				(uvwEnd[2] - uvwStart[2]) };
		getEdgeParametricDivision(EV);
		for (int ii = 1; ii < nDivs; ii++) {
			double uvw[] = { uvwStart[0] + EV.m_param_t[ii] * delta[0],
					uvwStart[1] + EV.m_param_t[ii] * delta[1], uvwStart[2]
							+ EV.m_param_t[ii] * delta[2] };
			double newCoords[3];
			getPhysCoordsFromParamCoords(uvw, newCoords);
			EV.m_verts[ii] = m_pMesh->addVert(newCoords);
//			printf("%3d %5f (%5f %5f %5f) (%8f %8f %8f)\n",
//					ii, EV.m_param_t[ii], uvw[0], uvw[1], uvw[2],
//					newCoords[0], newCoords[1], newCoords[2]);
		}
		vertsOnEdges.insert(std::make_pair(E, EV));
	} else {
//		printf("old\n");
		iterEdges->second.m_totalDihed += dihedral;
		EV = iterEdges->second;
		if (EV.m_totalDihed > (2 - 1.e-8) * M_PI) {
			vertsOnEdges.erase(iterEdges);
		}
	}
}

void getFaceParametricIntersectionPoint(const double uvL[2],
		const double uvR[2], const double uvB[2], const double uvT[2],
		double uv[2]) {
	double delta_uvTB[] = { uvT[0] - uvB[0], uvT[1] - uvB[1] };
	double delta_uvRL[] = { uvR[0] - uvL[0], uvR[1] - uvL[1] };

	double &duRL = delta_uvRL[0];
	double &dvRL = delta_uvRL[1];

	double &duTB = delta_uvTB[0];
	double &dvTB = delta_uvTB[1];

	double denom = duRL * dvTB - duTB * dvRL;

	uv[0] = (-duRL * duTB * uvB[1] + duRL * duTB * uvL[1] + duRL * dvTB * uvB[0]
			- duTB * dvRL * uvL[0]) / denom;
	uv[1] = (duRL * dvTB * uvL[1] - duTB * dvRL * uvB[1] + dvRL * dvTB * uvB[0]
			- dvRL * dvTB * uvL[0]) / denom;
}

bool CellDivider::isEdgeForwardForFace(const EdgeVerts &EV, emInt cornerStart,
		emInt cornerEnd) const {
	assert(EV.m_verts[0] == cornerEnd || EV.m_verts[0] == cornerStart);
	assert(EV.m_verts[nDivs] == cornerEnd || EV.m_verts[nDivs] == cornerStart);
	assert(cornerEnd != cornerStart);
	bool isForward = true;
	if (EV.m_verts[0] == cornerEnd && EV.m_verts[nDivs] == cornerStart)
		isForward = false;
	return isForward;
}

void CellDivider::initPerimeterParams(TriFaceVerts &TFV, const int face) const {
	// Need to identify which edge of the tri is which previously defined
	// edge, and use the edge parameter info to set up parameter info
	// for the triangle.

	// First edge

	// Grab the appropriate edge.  Indices are shifted by numEdges for
	// edges that are reversed in the face compared to the edge definition.
	{
		int actualEdge = faceEdgeIndices[face][0];
		assert(actualEdge >= 0 && actualEdge < numEdges);
		const EdgeVerts &EV = m_EV[actualEdge];
		emInt cornerStart = TFV.getCorner(0);
		emInt cornerEnd = TFV.getCorner(1);
		bool isForward = isEdgeForwardForFace(EV, cornerStart, cornerEnd);

		// Transcribe that edge's parametric division on the triangle.  This
		// is idiosyncratic enough for edges that I'm not trying to loop it.
		for (int pp = 0; pp <= nDivs; pp++) {
			double st[] = { EV.m_param_t[pp], 0 };
			int ii = pp, jj = 0;
			if (!isForward) {
				st[0] = 1 - st[0];
				ii = nDivs - ii;
			}
			TFV.setVertSTParams(ii, jj, st);
		}
	}
	// Second edge

	// Grab the appropriate edge.  Indices are shifted by numEdges for
	// edges that are reversed in the face compared to the edge definition.
	{
		int actualEdge = faceEdgeIndices[face][1];
		assert(actualEdge >= 0 && actualEdge < numEdges);
		const EdgeVerts &EV = m_EV[actualEdge];
		emInt cornerStart = TFV.getCorner(1);
		emInt cornerEnd = TFV.getCorner(2);
		bool isForward = isEdgeForwardForFace(EV, cornerStart, cornerEnd);

		// Transcribe that edge's parametric division on the triangle.  This
		// is idiosyncratic enough for edges that I'm not trying to loop it.
		for (int pp = 0; pp <= nDivs; pp++) {
			// This is the hypotenuse in the st parametric space, so
			// the edge runs from (s,t) = (1,0) to (0,1).
			double st[] = { 0, 0 };
			// Get t right.
			st[1] = EV.m_param_t[pp];
			int jj = pp;
			if (!isForward) {
				st[1] = 1 - st[1];
				jj = nDivs - jj;
			}
			// Now take advantage of knowing the s+t = 1
			st[0] = 1 - st[1];
			TFV.setVertSTParams(nDivs - jj, jj, st);
		}
	}

	// Third edge

	// Grab the appropriate edge.  Indices are shifted by numEdges for
	// edges that are reversed in the face compared to the edge definition.
	{
		int actualEdge = faceEdgeIndices[face][2];
		assert(actualEdge >= 0 && actualEdge < numEdges);
		const EdgeVerts &EV = m_EV[actualEdge];
		emInt cornerStart = TFV.getCorner(2);
		emInt cornerEnd = TFV.getCorner(0);
		bool isForward = isEdgeForwardForFace(EV, cornerStart, cornerEnd);

		// Transcribe that edge's parametric division on the triangle.  This
		// is idiosyncratic enough for edges that I'm not trying to loop it.
		for (int pp = 0; pp <= nDivs; pp++) {
			// This edge runs from (s,t) = (0,1) to (0,0)
			double st[] = { 0, 1 - EV.m_param_t[pp] };
			int jj = nDivs - pp;
			if (!isForward) {
				st[1] = 1 - st[1];
				jj = nDivs - jj;
			}
			TFV.setVertSTParams(0, jj, st);
		}
	}
}

TriFaceVerts CellDivider::getTriVerts(
exa_set<TriFaceVerts> &vertsOnTris, const int face) {
	int ind0 = faceVertIndices[face][0];
	int ind1 = faceVertIndices[face][1];
	int ind2 = faceVertIndices[face][2];

	emInt vert0 = cellVerts[ind0];
	emInt vert1 = cellVerts[ind1];
	emInt vert2 = cellVerts[ind2];

	const double uvw0[] = { uvwIJK[ind0][0], uvwIJK[ind0][1], uvwIJK[ind0][2] };
	const double uvw1[] = { uvwIJK[ind1][0], uvwIJK[ind1][1], uvwIJK[ind1][2] };
	const double uvw2[] = { uvwIJK[ind2][0], uvwIJK[ind2][1], uvwIJK[ind2][2] };

	double deltaUVWInS[] = { (uvw1[0] - uvw0[0]), (uvw1[1] - uvw0[1]), (uvw1[2]
			- uvw0[2]) };
	double deltaUVWInT[] = { (uvw2[0] - uvw0[0]), (uvw2[1] - uvw0[1]), (uvw2[2]
			- uvw0[2]) };

	TriFaceVerts TFV(nDivs, vert0, vert1, vert2);
	initPerimeterParams(TFV, face);

	// Find the existing iterator if the face has already been operated
	// on once.
	auto iterTris = vertsOnTris.find(TFV);
	bool newFace = (iterTris == vertsOnTris.end());
	int rotCase = 0;
	if (!newFace) {
		for (int cc = 0; cc < 3; cc++) {
			if (vert0 == iterTris->getCorner(cc)) {
				if (vert1 == iterTris->getCorner((cc+1)%3)) {
					// Oriented forward; bdry tri
					assert(vert2 == iterTris->getCorner((cc+2)%3));
					rotCase = cc+1;
				}
				else {
					assert(vert1 == iterTris->getCorner((cc+2)%3));
					assert(vert2 == iterTris->getCorner((cc+1)%3));
					rotCase = -(cc+1);
				}
			}
		}
		assert(rotCase != 0);
	}

	for (int jj = 1; jj <= nDivs - 2; jj++) {
		for (int ii = 1; ii <= nDivs - 1 - jj; ii++) {
			double st[] = { -100, -100 };
			emInt vert;
			if (newFace) {
				TFV.computeParaCoords(ii, jj, st);
				TFV.setVertSTParams(ii, jj, st);
			} else {
				iterTris->getVertAndST(ii, jj, vert, st, rotCase);
			}
			double &s = st[0];
			double &t = st[1];
			assert(s >= 0 && t >= 0 && (s + t) <= 1);
			double uvw[] = { uvw0[0] + deltaUVWInS[0] * s + deltaUVWInT[0] * t,
					uvw0[1] + deltaUVWInS[1] * s + deltaUVWInT[1] * t, uvw0[2]
							+ deltaUVWInS[2] * s + deltaUVWInT[2] * t };
			TFV.setVertUVWParams(ii, jj, uvw);
			if (newFace) {
				double newCoords[3];
				getPhysCoordsFromParamCoords(uvw, newCoords);
				vert = m_pMesh->addVert(newCoords);
			}
			TFV.setIntVertInd(ii, jj, vert);
		}
	} // Done looping over all interior verts for the triangle.
	if (newFace) {
		vertsOnTris.insert(TFV);
	} else {
		vertsOnTris.erase(iterTris); // Will never need this again.
	}
	return TFV;
}

void CellDivider::initPerimeterParams(QuadFaceVerts &QFV,
		const int face) const {
// Need to identify which edge of the quad is which previously defined
// edge, and use the edge parameter info to set up parameter info
// for the quad.

// First edge

// Grab the appropriate edge.
	{
		int actualEdge = faceEdgeIndices[face][0];
		assert(actualEdge >= 0 && actualEdge < numEdges);
		const EdgeVerts &EV = m_EV[actualEdge];
		emInt cornerStart = QFV.getCorner(0);
		emInt cornerEnd = QFV.getCorner(1);
		bool isForward = isEdgeForwardForFace(EV, cornerStart, cornerEnd);

		// Transcribe that edge's parametric division on the quad.  This
		// is idiosyncratic enough for edges that I'm not trying to loop it.
		for (int pp = 0; pp <= nDivs; pp++) {
			double st[] = { EV.m_param_t[pp], 0 };
			int ii = pp, jj = 0;
			if (!isForward) {
				st[0] = 1 - st[0];
				ii = nDivs - ii;
			}
			QFV.setVertSTParams(ii, jj, st);
		}
	}
// Second edge

// Grab the appropriate edge.
	{
		int actualEdge = faceEdgeIndices[face][1];
		assert(actualEdge >= 0 && actualEdge < numEdges);
		const EdgeVerts &EV = m_EV[actualEdge];
		emInt cornerStart = QFV.getCorner(1);
		emInt cornerEnd = QFV.getCorner(2);
		bool isForward = isEdgeForwardForFace(EV, cornerStart, cornerEnd);

		// Transcribe that edge's parametric division on the quad.  This
		// is idiosyncratic enough for edges that I'm not trying to loop it.
		for (int pp = 0; pp <= nDivs; pp++) {
			// This is the hypotenuse in the st parametric space, so
			// the edge runs from (s,t) = (1,0) to (0,1).
			double st[] = { 1, EV.m_param_t[pp] };
			int jj = pp;
			if (!isForward) {
				st[1] = 1 - st[1];
				jj = nDivs - pp;
			}
			QFV.setVertSTParams(nDivs, jj, st);
		}
	}

// Third edge

// Grab the appropriate edge.
	{
		int actualEdge = faceEdgeIndices[face][2];
		assert(actualEdge >= 0 && actualEdge < numEdges);
		const EdgeVerts &EV = m_EV[actualEdge];
		emInt cornerStart = QFV.getCorner(2);
		emInt cornerEnd = QFV.getCorner(3);
		bool isForward = isEdgeForwardForFace(EV, cornerStart, cornerEnd);

		// Transcribe that edge's parametric division on the quad.  This
		// is idiosyncratic enough for edges that I'm not trying to loop it.
		for (int pp = 0; pp <= nDivs; pp++) {
			// This edge runs from (s,t) = (1,1) to (0,1)
			double st[] = { 1 - EV.m_param_t[pp], 1 };
			int ii = nDivs - pp;
			if (!isForward) {
				st[0] = 1 - st[0];
				ii = pp;
			}
			QFV.setVertSTParams(ii, nDivs, st);
		}
	}

// Fourth edge

// Grab the appropriate edge.
	{
		int actualEdge = faceEdgeIndices[face][3];
		assert(actualEdge >= 0 && actualEdge < numEdges);
		const EdgeVerts &EV = m_EV[actualEdge];
		emInt cornerStart = QFV.getCorner(3);
		emInt cornerEnd = QFV.getCorner(0);
		bool isForward = isEdgeForwardForFace(EV, cornerStart, cornerEnd);

		// Transcribe that edge's parametric division on the quad.  This
		// is idiosyncratic enough for edges that I'm not trying to loop it.
		for (int pp = 0; pp <= nDivs; pp++) {
			// This edge runs from (s,t) = (0,1) to (0,0)
			double st[] = { 0, 1 - EV.m_param_t[pp] };
			int jj = nDivs - pp;
			if (!isForward) {
				st[1] = 1 - st[1];
				jj = nDivs - jj;
			}
			QFV.setVertSTParams(0, jj, st);
		}
	}
}

QuadFaceVerts CellDivider::getQuadVerts(exa_set<QuadFaceVerts> &vertsOnQuads,
		const int face) {
	int ind0 = faceVertIndices[face][0];
	int ind1 = faceVertIndices[face][1];
	int ind2 = faceVertIndices[face][2];
	int ind3 = faceVertIndices[face][3];

	emInt vert0 = cellVerts[ind0];
	emInt vert1 = cellVerts[ind1];
	emInt vert2 = cellVerts[ind2];
	emInt vert3 = cellVerts[ind3];
//	printf("%6d %6d %6d %6d ", vert0, vert1, vert2, vert3);

	const double uvw0[] = { uvwIJK[ind0][0], uvwIJK[ind0][1],
			uvwIJK[ind0][2] };
	const double uvw1[] = { uvwIJK[ind1][0], uvwIJK[ind1][1],
			uvwIJK[ind1][2] };
	const double uvw2[] = { uvwIJK[ind2][0], uvwIJK[ind2][1],
			uvwIJK[ind2][2] };
	const double uvw3[] = { uvwIJK[ind3][0], uvwIJK[ind3][1],
			uvwIJK[ind3][2] };

	const double deltaInS[] = { (uvw1[0] - uvw0[0]), (uvw1[1] - uvw0[1]), (uvw1[2]
			- uvw0[2]) };
	const double deltaInT[] = { (uvw3[0] - uvw0[0]), (uvw3[1] - uvw0[1]), (uvw3[2]
			- uvw0[2]) };
	const double crossDelta[] = { (uvw2[0] + uvw0[0] - uvw1[0] - uvw3[0]),
			(uvw2[1] + uvw0[1] - uvw1[1] - uvw3[1]), (uvw2[2] + uvw0[2]
					- uvw1[2] - uvw3[2]) };

	QuadFaceVerts QFV(nDivs, vert0, vert1, vert2, vert3);
	initPerimeterParams(QFV, face);

	auto iterQuads = vertsOnQuads.find(QFV);

	bool newFace = (iterQuads == vertsOnQuads.end());
	int rotCase = 0;
	if (!newFace) {
		for (int cc = 0; cc < 4; cc++) {
			if (vert0 == iterQuads->getCorner(cc)) {
				if (vert1 == iterQuads->getCorner((cc+1)%4)) {
					// Oriented forward; bdry quad
					assert(vert2 == iterQuads->getCorner((cc+2)%4));
					assert(vert3 == iterQuads->getCorner((cc+3)%4));
					rotCase = cc+1;
				}
				else {
					assert(vert1 == iterQuads->getCorner((cc+3)%4));
					assert(vert2 == iterQuads->getCorner((cc+2)%4));
					assert(vert3 == iterQuads->getCorner((cc+1)%4));
					rotCase = -(cc+1);
				}
			}
		}
		assert(rotCase != 0);
	}

	for (int jj = 1; jj <= nDivs - 1; jj++) {
		for (int ii = 1; ii <= nDivs - 1; ii++) {
			double st[] = { -100, -100 };
			double &s = st[0], &t = st[1];
			emInt vert = -100;
			if (newFace) {
				QFV.computeParaCoords(ii, jj, st);
				QFV.setVertSTParams(ii, jj, st);
			} else {
				iterQuads->getVertAndST(ii, jj, vert, st, rotCase);
			}
			assert(s >= 0 && s <= 1 && t >= 0 && t <= 1);
			double uvw[] = { uvw0[0] + deltaInS[0] * s + deltaInT[0] * t
					+ crossDelta[0] * s * t,
					uvw0[1] + deltaInS[1] * s + deltaInT[1] * t
					+ crossDelta[1] * s * t,
					uvw0[2] + deltaInS[2] * s + deltaInT[2] * t
						+ crossDelta[2] * s * t };
			QFV.setVertUVWParams(ii, jj, uvw);
			if (newFace) {
				double newCoords[3];
				getPhysCoordsFromParamCoords(uvw, newCoords);
				vert = m_pMesh->addVert(newCoords);
			}
			QFV.setIntVertInd(ii, jj, vert);
		}
	} // Done looping over all interior verts for the quad.
	if (newFace) {
		vertsOnQuads.insert(QFV);
	} else {
		vertsOnQuads.erase(iterQuads); // Will never need this again.
	}
	return QFV;
}

void CellDivider::divideEdges(exa_map<Edge, EdgeVerts> &vertsOnEdges) {
// Divide all the edges, including storing info about which new verts
// are on which edges
	for (int iE = 0; iE < numEdges; iE++) {

		EdgeVerts &EV = m_EV[iE];
		double dihedral = 0;
		getEdgeVerts(vertsOnEdges, iE, dihedral, EV);

		// Now transcribe these into the master table for this cell.
		emInt startIndex = 1000, endIndex = 1000;
		if (EV.m_verts[0] == cellVerts[edgeVertIndices[iE][0]]) {
			// Transcribe this edge forward.
			startIndex = edgeVertIndices[iE][0];
			endIndex = edgeVertIndices[iE][1];
		} else {
			startIndex = edgeVertIndices[iE][1];
			endIndex = edgeVertIndices[iE][0];
		}
		int startI = vertIJK[startIndex][0];
		int startJ = vertIJK[startIndex][1];
		int startK = vertIJK[startIndex][2];
		int incrI = (vertIJK[endIndex][0] - startI) / nDivs;
		int incrJ = (vertIJK[endIndex][1] - startJ) / nDivs;
		int incrK = (vertIJK[endIndex][2] - startK) / nDivs;

		double startU = uvwIJK[startIndex][0];
		double startV = uvwIJK[startIndex][1];
		double startW = uvwIJK[startIndex][2];
		double deltaU = (uvwIJK[endIndex][0] - startU);
		double deltaV = (uvwIJK[endIndex][1] - startV);
		double deltaW = (uvwIJK[endIndex][2] - startW);

		for (int ii = 0; ii <= nDivs; ii++) {
			int II = startI + ii * incrI;
			int JJ = startJ + ii * incrJ;
			int KK = startK + ii * incrK;
			assert(II >= 0 && II <= nDivs);
			assert(JJ >= 0 && JJ <= nDivs);
			assert(KK >= 0 && KK <= nDivs);
			localVerts[II][JJ][KK] = EV.m_verts[ii];

			// Now the param coords
			double u = startU + EV.m_param_t[ii] * deltaU;
			double v = startV + EV.m_param_t[ii] * deltaV;
			double w = startW + EV.m_param_t[ii] * deltaW;
			m_uvw[II][JJ][KK][0] = u;
			m_uvw[II][JJ][KK][1] = v;
			m_uvw[II][JJ][KK][2] = w;
		}
	}
}

void CellDivider::divideFaces(exa_set<TriFaceVerts> &vertsOnTris,
exa_set<QuadFaceVerts> &vertsOnQuads) {
// Divide all the faces, including storing info about which new verts
// are on which faces

// The quad faces are first.
	for (int iF = 0; iF < numQuadFaces; iF++) {
		QuadFaceVerts QFV = getQuadVerts(vertsOnQuads, iF);
		// Now extract info from the QFV and stuff it into the cell's point
		// array.

		// Critical first step: identify which vert is which.
		emInt corner[] = { 1000, 1000, 1000, 1000 };

		for (int iC = 0; iC < 4; iC++) {
			const emInt corn = QFV.getCorner(iC);
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

				localVerts[II][JJ][KK] = QFV.getIntVertInd(ii, jj);
				double uvw[3];
				QFV.getVertUVWParams(ii, jj, uvw);
				m_uvw[II][JJ][KK][0] = uvw[0];
				m_uvw[II][JJ][KK][1] = uvw[1];
				m_uvw[II][JJ][KK][2] = uvw[2];
			}
		}
	}

	for (int iF = numQuadFaces; iF < numQuadFaces + numTriFaces; iF++) {
		TriFaceVerts TFV = getTriVerts(vertsOnTris, iF);
		// Now extract info from the TFV and stuff it into the cell's point
		// array.

		// 1000 is way more points than cells have.
		emInt corner[] = { 1000, 1000, 1000 };
		// Critical first step: identify which vert is which.
		for (int iC = 0; iC < 3; iC++) {
			const emInt corn = TFV.getCorner(iC);
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

		for (int jj = 1; jj <= nDivs - 2; jj++) {
			for (int ii = 1; ii <= nDivs - 1 - jj; ii++) {
				int II = startI + incrIi * ii + incrIj * jj;
				int JJ = startJ + incrJi * ii + incrJj * jj;
				int KK = startK + incrKi * ii + incrKj * jj;
				assert(II >= 0 && II <= nDivs);
				assert(JJ >= 0 && JJ <= nDivs);
				assert(KK >= 0 && KK <= nDivs);

				localVerts[II][JJ][KK] = TFV.getIntVertInd(ii, jj);
				double uvw[3];
				TFV.getVertUVWParams(ii, jj, uvw);
				m_uvw[II][JJ][KK][0] = uvw[0];
				m_uvw[II][JJ][KK][1] = uvw[1];
				m_uvw[II][JJ][KK][2] = uvw[2];
			}
		}
	}
}

void CellDivider::computeParaCoords(const int ii, const int jj, const int kk,
		double uvw[3]) const {
	int iMin = minI(jj, kk);
	int iMax = maxI(jj, kk);

	int jMin = minJ(ii, kk);
	int jMax = maxJ(ii, kk);

	int kMin = minK(ii, jj);
	int kMax = maxK(ii, jj);

	double *uvwL = m_uvw[iMin][jj][kk];
	double *uvwR = m_uvw[iMax][jj][kk];

	double *uvwBot = m_uvw[ii][jMin][kk];
	double *uvwTop = m_uvw[ii][jMax][kk];

	double *uvwFr = m_uvw[ii][jj][kMin];
	double *uvwBa = m_uvw[ii][jj][kMax];

	getCellInteriorParametricIntersectionPoint(uvwL, uvwR, uvwBot, uvwTop, uvwFr,
			uvwBa, uvw);
	assert(uvw[0] >= 0 && uvw[0] <= 1);
	assert(uvw[1] >= 0 && uvw[1] <= 1);
	assert(uvw[2] >= 0 && uvw[2] <= 1);
}

void CellDivider::divideInterior() {
// Number of verts added:
//    Tets:      (nD-1)(nD-2)(nD-3)/2
//    Pyrs:      (nD-1)(nD-2)(2 nD-3)/6
//    Prisms:    (nD-1)(nD-2)nD/2
//    Hexes:     (nD-1)^3
	if (nDivs < getMinInteriorDivs())
	return;
	for (int kk = 1; kk < nDivs; kk++) {
		int jMax = maxJ(1, kk);
		for (int jj = 1; jj < jMax; jj++) {
			int iMax = maxI(jj, kk);
			for (int ii = 1; ii < iMax; ii++) {
				double uvw[3];
				// Now find uvw by finding the near-intersection point
				// of the lines of constant i, j, k
				computeParaCoords(ii, jj, kk, uvw);

				double &u = uvw[0];
				double &v = uvw[1];
				double &w = uvw[2];
				double coords[3];
				m_Map->computeTransformedCoords(uvw, coords);
				emInt vNew = m_pMesh->addVert(coords);
				localVerts[ii][jj][kk] = vNew;
				m_uvw[ii][jj][kk][0] = u;
				m_uvw[ii][jj][kk][1] = v;
				m_uvw[ii][jj][kk][2] = w;
//				printf("%3d %3d %3d %5f %5f %5f\n", ii, jj, kk, u, v, w);
			}
		}
	} // Done looping to create all verts inside the cell.
}

void getCellInteriorParametricIntersectionPoint(const double uvwA[3],
		const double uvwB[3], const double uvwC[3], const double uvwD[3],
		const double uvwE[3], const double uvwF[3], double uvw[3]) {
// Points A and B are opposite each other.
// Points C and D are opposite each other.
// Points E and F are opposite each other.

	double deltaBA[] = {uvwB[0] - uvwA[0], uvwB[1] - uvwA[1], uvwB[2] - uvwA[2]};
	double deltaDC[] = {uvwD[0] - uvwC[0], uvwD[1] - uvwC[1], uvwD[2] - uvwC[2]};
	double deltaFE[] = {uvwF[0] - uvwE[0], uvwF[1] - uvwE[1], uvwF[2] - uvwE[2]};

// The three lines in parametric space connecting these three
// pairs of points almost certainly don't intersect.  So find
// the closest point (in a least-squares sense) to the three
// lines.

// Along line AB, we have u = uA + (uB-uA) s1, and so on for other
// variables.  Along CD, the parameter is s2, and along EF, it's s3.
// So we'd like to have
// uA + (uB-uA) s1 = uC + (uD-uC) s2 = uE + (uF-uE) s3

// First, set up the least-squares system, with all nine of the
// equalities that we'd like to satisfy.

	double L2LHS[9][3] = { {deltaBA[0], -deltaDC[0], 0}, {deltaBA[1],
			-deltaDC[1], 0}, {deltaBA[2], -deltaDC[2], 0}, {0, deltaDC[0],
			-deltaFE[0]}, {0, deltaDC[1], -deltaFE[1]}, {0, deltaDC[2],
			-deltaFE[2]}, {-deltaBA[0], 0, deltaFE[0]}, {-deltaBA[1], 0,
			deltaFE[1]}, {-deltaBA[2], 0, deltaFE[2]}};
	double L2RHS[9] = {uvwC[0] - uvwA[0], uvwC[1] - uvwA[1], uvwC[2] - uvwA[2],
		uvwE[0] - uvwC[0], uvwE[1] - uvwC[1], uvwE[2] - uvwC[2], uvwA[0]
		- uvwE[0], uvwA[1] - uvwE[1], uvwA[2] - uvwE[2]};

// Now multiply L2LHS . s = L2RHS from the left by L2LHS^T
	double LHS[3][3], RHS[3];
	for (int ii = 0; ii < 3; ii++) {
		RHS[ii] = 0;
		for (int kk = 0; kk < 9; kk++) {
			RHS[ii] += L2RHS[kk] * L2LHS[kk][ii];
		}
		for (int jj = 0; jj < 3; jj++) {
			LHS[ii][jj] = 0;
			for (int kk = 0; kk < 9; kk++) {
				LHS[ii][jj] += L2LHS[kk][ii] * L2LHS[kk][jj];
			}
		}
	}

// Now we've got a 3x3 system to solve for the s's.
// LHS is symmetric: [aa bb cc]
//                   [bb dd ee]
//                   [cc ee ff]

	double &aa = LHS[0][0];
	double &bb = LHS[1][0];
	double &cc = LHS[2][0];
	double &dd = LHS[1][1];
	double &ee = LHS[2][1];
	double &ff = LHS[2][2];
	double LHSdet = aa * dd * ff - aa * ee * ee - ff * bb * bb - dd * cc * cc
	+ 2 * bb * ee * cc;
	double solnS[] = {(RHS[0] * (dd * ff - ee * ee) + RHS[1] * (cc * ee - bb * ff)
				+ RHS[2] * (bb * ee - cc * dd)) / LHSdet, (RHS[0] * (cc * ee - bb * ff)
				+ RHS[1] * (aa * ff - cc * cc) + RHS[2] * (bb * cc - aa * ee)) / LHSdet,
		(RHS[0] * (bb * ee - cc * dd) + RHS[1] * (bb * cc - aa * ee)
				+ RHS[2] * (aa * dd - bb * bb)) / LHSdet};

// All of these curve params had better fall in the range (0,1)
	assert(solnS[0] >= 0 && solnS[0] <= 1);
	assert(solnS[1] >= 0 && solnS[1] <= 1);
	assert(solnS[2] >= 0 && solnS[2] <= 1);

// Unfortunately, the three uvw triples that these s's imply aren't the same.
// So let's compute all three of them and then average.

	double uvwAB[] = {uvwA[0] + deltaBA[0] * solnS[0], uvwA[1]
		+ deltaBA[1] * solnS[0], uvwA[2] + deltaBA[2] * solnS[0]};
	double uvwCD[] = {uvwC[0] + deltaDC[0] * solnS[1], uvwC[1]
		+ deltaDC[1] * solnS[1], uvwC[2] + deltaDC[2] * solnS[1]};
	double uvwEF[] = {uvwE[0] + deltaFE[0] * solnS[2], uvwE[1]
		+ deltaFE[1] * solnS[2], uvwE[2] + deltaFE[2] * solnS[2]};

	uvw[0] = (uvwAB[0] + uvwCD[0] + uvwEF[0]) / 3;
	uvw[1] = (uvwAB[1] + uvwCD[1] + uvwEF[1]) / 3;
	uvw[2] = (uvwAB[2] + uvwCD[2] + uvwEF[2]) / 3;
}

void CellDivider::printAllPoints() {
	for (int kk = 0; kk <= nDivs; kk++) {
		for (int jj = minJ(0, kk); jj <= maxJ(0, kk); jj++) {
			for (int ii = minI(jj, kk); ii <= maxI(jj, kk); ii++) {
				double *uvw = m_uvw[ii][jj][kk];
				emInt point = localVerts[ii][jj][kk];
				printf("%3d %3d %3d (%5f %5f %5f) (%8f %8f %8f)\n", ii, jj, kk,
						uvw[0], uvw[1], uvw[2], m_pMesh->getX(point),
						m_pMesh->getY(point), m_pMesh->getZ(point));
			}
			printf("\n");
		}
		printf("\n");
	}
}
