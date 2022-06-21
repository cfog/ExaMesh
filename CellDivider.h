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
 * CellDivider.h
 *
 *  Created on: Jul. 3, 2019
 *      Author: cfog
 */

#ifndef SRC_CELLDIVIDER_H_
#define SRC_CELLDIVIDER_H_

#include <algorithm>
#include <cmath>

#include "exa-defs.h"
#include "ExaMesh.h"
#include "Mapping.h"
#include "UMesh.h"

class CellDivider {
protected:
	UMesh *m_pMesh;
	Mapping *m_Map;
	emInt (*localVerts)[MAX_DIVS + 1][MAX_DIVS + 1];
	double (*m_uvw)[MAX_DIVS+1][MAX_DIVS+1][3];
	int edgeVertIndices[12][2];
	int faceVertIndices[6][4];
	int faceEdgeIndices[6][4];
	EdgeVerts m_EV[12];
	int numTriFaces, numQuadFaces, numEdges, numVerts;
	int vertIJK[8][3];
	double uvwIJK[8][3];
	emInt cellVerts[8];
	int nDivs;

	// Used by both tets and pyramids.
	int checkOrient3D(const emInt verts[4]) const;
private:
	void getEdgeVerts(exa_map<Edge, EdgeVerts> &vertsOnEdges, const int edge,
			const double dihedral, EdgeVerts &EV);

	void getQuadVerts(exa_set<QuadFaceVerts> &vertsOnQuads, const int face,
			QuadFaceVerts &QFV);

	typename exa_set<TriFaceVerts>::iterator getTriVerts(
			exa_set<TriFaceVerts> &vertsOnTris,
			const int face,
			bool& shouldErase);
public:
	CellDivider(UMesh *pVolMesh, const emInt segmentsPerEdge) :
			m_pMesh(pVolMesh), m_Map(nullptr),
					numTriFaces(0), numQuadFaces(0), numEdges(0),
					numVerts(0), nDivs(segmentsPerEdge) {
		localVerts = new emInt[MAX_DIVS + 1][MAX_DIVS + 1][MAX_DIVS + 1];
		m_uvw = new double[MAX_DIVS + 1][MAX_DIVS + 1][MAX_DIVS + 1][3];
#ifndef NDEBUG
		for (int ii = 0; ii <= MAX_DIVS; ii++) {
			for (int jj = 0; jj <= MAX_DIVS; jj++) {
				for (int kk = 0; kk <= MAX_DIVS; kk++) {
					localVerts[ii][jj][kk] = 100;
					m_uvw[ii][jj][kk][0] =
							m_uvw[ii][jj][kk][1] =
									m_uvw[ii][jj][kk][2] = -1;
				}
			}
		}
#endif
	}
	virtual ~CellDivider() {
		delete[] localVerts;
		delete[] m_uvw;
		if (m_Map) delete m_Map;
	}
	void createDivisionVerts(exa_map<Edge, EdgeVerts> &vertsOnEdges,
			exa_set<TriFaceVerts> &vertsOnTris,
			exa_set<QuadFaceVerts> &vertsOnQuads) {
		divideEdges(vertsOnEdges);

		// Divide all the faces, including storing info about which new verts
		// are on which faces
		divideFaces(vertsOnTris, vertsOnQuads);

		// Divide the cell
		divideInterior();

//		printAllPoints();
	}
	void divideEdges(exa_map<Edge, EdgeVerts> &vertsOnEdges);
	void divideFaces(exa_set<TriFaceVerts> &vertsOnTris,
	exa_set<QuadFaceVerts> &vertsOnQuads);
	virtual void divideInterior();
	virtual void createNewCells() = 0;
	virtual void setupCoordMapping(const emInt verts[]) = 0;
	virtual void getPhysCoordsFromParamCoords(const double uvw[],
			double xyz[]) = 0;
	void getParamCoords(const int i, const int j, const int k,
			double uvw[]) {
		assert(i >= 0 && i <= MAX_DIVS);
		assert(j >= 0 && j <= MAX_DIVS);
		assert(k >= 0 && k <= MAX_DIVS);
		uvw[0] = m_uvw[i][j][k][0];
		uvw[1] = m_uvw[i][j][k][1];
		uvw[2] = m_uvw[i][j][k][2];
	}

	// The virtual functions here will be overridden in most subclasses.
	int minI(const int /*j*/, const int /*k*/) const {return 0;}
	int minJ(const int /*i*/, const int /*k*/) const {return 0;}
	virtual int minK(const int /*i*/, const int /*j*/) const {return 0;}

	virtual int maxI(const int /*j*/, const int /*k*/) const {return nDivs;}
	virtual int maxJ(const int /*i*/, const int /*k*/) const {return nDivs;}
	virtual int maxK(const int /*i*/, const int /*j*/) const {return nDivs;}

	virtual int getMinInteriorDivs() const {return 3;}
	void computeParaCoords(const int ii, const int jj, const int kk,
			double uvw[3]) const;
	// Output for diagnostic purposes
	void printAllPoints();
private:
	void getEdgeParametricDivision(EdgeVerts &EV) const;
	void initPerimeterParams(TriFaceVerts& TFV, const int face) const;
	void initPerimeterParams(QuadFaceVerts& QFV, const int face) const;
	bool isEdgeForwardForFace(const EdgeVerts &EV,
			emInt cornerStart, emInt cornerEnd) const;
};

void getFaceParametricIntersectionPoint(
		const double uvL[2], const double uvR[2],
		const double uvB[2], const double uvT[2],
		double uv[2]);

void getCellInteriorParametricIntersectionPoint(
		const double uvwA[3], const double uvwB[3],
		const double uvwC[3], const double uvwD[3],
		const double uvwE[3], const double uvwF[3],
		double uvw[3]);

#endif /* SRC_CELLDIVIDER_H_ */
