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
 * exa-defs.h
 *
 *  Created on: Oct. 24, 2019
 *      Author: cfog
 */

#ifndef SRC_EXA_DEFS_H_
#define SRC_EXA_DEFS_H_

#include <cmath>
#include <stdint.h>
#include <limits.h>
#include <assert.h>

#include "exa_config.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_ORDERED

#include <set>
#include <map>
#define exa_set std::set
#define exa_map std::map
#define exaMultiMap std::multi_map

#else

#include <unordered_set>
#include <unordered_map>

#define exa_set std::unordered_set
#define exa_map std::unordered_map
#define exa_multimap std::unordered_multimap

#endif

#undef PROFILE
#ifdef PROFILE
#include <valgrind/callgrind.h>
#else
#define CALLGRIND_TOGGLE_COLLECT
#endif

#define MAX_DIVS 50
#define FILE_NAME_LEN 1024

typedef uint32_t emInt;
#define EMINT_MAX UINT_MAX

#if (HAVE_CGNS == 0)
#define TRI_3 5
#define QUAD_4 7
#define TETRA_4 10
#define PYRA_5 12
#define PENTA_6 14
#define HEXA_8 17
#define TRI_10 26
#define QUAD_16 28
#define TETRA_20 30
#define PYRA_30 33
#define PENTA_40 36
#define HEXA_64 39
#endif

// Some vector operators

#define DIFF(a,b) {a[0] -b[0], a[1]-b[1],a[2]-b[2]}
#define SCALE(x, a, y) y[0]=(a)*x[0]; y[1]=(a)*x[1]; y[2]=(a)*x[2]
#define LEN(x) sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
#define CROSS(a,b,c) c[0] = a[1]*b[2]-a[2]*b[1], c[1] = a[2]*b[0]-a[0]*b[2], c[2] = a[0]*b[1]-a[1]*b[0]
#define DOT(a,b) (a[0]*b[0] + a[1]*b[1]+ a[2]*b[2])
#define NORMALIZE(a) {double tmp = 1./sqrt(DOT(a,a)); a[0]*=tmp; a[1]*=tmp; a[2]*=tmp;}

inline double safe_acos(const double arg) {
	if (arg < -1) return M_PI;
	else if (arg > 1) return 0;
	else return acos(arg);
}

inline double exaTime() {
#ifdef _OPENMP
	return omp_get_wtime();
#else
	return clock() / double(CLOCKS_PER_SEC);
#endif
}

class Edge {
private:
	emInt v0, v1;
public:
	Edge(emInt vA, emInt vB) {
		if (vA < vB) {
			v0 = vA;
			v1 = vB;
		}
		else {
			v0 = vB;
			v1 = vA;
		}
	}
	bool operator<(const Edge &E) const {
		return (v0 < E.v0 || (v0 == E.v0 && v1 < E.v1));
	}
	bool operator==(const Edge &E) const {
		return v0 == E.v0 && v1 == E.v1;
	}
	emInt getV0() const {
		return v0;
	}

	emInt getV1() const {
		return v1;
	}
};

struct EdgeVerts {
	emInt m_verts[MAX_DIVS + 1];
	double m_param_t[MAX_DIVS + 1];
	double m_totalDihed;
};

class FaceVerts {
protected:
	emInt m_corners[4], m_sorted[4];
	double m_cornerUVW[4][3];
	int m_nCorners, m_nDivs;
	emInt m_intVerts[MAX_DIVS - 1][MAX_DIVS - 1];
	double m_param_st[MAX_DIVS + 1][MAX_DIVS + 1][2];
	double m_param_uvw[MAX_DIVS + 1][MAX_DIVS + 1][3];
	emInt m_volElem, m_volElemType;
	bool m_bothSidesDone;
public:
	FaceVerts(const int nDivs, const emInt NC = 0) :
		m_nCorners(NC), m_nDivs(nDivs), m_volElem(EMINT_MAX), m_volElemType(0),
		m_bothSidesDone(false) {
		assert(NC == 3 || NC == 4);
	}
	virtual ~FaceVerts() {}
	bool isValidIJ(const int ii, const int jj) const {
		bool retVal = (ii >= 0) && (jj >= 0);
		if (m_nCorners == 3) {
			retVal = retVal && ((ii + jj) <= m_nDivs);
		}
		else {
			assert(m_nCorners == 4);
			retVal = retVal && (ii <= m_nDivs) && (jj <= m_nDivs);
		}
		return retVal;
	}
	bool isValidParam(const double param) const {
		// This isn't comprehensive, in that not all faces
		// can have this full parameter range.  But the more
		// accurate test requires significantly more information.
		return (param >= 0 && param <= 1);
	}
	virtual void setupSorted() = 0;
	emInt getSorted(const int ii) const {
		return m_sorted[ii];
	}
	void setIntVertInd(const int ii, const int jj, const emInt vert) {
		assert(isValidIJ(ii, jj));
		m_intVerts[ii][jj] = vert;
	}
	emInt getIntVertInd(const int ii, const int jj) const
	{
		assert(isValidIJ(ii, jj));
		return m_intVerts[ii][jj];
	}
	void setVertSTParams(const int ii, const int jj, const double st[2]){
		assert(isValidIJ(ii, jj));
		assert(isValidParam(st[0]));
		assert(isValidParam(st[1]));
		m_param_st[ii][jj][0] = st[0];
		m_param_st[ii][jj][1] = st[1];
	}
	void getVertSTParams(const int ii, const int jj, double st[2]) const {
		assert(isValidIJ(ii, jj));
		st[0] = m_param_st[ii][jj][0];
		st[1] = m_param_st[ii][jj][1];
		assert(isValidParam(st[0]));
		assert(isValidParam(st[1]));
	}
	void setVertUVWParams(const int ii, const int jj, const double uvw[3]){
		assert(isValidIJ(ii, jj));
		assert(isValidParam(uvw[0]));
		assert(isValidParam(uvw[1]));
		assert(isValidParam(uvw[2]));
		m_param_uvw[ii][jj][0] = uvw[0];
		m_param_uvw[ii][jj][1] = uvw[1];
		m_param_uvw[ii][jj][2] = uvw[2];
	}
	void getVertUVWParams(const int ii, const int jj, double uvw[3]) const {
		assert(isValidIJ(ii, jj));
		uvw[0] = m_param_uvw[ii][jj][0];
		uvw[1] = m_param_uvw[ii][jj][1];
		uvw[2] = m_param_uvw[ii][jj][2];
		assert(isValidParam(uvw[0]));
		assert(isValidParam(uvw[1]));
		assert(isValidParam(uvw[2]));
	}
	void setCorners(const emInt cA, const emInt cB, const emInt cC,
			const emInt cD = EMINT_MAX) {
		m_corners[0] = cA;
		m_corners[1] = cB;
		m_corners[2] = cC;
		m_corners[3] = cD;
		setupSorted();
	}
	emInt getCorner(const int ii) const {
		assert(ii >= 0 && ii < m_nCorners);
		return m_corners[ii];
	}
	virtual void computeParaCoords(const int ii, const int jj,
			double &s, double &t) const = 0;
	emInt getVolElement() const {
		return m_volElem;
	}

	emInt getVolElementType() const {
		return m_volElemType;
	}
};

class TriFaceVerts : public FaceVerts {
//	emInt m_corners[3], m_sorted[3];
//	volatile emInt (*m_intVerts)[MAX_DIVS - 2];
//	double (*m_intParam_st)[MAX_DIVS-2];
//	emInt volElement, volElementType;
public:
	TriFaceVerts(const int nDivs) : FaceVerts(nDivs, 3) {}
	TriFaceVerts(const int nDivs, emInt v0, const emInt v1, const emInt v2,
			const emInt type = 0, const emInt elemInd = EMINT_MAX);
	virtual ~TriFaceVerts() {}
//	void allocVertMemory() {
//		m_intVerts = new emInt[MAX_DIVS - 2][MAX_DIVS - 2];
//	}
//	void freeVertMemory() const {
//		if (m_intVerts) delete[] m_intVerts;
//	}
	virtual void computeParaCoords(const int ii, const int jj,
			double &s, double &t) const;
	virtual void setupSorted();
	friend bool operator<(const TriFaceVerts& a, const TriFaceVerts& b);
	friend bool operator==(const TriFaceVerts& a, const TriFaceVerts& b);
};

struct QuadFaceVerts : public FaceVerts {
//	emInt m_corners[4], m_sorted[4];
//	emInt m_intVerts[MAX_DIVS - 1][MAX_DIVS - 1];
//	emInt volElement, volElementType;
public:
	QuadFaceVerts(const int nDivs) : FaceVerts(nDivs, 4) {}
	QuadFaceVerts(const int nDivs, const emInt v0, const emInt v1, const emInt v2, const emInt v3,
			const emInt type = 0, const emInt elemInd = EMINT_MAX);
	virtual ~QuadFaceVerts() {}
	virtual void computeParaCoords(const int ii, const int jj,
			double &s, double &t) const;
	virtual void setupSorted();
	friend bool operator<(const QuadFaceVerts& a, const QuadFaceVerts& b);
	friend bool operator==(const QuadFaceVerts& a, const QuadFaceVerts& b);
};

namespace std {
	template<> struct hash<TriFaceVerts> {
		typedef TriFaceVerts argument_type;
		typedef std::size_t result_type;
		result_type operator()(const argument_type& TFV) const noexcept
		{
			const result_type h0 = TFV.getSorted(0);
			const result_type h1 = TFV.getSorted(1);
			const result_type h2 = TFV.getSorted(2);
			return (h0 ^ (h1 << 1)) ^ (h2 << 2);
		}
	};

	template<> struct hash<QuadFaceVerts> {
		typedef QuadFaceVerts argument_type;
		typedef std::size_t result_type;
		result_type operator()(const argument_type& QFV) const noexcept
		{
			const result_type h0 = QFV.getSorted(0);
			const result_type h1 = QFV.getSorted(1);
			const result_type h2 = QFV.getSorted(2);
			const result_type h3 = QFV.getSorted(3);
			return h0 ^ (h1 << 1) ^ (h2 << 2) ^ (h3 << 3);
		}
	};

	template<> struct hash<Edge> {
		typedef Edge argument_type;
		typedef std::size_t result_type;
		result_type operator()(const argument_type& E) const noexcept
		{
			const result_type h0 = E.getV0();
			const result_type h1 = E.getV1();
			return (h0 ^ (h1 << 1));
		}
	};
}

bool operator==(const TriFaceVerts& a, const TriFaceVerts& b);
bool operator<(const TriFaceVerts& a, const TriFaceVerts& b);
bool operator==(const QuadFaceVerts& a, const QuadFaceVerts& b);
bool operator<(const QuadFaceVerts& a, const QuadFaceVerts& b);

struct RefineStats {
	double refineTime, extractTime;
	emInt cells;
	size_t fileSize;
};

#endif /* SRC_EXA_DEFS_H_ */
