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

#define MAX_DIVS 50
#define FILE_NAME_LEN 1024

typedef uint32_t emInt;
#define EMINT_MAX UINT_MAX

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
	emInt verts[MAX_DIVS + 1];
	double m_totalDihed;
};

struct TriFaceVerts {
	emInt corners[3], sorted[3];
	volatile emInt (*intVerts)[MAX_DIVS - 2];
//	emInt intVerts[MAX_DIVS - 2][MAX_DIVS - 2];
	emInt volElement, volElementType;
	TriFaceVerts() :
			intVerts(nullptr), volElement(EMINT_MAX), volElementType(0) {
	}
	TriFaceVerts(const emInt v0, const emInt v1, const emInt v2,
			const emInt type = 0, const emInt elemInd = EMINT_MAX);
	~TriFaceVerts() {
	}
	void allocVertMemory() {
		intVerts = new emInt[MAX_DIVS - 2][MAX_DIVS - 2];
	}
	void freeVertMemory() const {
		// TODO  It's really stinky to do this to avoid a double free....
		if (intVerts) delete[] intVerts;
//		intVerts = nullptr;
	}
	void setupSorted();
};

struct QuadFaceVerts {
	emInt corners[4], sorted[4];
	emInt intVerts[MAX_DIVS - 1][MAX_DIVS - 1];
	emInt volElement, volElementType;
	QuadFaceVerts() :
			volElement(EMINT_MAX), volElementType(0) {
	}
	QuadFaceVerts(const emInt v0, const emInt v1, const emInt v2, const emInt v3,
			const emInt type = 0, const emInt elemInd = EMINT_MAX);
	void setupSorted();
};

namespace std {
	template<> struct hash<TriFaceVerts> {
		typedef TriFaceVerts argument_type;
		typedef std::size_t result_type;
		result_type operator()(const argument_type& TFV) const noexcept
		{
			const result_type h0 = TFV.sorted[0];
			const result_type h1 = TFV.sorted[1];
			const result_type h2 = TFV.sorted[2];
			return (h0 ^ (h1 << 1)) ^ (h2 << 2);
		}
	};

	template<> struct hash<QuadFaceVerts> {
		typedef QuadFaceVerts argument_type;
		typedef std::size_t result_type;
		result_type operator()(const argument_type& QFV) const noexcept
		{
			const result_type h0 = QFV.sorted[0];
			const result_type h1 = QFV.sorted[1];
			const result_type h2 = QFV.sorted[2];
			const result_type h3 = QFV.sorted[3];
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

enum {
	UNKNOWN, TRIANGLE, QUADRILATERAL, TETRAHEDRON, PYRAMID, PRISM, HEXAHEDRON
};

class TriFaceToCell {
	emInt m_corners[3];
	emInt m_cellLeftIndex, m_cellRightIndex;
	emInt m_cellLeftType, m_cellRightType;
	emInt m_refinedVerts[MAX_DIVS][MAX_DIVS];
	bool isRefined;
public:
	TriFaceToCell() :
			m_cellLeftIndex(EMINT_MAX), m_cellRightIndex(EMINT_MAX),
					m_cellLeftType(EMINT_MAX), m_cellRightType(EMINT_MAX),
					isRefined(false) {
	}
//	TriFaceToCell(const emInt corners[3], const emInt leftInd,
//			const emInt rightInd, const emInt leftType, const emInt rightType) :
//			m_cellLeftIndex(leftInd), m_cellRightIndex(rightInd),
//					m_cellLeftType(leftType), m_cellRightType(rightType), isRefined(false) {
//		m_corners[0] = corners[0];
//		m_corners[1] = corners[1];
//		m_corners[2] = corners[2];
//	}
	void init(const emInt corners[], const emInt leftInd, const emInt rightInd,
			const emInt leftType, const emInt rightType) {
		m_corners[0] = corners[0];
		m_corners[1] = corners[1];
		m_corners[2] = corners[2];
		m_cellLeftIndex = leftInd;
		m_cellRightIndex = rightInd;
		m_cellLeftType = leftType;
		m_cellRightType = rightType;
		isRefined = false;
	}
};

class CellInfo {
	emInt m_index, m_type, m_whichFace;
public:
	CellInfo(const emInt ind, const emInt type, const emInt face) :
			m_index(ind), m_type(type), m_whichFace(face) {
	}

	emInt getIndex() const {
		return m_index;
	}

	emInt getType() const {
		return m_type;
	}

	emInt getWhichFace() const {
		return m_whichFace;
	}
};
#endif /* SRC_EXA_DEFS_H_ */
