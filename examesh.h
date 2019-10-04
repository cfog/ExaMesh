/*
 * examesh.h
 *
 *  Created on: May 17, 2019
 *      Author: cfog
 */

#ifndef SRC_EXAMESH_H_
#define SRC_EXAMESH_H_

#define MAX_DIVS 50
#define FILE_NAME_LEN 1024

#include <stdint.h>
#include <limits.h>
#include <cmath>

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

class ExaMesh {
protected:
	double *m_lenScale;
	void setupLengthScales();

public:
	ExaMesh() :
			m_lenScale(nullptr) {
	}
	virtual ~ExaMesh() {
		if (m_lenScale) delete[] m_lenScale;
	}
	virtual double getX(const emInt vert) const = 0;
	virtual double getY(const emInt vert) const = 0;
	virtual double getZ(const emInt vert) const =0;
	virtual void getCoords(const emInt vert, double coords[3]) const = 0;

	virtual emInt numVerts() const = 0;
	virtual emInt numBdryTris() const = 0;
	virtual emInt numBdryQuads() const = 0;
	virtual emInt numTets() const = 0;
	virtual emInt numPyramids() const = 0;
	virtual emInt numPrisms() const = 0;
	virtual emInt numHexes() const = 0;

	const virtual emInt* getBdryTriConn(const emInt bdryTri) const=0;
	const virtual emInt* getBdryQuadConn(const emInt bdryQuad) const=0;
	const virtual emInt* getTetConn(const emInt tet) const=0;
	const virtual emInt* getPyrConn(const emInt pyr) const=0;
	const virtual emInt* getPrismConn(const emInt prism) const=0;
	const virtual emInt* getHexConn(const emInt hex) const=0;
	void printMeshSizeStats();
	double getLengthScale(const emInt vert) const {
		// TODO  re-add this after it's possible to get the number
		// of verts for a CubicMesh
		// assert(vert < m_nVerts);
		return m_lenScale[vert];
	}

};

struct MeshSize {
	emInt nBdryVerts, nVerts, nBdryTris, nBdryQuads, nTets, nPyrs, nPrisms,
			nHexes;
};

bool computeMeshSize(const struct MeshSize& MSIn, const emInt nDivs,
		struct MeshSize& MSOut);

// Defined elsewhere.
class UMesh;

emInt subdividePartMesh(const UMesh * const pVM_input, UMesh * const pVM_output,
		const int nDivs);

void sortVerts3(const emInt input[3], emInt output[3]);
void sortVerts4(const emInt input[4], emInt output[4]);

#endif /* SRC_EXAMESH_H_ */
