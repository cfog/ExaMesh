/*
 * Mapping.h
 *
 *  Created on: Oct. 1, 2019
 *      Author: cfog
 */

#ifndef SRC_MAPPING_H_
#define SRC_MAPPING_H_

#include "exa-defs.h"

class ExaMesh;

// Lagrange element mappings

class Mapping {
protected:
	const ExaMesh* m_pMesh;
	Mapping(const ExaMesh* const EM) :
			m_pMesh(EM) {
	}
public:
	virtual ~Mapping() {
	}
	virtual void setupCoordMapping(const emInt verts[]) = 0;
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const = 0;
	enum MappingType {
		Uniform, LengthScale, Lagrange, Invalid
	};
};

// This family of mappings just does linear / bilinear / trilinear interpolation
// of coordinates.
class UniformMapping: public Mapping {
protected:
	UniformMapping(const ExaMesh* const EM) :
			Mapping(EM) {
	}
	virtual ~UniformMapping() {
	}
};

class UniformTetMapping: public UniformMapping {
private:
	double A[3], dU[3], dV[3], dW[3];
public:
	UniformTetMapping(const ExaMesh* const EM) :
			UniformMapping(EM) {
	}
	void setupCoordMapping(const emInt verts[]);
	void computeTransformedCoords(const double uvw[3], double xyz[3]) const;
};

class UniformPyramidMapping: public UniformMapping {
private:
	double A[3], dU[3], dV[3], dUV[3], Apex[3];
public:
	UniformPyramidMapping(const ExaMesh* const EM) :
			UniformMapping(EM) {
	}
	void setupCoordMapping(const emInt verts[]);
	void computeTransformedCoords(const double uvw[3], double xyz[3]) const;
};

class UniformPrismMapping: public UniformMapping {
private:
	double A[3], dU[3], dV[3], dW[3], dUW[3], dVW[3];
public:
	UniformPrismMapping(const ExaMesh* const EM) :
			UniformMapping(EM) {
	}
	void setupCoordMapping(const emInt verts[]);
	void computeTransformedCoords(const double uvw[3], double xyz[3]) const;
};

class UniformHexMapping: public UniformMapping {
private:
	double A[3], dU[3], dV[3], dW[3], dUV[3], dUW[3], dVW[3], dUVW[3];
public:
	UniformHexMapping(const ExaMesh* const EM) :
			UniformMapping(EM) {
	}
	void setupCoordMapping(const emInt verts[]);
	void computeTransformedCoords(const double uvw[3], double xyz[3]) const;
};

class LengthScaleMapping: public Mapping {
protected:
	double getIsoLengthScale(const emInt vert);
	LengthScaleMapping(const ExaMesh* const EM) :
			Mapping(EM) {
	}
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const = 0;
public:
	virtual ~LengthScaleMapping() {
	}
};

class TetLengthScaleMapping: public LengthScaleMapping {
	double xyzOffset[3], uVec[3], vVec[3], wVec[3];
	double A[3]; // coeff for u^3
	double B[3]; // coeff for u^2 v
	double C[3]; // coeff for u v^2
	double E[3]; // coeff for v^3
	double F[3]; // coeff for v^2 w
	double G[3]; // coeff for v w^2
	double H[3]; // coeff for w^3
	double J[3]; // coeff for w^2 u
	double K[3]; // coeff for w u^2
	double L[3]; // coeff for u v w
	double M[3]; // coeff for u^2
	double P[3]; // coeff for v^2
	double R[3]; // coeff for w^2
	double T[3]; // coeff for constant
	double U[3]; // coeff for u
	double V[3]; // coeff for v
	double W[3]; // coeff for w
public:
	TetLengthScaleMapping(const ExaMesh * const mesh) :
			LengthScaleMapping(mesh) {
	}
	virtual ~TetLengthScaleMapping() {
	}
	// setPolyCoeffs is public for testing purposes only.
	void setPolyCoeffs(const double* xyz0, const double* xyz1, const double* xyz2,
			const double* xyz3, double uderiv0[3], double vderiv0[3],
			double wderiv0[3], double uderiv1[3], double vderiv1[3],
			double wderiv1[3], double uderiv2[3], double vderiv2[3],
			double wderiv2[3], double uderiv3[3], double vderiv3[3],
			double wderiv3[3]);
	virtual void setupCoordMapping(const emInt verts[]);
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const;
};

class LagrangeMapping: public Mapping {
	int m_numValues;
	virtual double computeBasisFunction(const int whichFunc,
			const double uvw[3]) const = 0;
	virtual void setModalValues() = 0;
protected:
	double (*m_nodalValues)[3];
public:
	LagrangeMapping(const ExaMesh* const EM, const int nVals);
	virtual ~LagrangeMapping();
	void setupCoordMapping(const emInt verts[]);
	// setNodalValues is public for test purposes only.
	void setNodalValues(double inputValues[][3]);
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const;
};

class LagrangeCubicMapping: public LagrangeMapping {
protected:
	double C[3], Cu[3], Cv[3], Cw[3], Cuu[3], Cuv[3], Cuw[3], Cvv[3], Cvw[3],
			Cww[3], Cuuu[3], Cuuv[3], Cuvv[3], Cuuw[3], Cuww[3], Cuvw[3], Cvvv[3],
			Cvvw[3], Cvww[3], Cwww[3];
	LagrangeCubicMapping(const ExaMesh* const EM, const emInt nDOFs) :
			LagrangeMapping(EM, nDOFs) {
	}
public:
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const;
private:
	virtual void setModalValues() = 0;
};

class LagrangeCubicTetMapping: public LagrangeCubicMapping {
public:
	LagrangeCubicTetMapping(const ExaMesh* const EM) :
			LagrangeCubicMapping(EM, 20) {
	}
	virtual ~LagrangeCubicTetMapping() {
	}
	// Public for testing purposes
	virtual double computeBasisFunction(const int whichFunc,
			const double uvw[3]) const;
	void setModalValues();
};



#endif /* SRC_MAPPING_H_ */
