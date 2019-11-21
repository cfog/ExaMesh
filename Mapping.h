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
	double A[3], dU[3], dV[3], dUV[3], dW[3], Apex[3];
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
			double xyz[3]) const = 0;
};

class LagrangeCubicMapping: public LagrangeMapping {
protected:
	LagrangeCubicMapping(const ExaMesh* const EM, const emInt nDOFs) :
			LagrangeMapping(EM, nDOFs) {
	}
public:
private:
	virtual void setModalValues() = 0;
};

class LagrangeCubicTetMapping: public LagrangeCubicMapping {
	double C[3], Cu[3], Cv[3], Cw[3], Cuu[3], Cuv[3], Cuw[3], Cvv[3], Cvw[3],
			Cww[3], Cuuu[3], Cuuv[3], Cuvv[3], Cuuw[3], Cuww[3], Cuvw[3], Cvvv[3],
			Cvvw[3], Cvww[3], Cwww[3];
public:
	LagrangeCubicTetMapping(const ExaMesh* const EM) :
			LagrangeCubicMapping(EM, 20) {
	}
	virtual ~LagrangeCubicTetMapping() {
	}
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const;

	// Public for testing purposes
	void setModalValues();
};

class LagrangeCubicPyramidMapping: public LagrangeCubicMapping {
	double C[3], Cu[3], Cv[3], Cw[3], Cuu[3], Cuv[3], Cuw[3], Cvv[3], Cvw[3],
			Cww[3], Cuuu[3], Cuuv[3], Cuvv[3], Cuuw[3], Cuww[3], Cuvw[3], Cvvv[3],
			Cvvw[3], Cvww[3], Cwww[3];
	double CuvOverw[3], Cu2vOverw[3], Cuv2Overw[3], Cu3vOverw[3], Cu2v2Overw[3],
			Cuv3Overw[3];
	double Cu2v2Overw2[3], Cu3v2Overw2[3], Cu2v3Overw2[3], Cu3v3Overw3[3];
	double Apex[3];
public:
	LagrangeCubicPyramidMapping(const ExaMesh* const EM) :
			LagrangeCubicMapping(EM, 30) {
	}
	virtual ~LagrangeCubicPyramidMapping() {
	}
	// Public for testing purposes
	void setModalValues();
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const;

};

class LagrangeCubicPrismMapping: public LagrangeCubicMapping {
	double C0[3], Cu0[3], Cv0[3], Cuu0[3], Cuv0[3], Cvv0[3], Cuuu0[3], Cuuv0[3],
			Cuvv0[3], Cvvv0[3];
	double C1[3], Cu1[3], Cv1[3], Cuu1[3], Cuv1[3], Cvv1[3], Cuuu1[3], Cuuv1[3],
			Cuvv1[3], Cvvv1[3];
	double C2[3], Cu2[3], Cv2[3], Cuu2[3], Cuv2[3], Cvv2[3], Cuuu2[3], Cuuv2[3],
			Cuvv2[3], Cvvv2[3];
	double C3[3], Cu3[3], Cv3[3], Cuu3[3], Cuv3[3], Cvv3[3], Cuuu3[3], Cuuv3[3],
			Cuvv3[3], Cvvv3[3];
public:
	LagrangeCubicPrismMapping(const ExaMesh* const EM) :
			LagrangeCubicMapping(EM, 40) {
	}
	virtual ~LagrangeCubicPrismMapping() {
	}
	// Public for testing purposes
	void setModalValues();
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const;

};

class LagrangeCubicHexMapping: public LagrangeCubicMapping {
	double C[3], Cu[3], Cv[3], Cw[3], Cu2[3], Cuv[3], Cv2[3], Cvw[3], Cw2[3],
			Cuw[3], Cu3[3], Cv3[3], Cw3[3], Cu2v[3], Cuv2[3], Cv2w[3], Cvw2[3],
			Cu2w[3], Cuw2[3], Cuvw[3], Cu3v[3], Cu3w[3], Cuv3[3], Cv3w[3], Cuw3[3],
			Cvw3[3], Cu2v2[3], Cu2w2[3], Cv2w2[3], Cu2vw[3], Cuv2w[3], Cuvw2[3],
			Cu3vw[3], Cuv3w[3], Cuvw3[3], Cu3v2[3], Cu3w2[3], Cu2v3[3], Cv3w2[3],
			Cu2w3[3], Cv2w3[3], Cu2v2w[3], Cu2vw2[3], Cuv2w2[3], Cu3v3[3], Cu3w3[3],
			Cv3w3[3], Cu3v2w[3], Cu3vw2[3], Cu2v3w[3], Cuv3w2[3], Cu2vw3[3],
			Cuv2w3[3], Cu2v2w2[3], Cu3v2w2[3], Cu2v3w2[3], Cu2v2w3[3], Cu3v3w[3],
			Cu3vw3[3], Cuv3w3[3], Cu3v3w2[3], Cu3v2w3[3], Cu2v3w3[3], Cu3v3w3[3];
public:
	LagrangeCubicHexMapping(const ExaMesh* const EM) :
			LagrangeCubicMapping(EM, 64) {
	}
	virtual ~LagrangeCubicHexMapping() {
	}
	// Public for testing purposes
	void setModalValues();
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const;

};

#endif /* SRC_MAPPING_H_ */
