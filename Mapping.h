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
	double getIsoLengthScale(const emInt vertInd);
	virtual void setupCoordMapping(const emInt verts[]) = 0;
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const = 0;
	enum MappingType {
		Uniform, LengthScale, Lagrange, Compound, Invalid
	};
};

// This family of mappings just does linear / bilinear / trilinear interpolation
// of coordinates.
class Q1Mapping: public Mapping {
protected:
	Q1Mapping(const ExaMesh* const EM) :
			Mapping(EM) {
	}
	virtual ~Q1Mapping() {
	}
};

class Q1TetMapping: public Q1Mapping {
private:
	double A[3], dU[3], dV[3], dW[3];
public:
	Q1TetMapping(const ExaMesh* const EM) :
			Q1Mapping(EM) {
	}
	void setupCoordMapping(const emInt verts[]);
	void computeTransformedCoords(const double uvw[3], double xyz[3]) const;
};

class Q1PyramidMapping: public Q1Mapping {
private:
	double A[3], dU[3], dV[3], dUV[3], dW[3], Apex[3];
public:
	Q1PyramidMapping(const ExaMesh* const EM) :
			Q1Mapping(EM) {
	}
	void setupCoordMapping(const emInt verts[]);
	void computeTransformedCoords(const double uvw[3], double xyz[3]) const;
};

class Q1PrismMapping: public Q1Mapping {
private:
	double A[3], dU[3], dV[3], dW[3], dUW[3], dVW[3];
public:
	Q1PrismMapping(const ExaMesh* const EM) :
			Q1Mapping(EM) {
	}
	void setupCoordMapping(const emInt verts[]);
	void computeTransformedCoords(const double uvw[3], double xyz[3]) const;
};

class Q1HexMapping: public Q1Mapping {
private:
	double A[3], dU[3], dV[3], dW[3], dUV[3], dUW[3], dVW[3], dUVW[3];
public:
	Q1HexMapping(const ExaMesh* const EM) :
			Q1Mapping(EM) {
	}
	void setupCoordMapping(const emInt verts[]);
	void computeTransformedCoords(const double uvw[3], double xyz[3]) const;
};

class Q1TriMapping: public Q1Mapping {
private:
	double A[3], dU[3], dV[3], dW[3], dUV[3], dUW[3], dVW[3], dUVW[3];
public:
	Q1TriMapping(const ExaMesh* const EM) :
			Q1Mapping(EM) {
	}
	void setupCoordMapping(const emInt verts[]);
	void computeTransformedCoords(const double uvw[3], double xyz[3]) const;
};

class Q1QuadMapping: public Q1Mapping {
private:
	double A[3], dU[3], dV[3], dW[3], dUV[3], dUW[3], dVW[3], dUVW[3];
public:
	Q1QuadMapping(const ExaMesh* const EM) :
			Q1Mapping(EM) {
	}
	void setupCoordMapping(const emInt verts[]) {
		assert(0);
	}
	void computeTransformedCoords(const double uvw[3], double xyz[3]) const {
		assert(0);
	}
};


class CompoundMapping: public Mapping {
	// A compound mapping starts from a uniform subdivision of a reference
	// element.  The first stage maps that to a non-uniform subdivision of
	// a reference element, based on the length scales at corners.  The
	// second stage maps those intermediate points to the physical space;
	// this mapping can be linear (Q1Mapping) or non-linear
	// (LagrangeMapping).  Single mappings (Uniform / Lagrange) are
	// equivalent to a CompoundMapping with an identity as the first stage.
protected:
	Mapping *Mapping_RefToInter, *Mapping_InterToPhysical;
public:
	CompoundMapping(Mapping* R2I, Mapping *I2P) : Mapping(nullptr) {
		Mapping_RefToInter = R2I;
		Mapping_InterToPhysical = I2P;
	}
	virtual void setupCoordMapping(const emInt verts[]) {
		Mapping_InterToPhysical->setupCoordMapping(verts);
	}
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const {
		double inter[3];
		Mapping_RefToInter->computeTransformedCoords(uvw, inter);
		Mapping_InterToPhysical->computeTransformedCoords(inter, xyz);
	}
};

class LengthScaleMapping: public Mapping {
protected:
	LengthScaleMapping(const ExaMesh* const EM) :
			Mapping(EM) {
	}
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const = 0;
public:
	virtual ~LengthScaleMapping() {
	}
};

class LengthScaleTetMapping: public LengthScaleMapping {
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
	LengthScaleTetMapping(const ExaMesh * const mesh) :
			LengthScaleMapping(mesh) {
	}
	virtual ~LengthScaleTetMapping() {
	}
	// setPolyCoeffs is public for testing purposes only.
	void setPolyCoeffs(const double xyz[][3],
			double uderiv[][3], double vderiv[][3], double wderiv[][3]);
	virtual void setupCoordMapping(const emInt verts[]);
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const;
};

class LengthScalePyramidMapping: public LengthScaleMapping {
	// This is all set up for a range of (u,v,w) in (0,1).
	// u and v still run from 0 to 1 for all w.
	double xyzOffset[3], uVec[3], vVec[3], wVec[3];
	// For the quad at the base (Hermite quad interpolation)
	double A[3]; // coeff for 1
	double B[3]; // coeff for u
	double C[3]; // coeff for v
	double E[3]; // coeff for u^2
	double F[3]; // coeff for u v
	double G[3]; // coeff for v^2
	double H[3]; // coeff for u^3
	double J[3]; // coeff for u^2 v
	double K[3]; // coeff for u v^2
	double L[3]; // coeff for v^3
	double M[3]; // coeff for u^3 v
	double N[3]; // coeff for u v^3

	double P[3]; // coeff for apex value

	// Bilinear interpolation of w deriv on the base
	double Q[3]; // coeff for 1
	double R[3]; // coeff for u
	double S[3]; // coeff for v
	double T[3]; // coeff for uv

	double U[3]; // coeff for u deriv at the apex
	double V[3]; // coeff for v deriv at the apex
	double W[3]; // coeff for w deriv at the apex
public:
	LengthScalePyramidMapping(const ExaMesh * const mesh) :
			LengthScaleMapping(mesh) {
	}
	virtual ~LengthScalePyramidMapping() {
	}
	// setPolyCoeffs is public for testing purposes only.
	void setPolyCoeffs(const double xyz[][3],
			const double uderiv[5][3], const double vderiv[5][3],
					const double wderiv[5][3]);
	virtual void setupCoordMapping(const emInt verts[]);
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const;
};

class LengthScalePrismMapping: public LengthScaleMapping {
	// Interpolation on the ends using a formulation that interpolates on the
	// end surfaces of the prism (cubic interpolation w/o uv term), and then uses
	// Hermite interpolation  along the sides.
	double xyzOffset[3], uVec[3], vVec[3], wVec[3];
	// These coefficients are for interpolation on the triangle at the
	// bottom of the prism.
	double Cb[3], Cbu[3], Cbv[3], Cbu2[3], Cbuv[3], Cbv2[3];
	double Cbu3[3], Cbu2v[3], Cbuv2[3], Cbv3[3];
	// These coefficients are for interpolation on the triangle at the
	// top of the prism.
	double Ct[3], Ctu[3], Ctv[3], Ctu2[3], Ctuv[3], Ctv2[3];
	double Ctu3[3], Ctu2v[3], Ctuv2[3], Ctv3[3];
	// These coefficients are for interpolating the gradient in the w-direction
	// across the bottom surface.
	double Gb[3], Gbu[3], Gbv[3];
	// And for gradient interpolation on the top surface.
	double Gt[3], Gtu[3], Gtv[3];
public:
	LengthScalePrismMapping(const ExaMesh * const mesh) :
			LengthScaleMapping(mesh) {
	}
	virtual ~LengthScalePrismMapping() {
	}
	// setPolyCoeffs is public for testing purposes only.
	void setPolyCoeffs(const double xyz[6][3], double uderiv[6][3], double vderiv[6][3],
			double wderiv[6][3]);
	virtual void setupCoordMapping(const emInt verts[]);
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const;
};

class LengthScaleHexMapping: public LengthScaleMapping {
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
	double N[3]; // coeff for u v
	double P[3]; // coeff for v^2
	double Q[3]; // coeff for v w
	double R[3]; // coeff for w^2
	double S[3]; // coeff for w u
	double T[3]; // coeff for constant
	double U[3]; // coeff for u
	double V[3]; // coeff for v
	double W[3]; // coeff for w

	double AA[3]; // coeff for u^3 v
	double BB[3]; // coeff for u v^3
	double CC[3]; // coeff for v^3 w
	double DD[3]; // coeff for v w^3
	double EE[3]; // coeff for w^3 u
	double FF[3]; // coeff for w u^3
	double GG[3]; // coeff for u^2 v w
	double HH[3]; // coeff for u v^2 w
	double II[3]; // coeff for u v w^2
	double JJ[3]; // coeff for u^3 v w
	double KK[3]; // coeff for u v^3 w
	double LL[3]; // coeff for u v w^3
public:
	LengthScaleHexMapping(const ExaMesh * const mesh) :
			LengthScaleMapping(mesh) {
	}
	virtual ~LengthScaleHexMapping() {
	}
	// setPolyCoeffs is public for testing purposes only.
	void setPolyCoeffs(const double xyz[8][3],
			 double uderiv[8][3], double vderiv[8][3],
			double wderiv[8][3]);
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

class LagrangeCubicTriMapping: public LagrangeCubicMapping {
	double C[3], Cu[3], Cv[3], Cu2[3], Cuv[3], Cv2[3], Cu3[3], Cv3[3],
		Cu2v[3], Cuv2[3];
public:
	LagrangeCubicTriMapping(const ExaMesh* const EM) :
			LagrangeCubicMapping(EM, 10) {
	}
	virtual ~LagrangeCubicTriMapping() {
	}
	// Public for testing purposes
	void setModalValues();
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const;

};

class LagrangeCubicQuadMapping: public LagrangeCubicMapping {
	double C[3], Cu[3], Cv[3], Cu2[3], Cuv[3], Cv2[3], Cu3[3], Cv3[3],
		Cu2v[3], Cuv2[3], Cu3v[3], Cuv3[3], Cu2v2[3], Cu3v2[3], Cu2v3[3],
			Cu3v3[3];
public:
	LagrangeCubicQuadMapping(const ExaMesh* const EM) :
			LagrangeCubicMapping(EM, 16) {
	}
	virtual ~LagrangeCubicQuadMapping() {
	}
	// Public for testing purposes
	void setModalValues()  {
		assert(0);
	}
	virtual void computeTransformedCoords(const double uvw[3],
			double xyz[3]) const {
		assert(0);
	}

};

#endif /* SRC_MAPPING_H_ */
