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
		Uniform,
		Lagrange,
		Invalid
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
