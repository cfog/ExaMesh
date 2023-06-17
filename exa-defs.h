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
#include <iostream>
#include <cmath>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <algorithm>
#include "exa_config.h"
//#include <mpi.h>
//#include <boost/serialization/access.hpp>
#include <boost/mpi/datatype.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
//#include  <boost/mpi/datatype.hpp>
//#include <boost/serialization/access.hpp>
//#include <boost/mpl/assert.hpp>
#include <boost/mpi.hpp> 
//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp> 
//#include <boost/mpl.hpp>


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
#define TOLTEST 1e-9

typedef int32_t emInt;
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

class FaceVerts
{
private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int /*version*/)
	{
		ar &m_remoteId;
		ar &m_global;
		ar &m_sortedGlobal;
		ar &m_remote;
		ar &m_sortedRemote; 
		ar &m_corners; 
		ar &m_sorted; 
		ar &m_cornerUVW; 
		ar &m_nCorners; 
		ar &m_nDivs; 
		ar &m_intVerts; 
		ar &m_param_st; 
		ar &m_param_uvw; 
		ar &m_volElem; 
		ar &m_volElemType; 
		ar &m_bothSidesDone; 
		ar &m_partId; 
		 
		ar &m_globalComparison; 

	}
protected:
	// CUATION: ANY CHANGE IN MEMBER SIZE OR DELETION OR ADDING SHOULD BE REFLECTED IN ITS MPI TYPE

	emInt m_global[4], m_sortedGlobal[4];
	emInt m_remote[4], m_sortedRemote[4];
	emInt m_corners[4], m_sorted[4];
	double m_cornerUVW[4][3];
	int m_nCorners, m_nDivs;
	emInt m_intVerts[MAX_DIVS - 1][MAX_DIVS - 1];
	double m_param_st[MAX_DIVS + 1][MAX_DIVS + 1][2];
	double m_param_uvw[MAX_DIVS + 1][MAX_DIVS + 1][3];
	emInt m_volElem, m_volElemType;
	bool m_bothSidesDone;
	emInt m_partId; 
	emInt m_remoteId; 
	bool m_globalComparison; 
public:
	 FaceVerts(){};
	FaceVerts(const int nDivs, const emInt NC = 0) : m_nCorners(NC), m_nDivs(nDivs), m_volElem(EMINT_MAX), m_volElemType(0),
													 m_bothSidesDone(false)
	{
		assert(NC == 3 || NC == 4);
	}
	//virtual ~FaceVerts() {};
	bool isValidIJ(const int ii, const int jj) const
	{
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
	//virtual void setupSorted() = 0;
	void setupSorted() {}; 
	emInt getSorted(const int ii) const
	{
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
	//virtual void getVertAndST(const int ii, const int jj, emInt &vert,
		//					  double st[2], const int rotCase = 0) const = 0;
	void setVertUVWParams(const int ii, const int jj, const double uvw[3])
	{
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
	// void setCorners(const emInt cA, const emInt cB, const emInt cC,
	// 				const emInt cD = EMINT_MAX)
	// {
	// 	m_corners[0] = cA;
	// 	m_corners[1] = cB;
	// 	m_corners[2] = cC;
	// 	m_corners[3] = cD;
	// 	setupSorted();
	// }
	void setGlobalCorners(const emInt cA, const emInt cB, const emInt cC,
			const emInt cD = EMINT_MAX) {
		m_global[0] = cA;
		m_global[1] = cB;
		m_global[2] = cC;
		m_global[3] = cD;
		// setupSorted();
	}
	emInt getCorner(const int ii) const {
		assert(ii >= 0 && ii < m_nCorners);
		return m_corners[ii];
	}
	emInt getGlobalCorner(const int ii) const{
		assert(ii >= 0 && ii < m_nCorners);
		return m_global[ii];
	}
	emInt getGlobalSorted(const int ii) const{
		assert(ii >= 0 && ii < m_nCorners);
		return m_sortedGlobal[ii];
	}
	emInt getRemoteIndices(const int ii) const{
		assert(ii >= 0 && ii < m_nCorners);
		return m_remote[ii];
	}
	emInt getSortedRemoteIndices(const int ii) const {
		assert(ii >= 0 && ii < m_nCorners);
		return m_sortedRemote[ii];
	}
	//virtual void computeParaCoords(const int ii, const int jj,
								//   double st[2]) const = 0;
	emInt getVolElement() const
	{
		return m_volElem;
	}

	emInt getVolElementType() const {
		return m_volElemType;
	}
	emInt getPartid () const {
		return m_partId; 
	}
	void setRemotePartID (const emInt remotePartid_){
		m_remoteId=remotePartid_; 
	}
	void setPartID(const emInt partID_){
		m_partId=partID_; 
	}
	// void setRemoteIndices(const emInt remote[3]){
	// 	for(auto i=0; i<3; i++){
	// 		remoteIndices[i]=remote[i]; 
	// 	}

	// }
	void setRemoteIndices(const emInt* remote){
		for(auto i=0; i<4; i++){
			m_remote[i]=*(remote + i); 
		}

	}
	emInt getRemotePartid ()const{
		return m_remoteId; 
	}

	bool getGlobalCompare() const {
		return m_globalComparison; 
	}
	emInt getNumDivs() const
	{
		return m_nDivs;
	}
	void setCompare(const bool compare)
	{
		m_globalComparison=compare; 
	}
	//friend MPI_Datatype register_mpi_type(FaceVerts const &);


};
//Notes from Documentation : 
// It may not be immediately obvious how this one template serves 
//for both saving data to an archive as well as loading data from the archive
// The key is that the & operator is defined as << for output archives and as >> input archives.
class TriFaceVerts : public FaceVerts
{
private: 

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int /*version*/)
  	{
    	ar & boost::serialization::base_object<FaceVerts>(*this);
   
  	}


public:
	TriFaceVerts(){}; // Dangerous! CHANGE IT LATER ON
	TriFaceVerts(const int nDivs, const emInt partID = -1, bool globalComparison = false) : FaceVerts(nDivs, 3) {
		m_partId=partID; 
		m_globalComparison=globalComparison;
	}
	TriFaceVerts(const int nDivs, emInt v0, const emInt v1, const emInt v2,
			const emInt type = 0, const emInt elemInd = EMINT_MAX, 
			const emInt partID=-1,bool globalComparison=false);

	TriFaceVerts(const int nDivs, const emInt local[3], 
	const emInt global[3],const emInt partid_=-1, const emInt remoteID=-1 ,
	const emInt type=0 ,const emInt elemInd=EMINT_MAX,bool globalComparison=false);

	// TriFaceVerts(const int nDivs,const emInt global[3],
	// const emInt partid_=-1, const emInt remoteID=-1 ,const emInt type=0 ,
	// const emInt elemInd=EMINT_MAX,const bool globalComparison=false);

	TriFaceVerts(const int nDivs, const emInt global[3],
				 const emInt partid_ = -1, const emInt remoteID = -1, const bool globalComparison = false,
				 const emInt type = 0,
				 const emInt elemInd = EMINT_MAX);

	TriFaceVerts(const int nDivs, const emInt local[3],
				 const emInt global[3], const emInt remoteIndices_[3], const emInt partid_,
				 const emInt remoteID,
				 const emInt type = 0, const emInt elemInd = EMINT_MAX, bool globalComparison = false);

	//virtual ~TriFaceVerts() {}; 
	//	void allocVertMemory() {
	//		m_intVerts = new emInt[MAX_DIVS - 2][MAX_DIVS - 2];
	//	}
	//	void freeVertMemory() const {
	//		if (m_intVerts) delete[] m_intVerts;
	//	}
	void computeParaCoords(const int ii, const int jj,
								   double st[2]) const;
	// virtual void computeParaCoords(const int ii, const int jj,
	// 							   double st[2]) const;
	//virtual void setupSorted();
	void setupSorted();
	void getVertAndST(const int ii, const int jj, emInt &vert,
					  double st[2], const int rotCase = 0) const;
	void getTrueIJ(const int ii, const int jj,
				   int &trueI, int &trueJ, const int rotCase = 0) const;
	friend bool operator<(const TriFaceVerts &a, const TriFaceVerts &b);
	friend bool operator==(const TriFaceVerts &a, const TriFaceVerts &b);
	friend inline MPI_Datatype register_mpi_type(TriFaceVerts const &);
	void setCorners(const emInt cA, const emInt cB, const emInt cC,
					const emInt cD = EMINT_MAX)
	{
		m_corners[0] = cA;
		m_corners[1] = cB;
		m_corners[2] = cC;
		m_corners[3] = cD;
		setupSorted();
	}
};

struct QuadFaceVerts : public FaceVerts
{
	private: 
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int /*version*/)
  	{
    	ar & boost::serialization::base_object<FaceVerts>(*this);
   
  	}

public:
	QuadFaceVerts(){}; // Dangerous! CHANGE IT LATER ON
	QuadFaceVerts(const int nDivs, const emInt partID = -1,
				  const emInt _remotePartid = -1, bool globalCompare = false) : FaceVerts(nDivs, 4) {
					m_partId=partID; 
					m_remoteId=_remotePartid; 
	 				m_globalComparison=globalCompare; 
				  }
	QuadFaceVerts(const int nDivs, const emInt v0, const emInt v1, const emInt v2, const emInt v3,

				  const emInt type = 0, const emInt elemInd = EMINT_MAX,
				  const emInt partID = -1, const emInt remoteID = -1, bool globalCompare = false);

	QuadFaceVerts(const int nDivs, const emInt local[4],
				  const emInt global[4], const emInt partid_ = -1, const emInt remoteID = -1, const emInt type = 0, const emInt elemInd = EMINT_MAX,
				  bool globalCompare = false);
	// QuadFaceVerts(const int nDivs,const emInt global[4],const emInt partid_=-1, const emInt remoteID=-1
	// ,const emInt type=0 ,const emInt elemInd=EMINT_MAX,
	// bool globalCompare=false);
	QuadFaceVerts(const int nDivs, const emInt global[4], const emInt partid_ = -1,
				  const emInt remoteID = -1, bool globalCompare = false, const emInt type = 0, const emInt elemInd = EMINT_MAX);
	QuadFaceVerts(const int nDivs, const emInt local[4],
				  const emInt global[4], const emInt remotelocal[4], const emInt partid_ = -1, const emInt remoteID = -1, const emInt type = 0, const emInt elemInd = EMINT_MAX,
				  bool globalCompare = false);

	//virtual ~QuadFaceVerts() {}
	void computeParaCoords(const int ii, const int jj,
								   double st[2]) const;
	// virtual void computeParaCoords(const int ii, const int jj,
	// 							   double st[2]) const;
	//virtual void setupSorted();
	void setupSorted();
	void getVertAndST(const int ii, const int jj, emInt &vert,
					  double st[2], const int rotCase = 0) const;
	void getTrueIJ(const int ii, const int jj,
				   int &trueI, int &trueJ, const int rotCase = 0) const;
	friend bool operator<(const QuadFaceVerts &a, const QuadFaceVerts &b);
	friend bool operator==(const QuadFaceVerts &a, const QuadFaceVerts &b);
	friend inline MPI_Datatype register_mpi_type(QuadFaceVerts const &);
	void setCorners(const emInt cA, const emInt cB, const emInt cC,
					const emInt cD = EMINT_MAX)
	{
		m_corners[0] = cA;
		m_corners[1] = cB;
		m_corners[2] = cC;
		m_corners[3] = cD;
		setupSorted();
	}
};

namespace std {
	template<> struct hash<TriFaceVerts> {
		typedef TriFaceVerts argument_type;
		typedef std::size_t result_type;
		result_type operator()(const argument_type& TFV) const noexcept
		{
			if(TFV.getGlobalCompare()==false){
				const result_type h0 = TFV.getSorted(0);
				const result_type h1 = TFV.getSorted(1);
				const result_type h2 = TFV.getSorted(2);
				return (h0 ^ (h1 << 1)) ^ (h2 << 2);
			}if(TFV.getGlobalCompare()==true){
				const result_type h0 = TFV.getGlobalSorted(0);
				const result_type h1 = TFV.getGlobalSorted(1);
				const result_type h2 = TFV.getGlobalSorted(2);
				return (h0 ^ (h1 << 1)) ^ (h2 << 2);
			}

			
		}
	};

	template<> struct hash<QuadFaceVerts> {
		typedef QuadFaceVerts argument_type;
		typedef std::size_t result_type;
		result_type operator()(const argument_type& QFV) const noexcept
		{
			if(QFV.getGlobalCompare()==false){
				const result_type h0 = QFV.getSorted(0);
				const result_type h1 = QFV.getSorted(1);
				const result_type h2 = QFV.getSorted(2);
				const result_type h3 = QFV.getSorted(3);
				return h0 ^ (h1 << 1) ^ (h2 << 2) ^ (h3 << 3);
			}else{
				const result_type h0 = QFV.getGlobalSorted(0);
				const result_type h1 = QFV.getGlobalSorted(1);
				const result_type h2 = QFV.getGlobalSorted(2);
				const result_type h3 = QFV.getGlobalSorted(3);
				return h0 ^ (h1 << 1) ^ (h2 << 2) ^ (h3 << 3);

			}

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
inline void printTris(const exa_set<TriFaceVerts>  &tris, emInt nDivs){
	
	//cout<<"checking for tris: "<<endl; 
	//for(auto i=0 ; i<tris.size(); i++){
		
		std::cout<<"size of set: "<< tris.size()<<std::endl; 
		std::cout<<"-----------------------------------------------------------"<<std::endl; 
		for(auto itr=tris.begin(); itr!=tris.end();itr++){
			std::cout<<"Part: "<< itr->getPartid()<<
			" local indices: "<<itr->getCorner(0)<<" "<<
			itr->getCorner(1)<<" "<<itr->getCorner(2)<<
			// " global: "<<itr->getGlobalSorted(0)<<" "<<
			// itr->getGlobalSorted(1)<<" "<<itr->getGlobalSorted(2)<<
			" Unsorted global: "<<itr->getGlobalCorner(0)<<" "<<
			 itr->getGlobalCorner(1)<<" "<<itr->getGlobalCorner(2)<<
		 	" Remote ID: "<<itr->getRemotePartid()<<
			" Remote Indices: "<<itr->getRemoteIndices(0)<<" "<<
			itr->getRemoteIndices(1)<<" "<<
			itr->getRemoteIndices(2)<<
			" boolean value: "<<itr->getGlobalCompare()<<
			std::endl;
			std::cout<<"Refined verts: "<<std::endl; 
			for (int ii = 0; ii <= nDivs ; ii++) {
	 			for (int jj = 0; jj <= nDivs-ii ; jj++) {

	 				 std::cout<<itr->getIntVertInd(ii,jj)<<" "; 
	 			}
			}
			
			std::cout<<std::endl;
		}
	//}

	
}
inline void printQuads(const exa_set<QuadFaceVerts>  &quads, emInt nDivs){
	
	//cout<<"checking for tris: "<<endl; 
	//for(auto i=0 ; i<tris.size(); i++){
		
		std::cout<<"size of set: "<< quads.size()<<std::endl; 
		std::cout<<"-----------------------------------------------------------"<<std::endl; 
		for(auto itr=quads.begin(); itr!=quads.end();itr++){
			std::cout<<"Part: "<< itr->getPartid()<<
			" local indices: "<<itr->getCorner(0)<<" "<<
			itr->getCorner(1)<<" "<<itr->getCorner(2)<<" "<<itr->getCorner(3)<<
			
			// " global: "<<itr->getGlobalSorted(0)<<" "<<
			// itr->getGlobalSorted(1)<<" "<<itr->getGlobalSorted(2)<<
			" Unsorted global: "<<itr->getGlobalCorner(0)<<" "<<
			 itr->getGlobalCorner(1)<<" "<<itr->getGlobalCorner(2)<<" "<<itr->getGlobalCorner(3)<<
		 	" Remote ID: "<<itr->getRemotePartid()<<
			" Remote Indices: "<<itr->getRemoteIndices(0)<<" "<<
			itr->getRemoteIndices(1)<<" "<<
			itr->getRemoteIndices(2)<<" "<<itr->getRemoteIndices(3)<<
			" boolean value: "<<itr->getGlobalCompare()<<
			std::endl;
			std::cout<<"Refined verts: "<<std::endl; 
			for (int ii = 0; ii <= nDivs ; ii++) {
	 			for (int jj = 0; jj <= nDivs ; jj++) {

	 				 std::cout<<itr->getIntVertInd(ii,jj)<<" "; 
	 			}
			}
			
			std::cout<<std::endl;
		}
	//}
}
inline emInt getTriRotation(const TriFaceVerts &localTri, 
const exa_set<TriFaceVerts> &remote,emInt nDivs){

	//emInt Lvert0 = localTri.getRemoteIndices(0);
	//emInt Lvert1 = localTri.getRemoteIndices(1);
	//emInt Lvert2 = localTri.getRemoteIndices(2);

	emInt global[3]=
	{
		localTri.getGlobalCorner(0), 
		localTri.getGlobalCorner(1), 
		localTri.getGlobalCorner(2)
	};

	TriFaceVerts TF(nDivs,global); 
	TF.setCompare(true); 
	//TriFaceVerts TF(nDivs,Lvert0,Lvert1,Lvert2);

	auto iterTris=remote.find(TF); 
	assert(iterTris!=remote.end()); 
	emInt vert0= localTri.getGlobalCorner(0); 
	emInt vert1= localTri.getGlobalCorner(1); 
	emInt vert2= localTri.getGlobalCorner(2); 
	// emInt globalIterTris [3]={
	// 	iterTris->getGlobalCorner(0),
	// 	iterTris->getGlobalCorner(1),
	// 	iterTris->getGlobalCorner(2)

	// };
	int rotCase = 0;
		for (int cc = 0; cc < 3; cc++) {
			if (vert0 == iterTris->getGlobalCorner(cc)) {
				if (vert1 == iterTris->getGlobalCorner((cc+1)%3)){
					assert(vert2 == iterTris->getGlobalCorner((cc+2)%3));
					rotCase = cc+1;
				}
				else {
					assert(vert1 == iterTris->getGlobalCorner((cc+2)%3));
					assert(vert2 == iterTris->getGlobalCorner((cc+1)%3));
					rotCase = -(cc+1);
				}
			}
		}
		assert(rotCase != 0);
		return rotCase; 

}
inline emInt getQuadRotation(const QuadFaceVerts &localQuad, 
const exa_set<QuadFaceVerts> &remote,emInt nDivs){

	// emInt Lvert0 = localQuad.getRemoteIndices(0);
	// emInt Lvert1 = localQuad.getRemoteIndices(1);
	// emInt Lvert2 = localQuad.getRemoteIndices(2);
	// emInt Lvert3 = localQuad.getRemoteIndices(3); 
	emInt global [4]=
	{
		localQuad.getGlobalCorner(0),
		localQuad.getGlobalCorner(1), 
		localQuad.getGlobalCorner(2), 
		localQuad.getGlobalCorner(3)
	}; 


	QuadFaceVerts QF (nDivs,global); 
	QF.setCompare(true); 
	//printQuads(remote,nDivs); 



	//QuadFaceVerts QF(nDivs,Lvert0,Lvert1,Lvert2,Lvert3);

	auto iterQuads=remote.find(QF); 

	assert(iterQuads!=remote.end()); 

	emInt vert0= localQuad.getGlobalCorner(0); 
	emInt vert1= localQuad.getGlobalCorner(1); 
	emInt vert2= localQuad.getGlobalCorner(2); 
	emInt vert3= localQuad.getGlobalCorner(3); 
	// emInt globalIterTris [4]={
	// 	iterQuads->getGlobalCorner(0),
	// 	iterQuads->getGlobalCorner(1),
	// 	iterQuads->getGlobalCorner(2),
	// 	iterQuads->getGlobalCorner(3)

	// };
	int rotCase = 0;
	
	for (int cc = 0; cc < 4; cc++){
		if (vert0 == iterQuads->getGlobalCorner(cc)) {
			if (vert1 == iterQuads->getGlobalCorner((cc+1)%4)) {
				// Oriented forward; bdry quad
				assert(vert2 == iterQuads->getGlobalCorner((cc+2)%4));
				assert(vert3 == iterQuads->getGlobalCorner((cc+3)%4));
				rotCase = cc+1;
			}
			else {
				assert(vert1 == iterQuads->getGlobalCorner((cc+3)%4));
				assert(vert2 == iterQuads->getGlobalCorner((cc+2)%4));
				assert(vert3 == iterQuads->getGlobalCorner((cc+1)%4));
				rotCase = -(cc+1);
			}
		}
	}
	assert(rotCase != 0);

	return rotCase; 

}; 

inline void comPareTri(const TriFaceVerts &localTri, 
const emInt rotation, const emInt nDivs, 
const exa_set<TriFaceVerts> &remoteTriSet
,std::unordered_map<emInt, emInt> &localRemote){
	

	// emInt vert0 = localTri.getRemoteIndices(0);
	// emInt vert1 = localTri.getRemoteIndices(1);
	// emInt vert2 = localTri.getRemoteIndices(2);

	emInt global[3]=
	{
		localTri.getGlobalCorner(0), 
		localTri.getGlobalCorner(1), 
		localTri.getGlobalCorner(2)
	}; 

	TriFaceVerts TF(nDivs,global); 
	TF.setCompare(true); 


	//TriFaceVerts TF(nDivs,vert0,vert1,vert2);

	auto itRemote= remoteTriSet.find(TF); 
	assert(itRemote!=remoteTriSet.end()); 
	assert(localTri.getPartid()==itRemote->getRemotePartid()); 
	assert(localTri.getRemotePartid()==itRemote->getPartid()); 
	// for(auto i=0; i<3; i++){
	// 	assert(localTri.getCorner(i)==
	// 	itRemote->getRemoteIndices(i)); 
	// 	assert(itRemote->getCorner(i)==localTri.getRemoteIndices(i)); 
	// }


	for (int ii = 0; ii <= nDivs ; ii++) {
	 	for (int jj = 0; jj <= nDivs-ii ; jj++) {

			int trueI; 
			int trueJ; 
			itRemote->getTrueIJ(ii,jj,trueI,trueJ,rotation); 
			
			emInt vertLocal=localTri.getIntVertInd(trueI,trueJ); 
			emInt vertRemote=itRemote->getIntVertInd(ii,jj); 
			localRemote.insert({vertLocal,vertRemote}); 
			
			// std::cout<<"ii: "<<ii<<" jj: "<<jj<<"--- ( "<<
			// vertLocal<<" "<<
			// vertRemote<<" )"
			// <<"-> ("<<local->getX(vertLocal)<<", "<<
			// local->getY(vertLocal)<<", "<<local->getZ(vertLocal)<<" ) "<<" ( "<<
			// remote->getX(vertRemote)<<", "<<remote->getY(vertRemote)<<", "<<
			// remote->getZ(vertRemote)<<" ) "
			//   <<std::endl; 
 			
		}
	}
	//std::cout<<std::endl; 



}
inline void comPareQuad(const QuadFaceVerts &localQuad, 
const emInt rotation, const emInt nDivs, 
const exa_set<QuadFaceVerts> &remoteQuadSet,
std::unordered_map<emInt, emInt> &localRemote){
	//exa_set<TriFaceVerts> remoteTriSet= remote->getRefinedPartTris();

	// emInt vert0 = localQuad.getRemoteIndices(0);
	// emInt vert1 = localQuad.getRemoteIndices(1);
	// emInt vert2 = localQuad.getRemoteIndices(2);
	// emInt vert3 = localQuad.getRemoteIndices(3); 

	emInt global [4]=
	{
		localQuad.getGlobalCorner(0), 
		localQuad.getGlobalCorner(1), 
		localQuad.getGlobalCorner(2), 
		localQuad.getGlobalCorner(3)
	}; 
	QuadFaceVerts QF(nDivs,global); 
	QF.setCompare(true); 
	

	//QuadFaceVerts QF(nDivs,vert0,vert1,vert2,vert3);
	

	auto itRemote= remoteQuadSet.find(QF); 
	assert(itRemote!=remoteQuadSet.end()); 
	assert(localQuad.getPartid()==itRemote->getRemotePartid()); 
	assert(localQuad.getRemotePartid()==itRemote->getPartid()); 
	// for(auto i=0; i<4; i++){
	// 	assert(localQuad.getCorner(i)==
	// 	itRemote->getRemoteIndices(i)); 
	// 	assert(itRemote->getCorner(i)==localQuad.getRemoteIndices(i)); 
	// }
	for (int ii = 0; ii <= nDivs ; ii++) 
	{
	 	for (int jj = 0; jj <= nDivs ; jj++) 
		{

			int trueI; 
			int trueJ; 
			itRemote->getTrueIJ(ii,jj,trueI,trueJ,rotation); 
			emInt vertLocal=localQuad.getIntVertInd(ii,jj); 
			emInt vertRemote=itRemote->getIntVertInd(trueI,trueJ); 
			//emInt vertLocal  = localQuad.getIntVertInd(trueI,trueJ); 
			//emInt vertRemote = itRemote->getIntVertInd(ii,jj);
			//std::cout<<vertLocal<<" "<<vertRemote<<std::endl;
			localRemote.insert({vertLocal,vertRemote});  
			
			// std::cout<<"ii: "<<ii<<" jj: "<<jj<<"--- ( "<<
			// vertLocal<<" "<<
			// vertRemote<<" )"
			// <<"-> ("<<local->getX(vertLocal)<<", "<<
			// local->getY(vertLocal)<<", "<<local->getZ(vertLocal)<<" ) "<<" ( "<<
			// remote->getX(vertRemote)<<", "<<remote->getY(vertRemote)<<", "<<
			// remote->getZ(vertRemote)<<" ) "
			//   <<std::endl; 

	 			
		}
	}	
}


inline void printTris(const TriFaceVerts &tris)
{
	std::cout<<"Num divison: "<<tris.getNumDivs()<<" "<<
	// for(auto itr=tris.begin(); itr!=tris.end();itr++){
	"Part: " << tris.getPartid() <<
		// " local indices: "<<itr->getCorner(0)<<" "<<
		// itr->getCorner(1)<<" "<<itr->getCorner(2)<<
		// " global: "<<itr->getGlobalSorted(0)<<" "<<
		// itr->getGlobalSorted(1)<<" "<<itr->getGlobalSorted(2)<<
		" Unsorted global: " << tris.getGlobalCorner(0) << " " << tris.getGlobalCorner(1) << " " << tris.getGlobalCorner(2) << " Remote ID: " << tris.getRemotePartid() << std::endl;
	//<<
	// " Remote Indices: "<<itr->getRemoteIndices(0)<<" "<<
	// itr->getRemoteIndices(1)<<" "<<
	// itr->getRemoteIndices(2)<<
	//" boolean value: "<<itr->getGlobalCompare()<<
	// std::endl;
	// std::cout<<"Refined verts: "<<std::endl;
	// for (int ii = 0; ii <= nDivs ; ii++) {
	// 	for (int jj = 0; jj <= nDivs-ii ; jj++) {

	// 		 std::cout<<itr->getIntVertInd(ii,jj)<<" ";
	// 	}
	// }

	std::cout << std::endl;
	//}
	//}
}
template <typename T>
inline void SetToVector(const std::unordered_set<T>& sourceSet, 
std::vector<T>& destinationVector) {
    destinationVector.clear();
	// Change to the assignment 
   	for (const auto& element : sourceSet) {
        destinationVector.push_back(element);
    }
}

template <typename T>
inline void vectorToSet(const std::vector<T>& sourceVector, 
std::unordered_set<T>& destinationSet) {
    destinationSet.clear();
  
    for (const auto& element : sourceVector) {
        destinationSet.insert(element);
    }
	assert(destinationSet.size()==sourceVector.size()); 
}
using triHash               = std::unordered_set<TriFaceVerts>; 
using quadHash              = std::unordered_set<QuadFaceVerts>; 
using vecTriHash            = std::vector<triHash>; 
using vecQuadHash           = std::vector<quadHash>;
using triVec                = std::vector<TriFaceVerts>; 
using quadVec               = std::vector<QuadFaceVerts>; 
using vecTriVec             = std::vector<triVec>  ; 
using vecQuadVec            = std::vector<quadVec> ; 

// using vecPartPtr               = std::unistd::vector<Part>; 
// using vecCellPartDataPtr       = std::vector<CellPartData>; 

BOOST_SERIALIZATION_ASSUME_ABSTRACT(FaceVerts)

namespace boost { namespace mpi {
	template <>
struct is_mpi_datatype<TriFaceVerts> : mpl::true_ { };
} }
namespace boost { namespace mpi {
	template <>
struct is_mpi_datatype<QuadFaceVerts> : mpl::true_ { };
} }
//using BoostMPI     = boost::mpi; 
#endif /* SRC_EXA_DEFS_H_ */
