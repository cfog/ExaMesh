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
#include <boost/mpi/datatype.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/mpi.hpp> 



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

#include <set> 
#include <map>

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

#define MAXADJ 6
#define MAX_DIVS 50
#define FILE_NAME_LEN 1024
#define TOLTEST 1e-9
#define MASTER 0 
typedef uint32_t emInt;
#define EMINT_MAX UINT_MAX

#if (HAVE_CGNS == 0)
#define CGNS_ENUMV(a) a
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
#else
#include <cgnslib.h>
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


struct RefineStats {
	double refineTime, extractTime;
	emInt cells;
	size_t fileSize;
};
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
struct hashFunctionCell2Cell 
{
    std::size_t operator()(const std::pair<emInt, emInt>& p) const 
	{
        // Combine the hash values of the pair's components
        std::size_t hashValue = std::hash<emInt>()(p.first);
        hashValue ^= std::hash<emInt>()(p.second) + 0x9e3779b9 + (hashValue << 6) + (hashValue >> 2);
        return hashValue;
    }
};



struct hashFunctionFace2Cell 
{
    std::size_t operator()(const std::vector<emInt>& s) const {
        // Implement a custom hash function for std::set<int>
        // You might want to combine the hash values of individual elements in the set
        // to create a hash value for the entire set.
        // Example:
        std::size_t hash_value = 0;
        for (int element : s) {
            hash_value ^= std::hash<emInt>()(element) + 0x9e3779b9 + (hash_value << 6) + (hash_value >> 2);
        }
        return hash_value;
    }
};


struct pairHash {
    template <class T1, class T2>
    std::size_t operator () (std::pair<T1,T2> const& pair) const {
        auto first = std::min(pair.first, pair.second);
        auto second = std::max(pair.first, pair.second);

        std::size_t h1 = std::hash<T1>{}(first);
        std::size_t h2 = std::hash<T2>{}(second);

        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, use sth. like boost.hash_combine
        return h1 ^ h2;  
    }
};

using vecSizes                   = std::vector<std::size_t>;
using multimpFace2Cell			 = std::unordered_multimap < std::vector<emInt>, std::pair<emInt,emInt>, hashFunctionFace2Cell>;
using TableCell2Cell			 = std::unordered_map < std::pair<emInt,emInt>, std::set<emInt>, hashFunctionCell2Cell>;
using VecVecReqs                 = std::vector<std::vector<boost::mpi::request>>; 
using vecPairInt2Int			 = std::vector<std::pair<int, int>>;
using Request                    = boost::mpi::request;
using vecReqs                    = std::vector<Request>;
using vecVecReqs                 = std::vector<vecReqs>;
		
struct timeResults
{
	
	double preProcessing;
	double partition;
	double InitialFaceMatching; 
	double broadcasting; 
	double serial; 
	double extract; 
	double refine;
	double faceExchange; 
	double matchtris;
	double matchquads; 
	double totalMatch; 
	double waitTri; 
	double waitQuad;
	double totalFacesWait; 
	double total; 

};

inline 
FILE*
openFile (std::string fileName)
{
	FILE *file = fopen(fileName.c_str(), "a");
    if (file == NULL) 
	{
        fprintf(stderr, "Error opening the file!\n");
    }
    return file; 
}
inline 
long 
getSeekFile (FILE *file)
{
	fseek(file, 0, SEEK_END);
	long size = ftell(file);
	return size; 
}


#endif /* SRC_EXA_DEFS_H_ */
