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

typedef uint32_t emInt;
#define EMINT_MAX UINT_MAX

class UMesh;

struct MeshSize {
	emInt nBdryVerts, nVerts, nBdryTris, nBdryQuads, nTets, nPyrs, nPrisms,
			nHexes;
};

bool computeMeshSize(const struct MeshSize& MSIn, const emInt nDivs,
		struct MeshSize& MSOut);

emInt subdividePartMesh(const UMesh * const pVM_input, UMesh * const pVM_output,
		const int nDivs);


#endif /* SRC_EXAMESH_H_ */
