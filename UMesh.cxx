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
 * UMesh.cxx
 *
 *  Created on: Jul. 2, 2019
 *      Author: cfog
 */

#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <memory>
#include <set>
#include <vector>
#include <string.h>
#include <map>
#include <fstream>
#include <chrono>
#include <thread>

// This include file is deliberately before ExaMesh headers so
// there aren't warnings about standard autoconf things being
// redefined.
#include "GMGW_FileWrapper.hxx"
#undef PACKAGE_TARNAME
#undef PACKAGE_NAME
#undef PACKAGE_VERSION
#undef PACKAGE_STRING
#undef PACKAGE_BUGREPORT
#undef PACKAGE_URL

#include <mpi.h>

#include "BdryTriDivider.h"
#include "ExaMesh.h"
#include "exa-defs.h"
#include "UMesh.h"
#include <cstdint>
#include <cstddef>
#include "resultGenerator.cxx"
#if (HAVE_CGNS == 1)
#include <cgnslib.h>
#endif

using std::cout;
using std::endl;

#ifndef NDEBUG

static bool memoryCheck(void *address, int nBytes) {
	char *checkPtr = reinterpret_cast<char*>(address);
	bool retVal = true;
	for (int ii = 0; ii < nBytes; ii++) {
		retVal = retVal && (checkPtr[ii] == 0x00);
	}
	return retVal;
}
#endif

void UMesh::init(const emInt nVerts, const emInt nBdryVerts,
		const emInt nBdryTris, const emInt nBdryQuads, const emInt nTets,
		const emInt nPyramids, const emInt nPrisms, const emInt nHexes) {
	m_nVerts = nVerts;
	m_nBdryVerts = nBdryVerts;
	m_nTris = nBdryTris;
	m_nQuads = nBdryQuads;
	m_nTets = nTets;
	m_nPyrs = nPyramids;
	m_nPrisms = nPrisms;
	m_nHexes = nHexes;

	// All sizes are computed in bytes.
	// Work out buffer size, including padding to ensure 8-byte alignment for the coordinates.
	size_t intSize = sizeof(emInt);
	size_t headerSize = 7 * intSize;
	// How many bytes to add to get eight-byte alignment for coordinates,
	// assuming eight-byte alignment for the buffer overall.
	size_t slack1Size = (intSize == 4) ? 4 : 0;
	size_t coordSize = 3 * sizeof(double) * nVerts;
	size_t connSize = (3 * nBdryTris + 4 * nBdryQuads + 4 * size_t(nTets)
			+ 5 * size_t(nPyramids) + 6 * size_t(nPrisms) + 8 * size_t(nHexes))
			* intSize;
	size_t BCSize = (nBdryTris + nBdryQuads) * intSize;

	// How many bytes to add to fill up the last eight-byte chunk?
	size_t slack2Size = ((((connSize + BCSize) / 8 + 1) * 8)
			- (connSize + BCSize)) % 8;
	size_t bufferBytes = headerSize + coordSize + connSize + BCSize + slack1Size
			+ slack2Size;
	assert((headerSize + slack1Size) % 8 == 0);
	assert((connSize + BCSize + slack2Size) % 8 == 0);
	assert(bufferBytes % 8 == 0);
	size_t bufferWords = bufferBytes / 8;
	// Use words instead of bytes to ensure 8-byte alignment.
	m_buffer = reinterpret_cast<char*>(calloc(bufferWords, 8));
	assert(m_buffer!=NULL);

	// The pointer arithmetic here is made more complicated because the pointers aren't
	// compatible with each other.
	m_header = reinterpret_cast<emInt*>(m_buffer + slack1Size);
	assert(m_header!=NULL);
	std::fill(m_header, m_header + 7, 0);
	m_coords = reinterpret_cast<double (*)[3]>(m_buffer + headerSize
			+ slack1Size);
	m_TriConn = reinterpret_cast<emInt (*)[3]>(m_buffer + headerSize
			+ slack1Size + coordSize);
	m_QuadConn = reinterpret_cast<emInt (*)[4]>(m_TriConn + nBdryTris);
	m_TriBC = reinterpret_cast<emInt*>(m_QuadConn + nBdryQuads);
	m_QuadBC = m_TriBC + nBdryTris;
	m_TetConn = reinterpret_cast<emInt (*)[4]>(m_QuadBC + nBdryQuads);
	m_PyrConn = reinterpret_cast<emInt (*)[5]>(m_TetConn + nTets);
	m_PrismConn = reinterpret_cast<emInt (*)[6]>(m_PyrConn + nPyramids);
	m_HexConn = reinterpret_cast<emInt (*)[8]>(m_PrismConn + nPrisms);
	m_fileImage = m_buffer + slack1Size;
	m_fileImageSize = bufferBytes - slack1Size - slack2Size;

	//	printf("Diagnostics for UMesh data struct:\n");
	//	printf("Buffer size, in bytes:     %lu\n", bufferBytes);
	//	printf("File image size, in bytes: %lu\n", m_fileImageSize);
	//	printf("Num verts: %10u\n", nVerts);
	//	printf("Num tris:  %10u\n", nBdryTris);
	//	printf("Num quads: %10u\n", nBdryQuads);
	//	printf("Num tets:  %10u\n", nTets);
	//	printf("Num pyrs:  %10u\n", nPyramids);
	//	printf("Num prisms: %8u\n", nPrisms);
	//	printf("Num hexes: %10u\n", nHexes);
	//	printf(
	//			"Coord offset: %10lu\n",
	//			reinterpret_cast<char*>(m_coords) - reinterpret_cast<char*>(m_header));
	//	printf(
	//			"Tri conn offset: %10lu\n",
	//			reinterpret_cast<char*>(m_TriConn) - reinterpret_cast<char*>(m_header));
	//	printf("Tri BC offset: %10lu\n",
	//					reinterpret_cast<char*>(m_TriBC) - reinterpret_cast<char*>(m_header));
	//	printf(
	//			"Tet conn offset: %10lu\n",
	//			reinterpret_cast<char*>(m_TetConn) - reinterpret_cast<char*>(m_header));

	m_lenScale = new double[m_nVerts];
	for (emInt ii = 0; ii < m_nVerts; ii++) {
		m_lenScale[ii] = 0;
	}
}

UMesh::UMesh(const emInt nVerts, const emInt nBdryVerts, const emInt nBdryTris,
		const emInt nBdryQuads, const emInt nTets, const emInt nPyramids,
		const emInt nPrisms, const emInt nHexes) :
		m_nVerts(nVerts), m_nBdryVerts(nBdryVerts), m_nTris(nBdryTris), m_nQuads(
				nBdryQuads), m_nTets(nTets), m_nPyrs(nPyramids), m_nPrisms(
				nPrisms), m_nHexes(nHexes), m_fileImageSize(0), m_header(
				nullptr), m_coords(nullptr), m_TriConn(nullptr), m_QuadConn(
				nullptr), m_TetConn(nullptr), m_PyrConn(nullptr), m_PrismConn(
				nullptr), m_HexConn(nullptr), m_buffer(nullptr), m_fileImage(
				nullptr) {

	// All sizes are computed in bytes.

	// Work out buffer size, including padding to ensure 8-byte alignment for the coordinates.
	init(nVerts, nBdryVerts, nBdryTris, nBdryQuads, nTets, nPyramids, nPrisms,
			nHexes);
}

UMesh::UMesh(const emInt nVerts, const emInt nBdryVerts, const emInt nBdryTris,
		const emInt nBdryQuads, const emInt nTets, const emInt nPyramids,
		const emInt nPrisms, const emInt nHexes,
		const std::vector<std::vector<emInt>> &tetConns,
		const std::vector<emInt> &header,
		const std::vector<std::pair<emInt, emInt>> &cellID2cellTypeLocal,
		const std::vector<std::vector<emInt>> &vTriConns,
		const std::vector<double> &vLengthSclae,
		const std::vector<std::vector<emInt>> &vQuadConn,
		const std::vector<std::vector<emInt>> &vPyrConn,
		const std::vector<std::vector<emInt>> &vPrsimConn,
		const std::vector<std::vector<emInt>> &vHexConn) :
		m_nVerts(nVerts), m_nBdryVerts(nBdryVerts), m_nTris(nBdryTris), m_nQuads(
				nBdryQuads), m_nTets(nTets), m_nPyrs(nPyramids), m_nPrisms(
				nPrisms), m_nHexes(nHexes), m_fileImageSize(0), m_header(
				nullptr), m_coords(nullptr), m_TriConn(nullptr), m_QuadConn(
				nullptr), m_TetConn(nullptr), m_PyrConn(nullptr), m_PrismConn(
				nullptr), m_HexConn(nullptr), m_buffer(nullptr), m_fileImage(
				nullptr) {

	// All sizes are computed in bytes.

	// Work out buffer size, including padding to ensure 8-byte alignment for the coordinates.
	init(nVerts, nBdryVerts, nBdryTris, nBdryQuads, nTets, nPyramids, nPrisms,
			nHexes);

	for (auto i = 0; i < 7; i++) {
		m_header[i] = header[i];
	}

	for (emInt iTet = 0; iTet < nTets; iTet++) {
		for (auto iConn = 0; iConn < 4; iConn++) {
			m_TetConn[iTet][iConn] = tetConns[iTet][iConn];
		}
	}

	cellID2cellTypeLocalID = cellID2cellTypeLocal;

	for (emInt iTri = 0; iTri < header[eTri]; iTri++) {
		for (auto j = 0; j < 3; j++) {
			m_TriConn[iTri][j] = vTriConns[iTri][j];
		}
	}

	for (emInt iVert = 0; iVert < header[eVert]; iVert++) {
		m_lenScale[iVert] = vLengthSclae[iVert];
	}
	for (emInt iQuad = 0; iQuad < header[eQuad]; iQuad++) {
		for (auto j = 0; j < 4; j++) {
			m_QuadConn[iQuad][j] = vQuadConn[iQuad][j];
		}
	}
	for (emInt iPyr = 0; iPyr < header[ePyr]; iPyr++) {
		for (auto j = 0; j < 5; j++) {
			m_PyrConn[iPyr][j] = vPyrConn[iPyr][j];
		}
	}
	for (emInt iPrism = 0; iPrism < header[ePrism]; iPrism++) {
		for (auto j = 0; j < 6; j++) {
			m_PrismConn[iPrism][j] = vPrsimConn[iPrism][j];
		}
	}
	for (emInt iHex = 0; iHex < header[eHex]; iHex++) {
		for (auto j = 0; j < 8; j++) {
			m_HexConn[iHex][j] = vHexConn[iHex][j];
		}
	}

}

emInt UMesh::addVert(const double newCoords[3]) {
	assert(m_header[eVert] < m_nVerts);
	assert(memoryCheck(m_coords[m_header[eVert]], 24));
	m_coords[m_header[eVert]][0] = newCoords[0];
	m_coords[m_header[eVert]][1] = newCoords[1];
	m_coords[m_header[eVert]][2] = newCoords[2];
	return (m_header[eVert]++);
}

emInt UMesh::addBdryTri(const emInt verts[3]) {
	assert(memoryCheck(m_TriConn[m_header[eTri]], 3 * sizeof(emInt)));
	for (int ii = 0; ii < 3; ii++) {
		assert(verts[ii] < m_header[eVert]);
		m_TriConn[m_header[eTri]][ii] = verts[ii];
	}
	//vTriConns.push_back({verts[0],verts[1],verts[2]}); 
	return (m_header[eTri]++);
}

emInt UMesh::addBdryQuad(const emInt verts[4]) {
	assert(memoryCheck(m_QuadConn[m_header[eQuad]], 4 * sizeof(emInt)));
	for (int ii = 0; ii < 4; ii++) {
		assert(verts[ii] < m_header[eVert]);
		m_QuadConn[m_header[eQuad]][ii] = verts[ii];
	}
	//vQuadConns.push_back({verts[0],verts[1],verts[2],verts[3]});
	return (m_header[eQuad]++);
}

emInt UMesh::addTet(const emInt verts[4]) {
#ifndef OLD_ADD_ELEMENT
	emInt thisTetInd = m_header[eTet]++;
	emInt *thisConn = m_TetConn[thisTetInd];
	assert(memoryCheck(thisConn, 4 * sizeof(emInt)));
	std::copy(verts, verts + 4, thisConn);
	//vTetConns.push_back({verts[0],verts[1],verts[2],verts[3]});
	return thisTetInd;
#else
	assert(memoryCheck(m_TetConn[m_header[eTet]], 4 * sizeof(emInt)));
	for (int ii = 0; ii < 4; ii++)
	{
		assert(verts[ii] < m_header[eVert]);
		m_TetConn[m_header[eTet]][ii] = verts[ii];
	}
	vTetConns.push_back({verts[0],verts[1],verts[2],verts[3]});
	return (m_header[eTet]++);
#endif
}

emInt UMesh::addPyramid(const emInt verts[5]) {
#ifndef OLD_ADD_ELEMENT
	emInt thisPyrInd = m_header[ePyr]++;
	emInt *thisConn = m_PyrConn[thisPyrInd];
	assert(memoryCheck(thisConn, 5 * sizeof(emInt)));
	std::copy(verts, verts + 5, thisConn);
	//vPyrmConns.push_back({verts[0],verts[1],verts[2],verts[3],verts[4]});
	return thisPyrInd;
#else
	assert(memoryCheck(m_PyrConn[m_header[ePyr]], 5 * sizeof(emInt)));
	for (int ii = 0; ii < 5; ii++)
	{
		assert(verts[ii] < m_header[eVert]);
		m_PyrConn[m_header[ePyr]][ii] = verts[ii];
	}
	vPyrmConns.push_back({verts[0],verts[1],verts[2],verts[3],verts[4]});
	return (m_header[ePyr]++);
#endif
}

emInt UMesh::addPrism(const emInt verts[6]) {
#ifndef OLD_ADD_ELEMENT
	emInt thisPrismInd = m_header[ePrism]++;
	emInt *thisConn = m_PrismConn[thisPrismInd];
	assert(memoryCheck(thisConn, 6 * sizeof(emInt)));
	std::copy(verts, verts + 6, thisConn);
	//vPrsimConns.push_back({verts[0],verts[1],verts[2],verts[3],verts[4],verts[5]});
	return thisPrismInd;
#else
	assert(memoryCheck(m_PrismConn[m_header[ePrism]], 6 * sizeof(emInt)));
	for (int ii = 0; ii < 6; ii++)
	{
		assert(verts[ii] < m_header[eVert]);
		m_PrismConn[m_header[ePrism]][ii] = verts[ii];
	}
	vPrsimConns.push_back({verts[0],verts[1],verts[2],verts[3],verts[4],verts[5]});
	return (m_header[ePrism]++);
#endif
}
#define OLD_ADD_ELEMENT
emInt UMesh::addHex(const emInt verts[8]) {
#ifndef OLD_ADD_ELEMENT
	emInt thisHexInd = m_header[eHex]++;
	emInt *thisConn = m_HexConn[thisHexInd];
	assert(memoryCheck(thisConn, 8 * sizeof(emInt)));
	std::copy(verts, verts + 8, thisConn);
	vHexConns.push_back({verts[0],verts[1],verts[2],verts[3],verts[4],verts[5],verts[6],verts[7]}); 
	return thisHexInd;
#else
	assert(memoryCheck(m_HexConn[m_header[eHex]], 8 * sizeof(emInt)));
	for (int ii = 0; ii < 8; ii++) {
		assert(verts[ii] < m_header[eVert]);
		m_HexConn[m_header[eHex]][ii] = verts[ii];
	}
	//vHexConns.push_back({verts[0],verts[1],verts[2],verts[3],verts[4],verts[5],verts[6],verts[7]}); 
	return (m_header[eHex]++);
#endif
}

UMesh::~UMesh() {
	free(m_buffer);
}

void checkConnectivitySize(const char cellType, const emInt nVerts) {
	emInt expected[] = { 0, 0, 0, 0, 0, 3, 0, 4, 0, 0, 4, 0, 5, 0, 6, 0, 0, 8 };
	if (expected[int(cellType)] != nVerts) {
		fprintf(
		stderr,
				"Error reading mesh file.  Cell type %d expects %u verts; found %u.\n",
				cellType, expected[int(cellType)], nVerts);
		exit(1);
	}
}

class vertTriple {
	emInt corners[3];

public:
	vertTriple(const emInt vA, const emInt vB, const emInt vC) {
		corners[0] = vA;
		corners[1] = vB;
		corners[2] = vC;
	}
	bool operator<(const vertTriple &that) const {
		emInt thisTemp[3], thatTemp[3];
		sortVerts3(corners, thisTemp);
		sortVerts3(that.corners, thatTemp);
		return (thisTemp[0] < thatTemp[0]
				|| (thisTemp[0] == thatTemp[0] && thisTemp[1] < thatTemp[1])
				|| (thisTemp[0] == thatTemp[0] && thisTemp[1] == thatTemp[1]
						&& thisTemp[2] < thatTemp[2]));
	}
	const emInt* getCorners() const {
		return corners;
	}
};

class vertQuadruple {
	emInt corners[4];

public:
	vertQuadruple(const emInt vA, const emInt vB, const emInt vC,
			const emInt vD) {
		corners[0] = vA;
		corners[1] = vB;
		corners[2] = vC;
		corners[3] = vD;
	}
	bool operator<(const vertQuadruple &that) const {
		emInt thisTemp[4], thatTemp[4];
		sortVerts4(corners, thisTemp);
		sortVerts4(that.corners, thatTemp);
		return (thisTemp[0] < thatTemp[0]
				|| (thisTemp[0] == thatTemp[0] && thisTemp[1] < thatTemp[1])
				|| (thisTemp[0] == thatTemp[0] && thisTemp[1] == thatTemp[1]
						&& thisTemp[2] < thatTemp[2])
				|| (thisTemp[0] == thatTemp[0] && thisTemp[1] == thatTemp[1]
						&& thisTemp[2] == thatTemp[2]
						&& thisTemp[3] < thatTemp[3]));
	}
	const emInt* getCorners() const {
		return corners;
	}
};

void updateTriSet(std::set<vertTriple> &triSet, const emInt v0, const emInt v1,
		const emInt v2) {
	vertTriple VT(v0, v1, v2);
	typename std::set<vertTriple>::iterator vertIter, VIend = triSet.end();

	//vertIter = triSet.find(VT);
	//if (vertIter == VIend)
	//{
	auto inserResult = triSet.insert(VT);
	if (!inserResult.second) {
		triSet.erase(inserResult.first);
	}
	//}
	//else
	//{
	//triSet.erase(vertIter);
	//}
}

void updateQuadSet(std::set<vertQuadruple> &quadSet, const emInt v0,
		const emInt v1, const emInt v2, const emInt v3) {
	vertQuadruple VQ(v0, v1, v2, v3);
	typename std::set<vertQuadruple>::iterator vertIter, VQend = quadSet.end();

	//vertIter = quadSet.find(VQ);
	//if (vertIter == VQend)
	//{
	auto inserResult = quadSet.insert(VQ);
	if (!inserResult.second) {
		quadSet.erase(inserResult.first);
	}
	//}
	//else
	//{
	// If the elmenet is already in the set, remove it.
	//quadSet.erase(vertIter);
	//}
}

UMesh::UMesh(const std::string baseFileName, const std::string fileSuffix,
		const std::string ugridInfix) :
		m_nVerts(0), m_nBdryVerts(0), m_nTris(0), m_nQuads(0), m_nTets(0), m_nPyrs(
				0), m_nPrisms(0), m_nHexes(0), m_fileImageSize(0), m_header(
				nullptr), m_coords(nullptr), m_TriConn(nullptr), m_QuadConn(
				nullptr), m_TetConn(nullptr), m_PyrConn(nullptr), m_PrismConn(
				nullptr), m_HexConn(nullptr), m_buffer(nullptr), m_fileImage(
				nullptr) {
	// Use the same IO routines as the mesh analyzer code from GMGW.
	FileWrapper *reader = FileWrapper::factory(baseFileName.c_str(),
			fileSuffix.c_str(), ugridInfix.c_str());

	reader->scanFile();

	// Identify any bdry tris and quads that aren't in the file

	emInt nBdryTris = reader->getNumBdryTris();
	emInt nBdryQuads = reader->getNumBdryQuads();

	std::set<vertTriple> setTris;
	std::set<vertQuadruple> setQuads;

	reader->seekStartOfConnectivity();

	vcellID2type.resize(reader->getNumCells());

	for (emInt ii = 0; ii < reader->getNumCells(); ii++) {
		char cellType = reader->getCellType(ii);
		vcellID2type[ii] = cellType;
		unsigned nConn, connect[8];
		reader->getNextCellConnectivity(nConn, connect);
		checkConnectivitySize(cellType, nConn);

		switch (cellType) {
		case CGNS_ENUMV(TRI_3):
			updateTriSet(setTris, connect[0], connect[1], connect[2]);

			break;
		case CGNS_ENUMV(QUAD_4):
			updateQuadSet(setQuads, connect[0], connect[1], connect[2],
					connect[3]);

			break;
		case CGNS_ENUMV(TETRA_4):
			updateTriSet(setTris, connect[0], connect[1], connect[2]);
			updateTriSet(setTris, connect[0], connect[1], connect[3]);
			updateTriSet(setTris, connect[1], connect[2], connect[3]);
			updateTriSet(setTris, connect[2], connect[0], connect[3]);
			break;
		case CGNS_ENUMV(PYRA_5):
			updateTriSet(setTris, connect[0], connect[1], connect[4]);
			updateTriSet(setTris, connect[1], connect[2], connect[4]);
			updateTriSet(setTris, connect[2], connect[3], connect[4]);
			updateTriSet(setTris, connect[3], connect[0], connect[4]);
			updateQuadSet(setQuads, connect[0], connect[1], connect[2],
					connect[3]);
			break;
		case CGNS_ENUMV(PENTA_6):
			updateTriSet(setTris, connect[0], connect[1], connect[2]);
			updateTriSet(setTris, connect[3], connect[4], connect[5]);
			updateQuadSet(setQuads, connect[0], connect[1], connect[4],
					connect[3]);
			updateQuadSet(setQuads, connect[1], connect[2], connect[5],
					connect[4]);
			updateQuadSet(setQuads, connect[2], connect[0], connect[3],
					connect[5]);
			break;
		case CGNS_ENUMV(HEXA_8):
			updateQuadSet(setQuads, connect[0], connect[1], connect[2],
					connect[3]);
			updateQuadSet(setQuads, connect[4], connect[5], connect[6],
					connect[7]);
			updateQuadSet(setQuads, connect[0], connect[1], connect[5],
					connect[4]);
			updateQuadSet(setQuads, connect[1], connect[2], connect[6],
					connect[5]);
			updateQuadSet(setQuads, connect[2], connect[3], connect[7],
					connect[6]);
			updateQuadSet(setQuads, connect[3], connect[0], connect[4],
					connect[7]);
			break;
		default:
			assert(0);
		}
	}
#ifndef NDEBUG
	//testCell2FaceConn(reader->getNumCells()); 
#endif
	//double start = exaTime(); 
	//double time  = exaTime()-start; 

	//std::cout<<"Building cell2cell connectivity took: "<<time<<std::endl; 

	nBdryTris += setTris.size();
	nBdryQuads += setQuads.size();

	//BdryVertsTrisQuads.push_back(reader->getNumBdryVerts());
	//BdryVertsTrisQuads.push_back(numBdryTris); 
	//BdryVertsTrisQuads.push_back(numBdryQuads); 

	init(reader->getNumVerts(), reader->getNumBdryVerts(), nBdryTris,
			nBdryQuads, reader->getNumTets(), reader->getNumPyramids(),
			reader->getNumPrisms(), reader->getNumHexes());

	reader->seekStartOfCoords();
	for (emInt ii = 0; ii < m_nVerts; ii++) {
		double coords[3];
		reader->getNextVertexCoords(coords[0], coords[1], coords[2]);
		addVert(coords);
	}

	reader->seekStartOfConnectivity();

	// Initializing counters to be zero 
	// Transforming enum 
	for (emInt ii = 0; ii < reader->getNumCells(); ii++) {

		char cellType = reader->getCellType(ii);
		emInt nConn, connect[8];
		reader->getNextCellConnectivity(nConn, connect);
		checkConnectivitySize(cellType, nConn);
		switch (cellType) {
		case CGNS_ENUMV(TRI_3):
			addBdryTri(connect);
			break;
		case CGNS_ENUMV(QUAD_4):
			addBdryQuad(connect);
			break;
		case CGNS_ENUMV(TETRA_4):
			addTet(connect);
			break;
		case CGNS_ENUMV(PYRA_5):
			addPyramid(connect);
			break;
		case CGNS_ENUMV(PENTA_6):
			addPrism(connect);
			break;
		case CGNS_ENUMV(HEXA_8):
			addHex(connect);
			break;
		default:
			assert(0);
		}
	}

	for (auto VT : setTris) {
		const emInt *const corners = VT.getCorners();
		addBdryTri(corners);
	}
	for (auto VQ : setQuads) {
		const emInt *const corners = VQ.getCorners();
		addBdryQuad(corners);
	}

	// Now tag all bdry verts
	bool *isBdryVert = new bool[m_nVerts];
	for (emInt ii = 0; ii < m_nVerts; ii++) {
		isBdryVert[ii] = false;
	}
	for (emInt iTri = 0; iTri < m_nTris; iTri++) {
		isBdryVert[m_TriConn[iTri][0]] = true;
		isBdryVert[m_TriConn[iTri][1]] = true;
		isBdryVert[m_TriConn[iTri][2]] = true;
	}
	for (emInt iQuad = 0; iQuad < m_nQuads; iQuad++) {
		isBdryVert[m_QuadConn[iQuad][0]] = true;
		isBdryVert[m_QuadConn[iQuad][1]] = true;
		isBdryVert[m_QuadConn[iQuad][2]] = true;
		isBdryVert[m_QuadConn[iQuad][3]] = true;
	}
	m_nBdryVerts = 0;
	for (emInt ii = 0; ii < m_nVerts; ii++) {
		if (isBdryVert[ii]) {
			m_nBdryVerts++;
		}
	}
	delete[] isBdryVert;

	buildCellToCellConnectivity();

	// If any of these fail, your file was invalid.
	assert(m_nVerts == m_header[eVert]);
	assert(m_nTris == m_header[eTri]);
	assert(m_nQuads == m_header[eQuad]);
	assert(m_nTets == m_header[eTet]);
	assert(m_nPyrs == m_header[ePyr]);
	assert(m_nPrisms == m_header[ePrism]);
	assert(m_nHexes == m_header[eHex]);

	//vheader.assign(m_header, m_header+7);

	delete reader;

	setupLengthScales();

	//vLengthScale.assign(m_lenScale,m_lenScale+m_header[0]);
}

std::unique_ptr<UMesh> UMesh::subdivideMesh(const emInt nDivs,
		const emInt partID) const {
	setlocale(LC_ALL, "");
	size_t totalInputCells = size_t(m_nTets) + m_nPyrs + m_nPrisms + m_nHexes;
	fprintf(
	stderr,
			"Initial mesh has:\n %'15u verts,\n %'15u bdry tris,\n %'15u bdry quads,\n %'15u tets,\n %'15u pyramids,\n %'15u prisms,\n %'15u hexes,\n%'15lu cells total\n",
			m_nVerts, m_nTris, m_nQuads, m_nTets, m_nPyrs, m_nPrisms, m_nHexes,
			totalInputCells);

	MeshSize MSOut = computeFineMeshSize(nDivs);
	auto outMesh = std::make_unique<UMesh>(MSOut.nVerts, MSOut.nBdryVerts,
			MSOut.nBdryTris, MSOut.nBdryQuads, MSOut.nTets, MSOut.nPyrs,
			MSOut.nPrisms, MSOut.nHexes);
	// Copy length scale data from the other mesh.
	for (emInt vv = 0; vv < m_nVerts; vv++) {
		outMesh->setLengthScale(vv, m_lenScale[vv]);
	}

	subdividePartMesh(this, outMesh.get(), nDivs, partID);
#ifndef NDEBUG	
	setlocale(LC_ALL, "");
	fprintf(
	stderr,
			"Final mesh has:\n %'15u verts,\n %'15u bdry tris,\n %'15u bdry quads,\n %'15u tets,\n %'15u pyramids,\n %'15u prisms,\n %'15u hexes,\n%'15u cells total\n",
			outMesh->numVerts(), outMesh->numBdryTris(),
			outMesh->numBdryQuads(), outMesh->numTets(), outMesh->numPyramids(),
			outMesh->numPrisms(), outMesh->numHexes(), outMesh->numCells());
#endif		
	return outMesh;
}

bool UMesh::writeVTKFile(const char fileName[]) {
	double timeBefore = exaTime();

	FILE *outFile = fopen(fileName, "w");
	if (!outFile) {
		fprintf(stderr, "Couldn't open file %s for writing.  Bummer!\n",
				fileName);
		return false;
	}

	fprintf(outFile, "# vtk DataFile Version 1.0\n");
	fprintf(outFile, "GRUMMP Tetra example\n");
	fprintf(outFile, "ASCII\n");
	fprintf(outFile, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(outFile, "POINTS %d float\n", m_header[eVert]);

	//-------------------------------------
	// write 3d vertex data
	//-------------------------------------
	for (emInt i = 0; i < m_header[eVert]; ++i) {
		fprintf(outFile, "%16.8f %16.8f %16.8f\n", getX(i), getY(i), getZ(i));
	}

	const size_t nTris = numBdryTris();
	const size_t nQuads = numBdryQuads();
	const size_t nTets = numTets();
	const size_t nPyrs = numPyramids();
	const size_t nPrisms = numPrisms();
	const size_t nHexes = numHexes();

	const size_t numEnts = nTris + nQuads + nTets + nPyrs + nPrisms + nHexes;
	const size_t dataSize = 4 * nTris + 5 * (nQuads + nTets) + 6 * nPyrs
			+ 7 * nPrisms + 9 * nHexes;

	fprintf(outFile, "CELLS %lu %lu\n", numEnts, dataSize);

	// Write all the bdry tris
	for (std::size_t i = 0; i < nTris; i++) {
		const emInt *verts = getBdryTriConn(i);
		fprintf(outFile, "3 %d %d %d\n", verts[0], verts[1], verts[2]);
	}

	// Write all the bdry quads
	for (std::size_t i = 0; i < nQuads; i++) {
		const emInt *verts = getBdryQuadConn(i);
		fprintf(outFile, "4 %d %d %d %d\n", verts[0], verts[1], verts[2],
				verts[3]);
	}

	// Write all the tets
	for (std::size_t i = 0; i < nTets; i++) {
		const emInt *verts = getTetConn(i);
		fprintf(outFile, "4 %d %d %d %d\n", verts[0], verts[1], verts[2],
				verts[3]);
	}

	// Write all the pyramids
	for (std::size_t i = 0; i < nPyrs; i++) {
		const emInt *verts = getPyrConn(i);
		fprintf(outFile, "5 %d %d %d %d %d\n", verts[0], verts[1], verts[2],
				verts[3], verts[4]);
	}

	// Write all the prisms
	for (std::size_t i = 0; i < nPrisms; i++) {
		const emInt *verts = getPrismConn(i);
		fprintf(outFile, "6 %d %d %d %d %d %d\n", verts[0], verts[1], verts[2],
				verts[3], verts[4], verts[5]);
	}

	// Write all the hexes
	for (std::size_t i = 0; i < nHexes; i++) {
		const emInt *verts = getHexConn(i);
		fprintf(outFile, "8 %d %d %d %d %d %d %d %d\n", verts[0], verts[1],
				verts[2], verts[3], verts[4], verts[5], verts[6], verts[7]);
	}

	//-------------------------------------
	// write cell type (VTK_TRIANGLE = 5, VTK_QUAD = 9,
	// VTK_TETRA = 10, VTK_HEXAHEDRON = 12, VTK_WEDGE = 13,
	// VTK_PYRAMID = 14)
	//-------------------------------------
	fprintf(outFile, "CELL_TYPES %lu\n", numEnts);
	for (std::size_t ct = 0; ct < nTris; ++ct)
		fprintf(outFile, "5\n");
	for (std::size_t ct = 0; ct < nQuads; ++ct)
		fprintf(outFile, "9\n");
	for (std::size_t ct = 0; ct < nTets; ++ct)
		fprintf(outFile, "10\n");
	for (std::size_t ct = 0; ct < nPyrs; ++ct)
		fprintf(outFile, "14\n");
	for (std::size_t ct = 0; ct < nPrisms; ++ct)
		fprintf(outFile, "13\n");
	for (std::size_t ct = 0; ct < nHexes; ++ct)
		fprintf(outFile, "12\n");

	fclose(outFile);
	double timeAfter = exaTime();
	double elapsed = timeAfter - timeBefore;
	size_t totalCells = size_t(m_nTets) + m_nPyrs + m_nPrisms + m_nHexes;
	fprintf(stderr, "CPU time for VTK file write = %5.2F seconds\n", elapsed);
	fprintf(stderr, "                          %5.2F million cells / minute\n",
			(totalCells / 1000000.) / (elapsed / 60));

	return true;
}

void UMesh::incrementVertIndices(emInt *conn, emInt size, int inc) {
	for (emInt ii = 0; ii < size; ii++) {
		conn[ii] += inc;
	}
}

bool UMesh::writeUGridFile(const char fileName[]) {
	double timeBefore = exaTime();

	// Need to increment all vert indices, because UGRID files are 1-based.
	emInt size = m_nTris * 3 + m_nQuads * 4;
	incrementVertIndices(reinterpret_cast<emInt*>(m_TriConn), size, 1);
	size = m_nTets * 4 + m_nPyrs * 5 + m_nPrisms * 6 + m_nHexes * 8;
	incrementVertIndices(reinterpret_cast<emInt*>(m_TetConn), size, 1);

	// Also need to swap verts 2 and 4 for pyramids, because UGRID treats
	// pyramids as prisms with the edge from 2 to 5 collapsed.  Compared
	// with the ordering the rest of the world uses, this has the effect
	// of switching verts 2 and 4.
	for (emInt ii = 0; ii < m_nPyrs; ii++) {
		std::swap(m_PyrConn[ii][2], m_PyrConn[ii][4]);
	}

	FILE *outFile = fopen(fileName, "w");
	if (!outFile) {
		fprintf(stderr, "Couldn't open file %s for writing.  Bummer!\n",
				fileName);
		return false;
	}

	fwrite(m_fileImage, m_fileImageSize, 1, outFile);
	fclose(outFile);

	// Need to undo the increment for future use
	size = m_nTris * 3 + m_nQuads * 4;
	incrementVertIndices(reinterpret_cast<emInt*>(m_TriConn), size, -1);
	size = m_nTets * 4 + m_nPyrs * 5 + m_nPrisms * 6 + m_nHexes * 8;
	incrementVertIndices(reinterpret_cast<emInt*>(m_TetConn), size, -1);

	// Need to swap verts 2 and 4 back for future use.
	for (emInt ii = 0; ii < m_nPyrs; ii++) {
		std::swap(m_PyrConn[ii][2], m_PyrConn[ii][4]);
	}

	double timeAfter = exaTime();
	double elapsed = timeAfter - timeBefore;
	size_t totalCells = size_t(m_nTets) + m_nPyrs + m_nPrisms + m_nHexes;
	fprintf(stderr, "CPU time for UGRID file write = %5.2F seconds\n", elapsed);
	fprintf(stderr, "                          %5.2F million cells / minute\n",
			(totalCells / 1000000.) / (elapsed / 60));

	return true;
}

static void remapIndices(const emInt nPts, const std::vector<emInt> &newIndices,
		const emInt *conn, emInt *newConn) {
	for (emInt jj = 0; jj < nPts; jj++) {
		newConn[jj] = newIndices[conn[jj]];
	}
}

std::unique_ptr<UMesh> UMesh::createFineUMesh(const emInt numDivs, Part &P,
		std::vector<CellPartData> &vecCPD, struct RefineStats &RS) const {
	// Create a coarse
	double start = exaTime();
	auto coarse = extractCoarseMeshPseudoParallel(P, vecCPD, numDivs);
	double middle = exaTime();
	RS.extractTime = middle - start;

	auto UUM = coarse->subdivideMesh(numDivs);
	RS.cells = UUM->numCells();
	RS.refineTime = exaTime() - middle;
	return UUM;
}

void UMesh::setupCellDataForPartitioning(std::vector<CellPartData> &vecCPD,
		double &xmin, double &ymin, double &zmin, double &xmax, double &ymax,
		double &zmax) const {
	// Partitioning only cells, not bdry faces.  Also, currently no
	// cost differential for different cell types.
	for (emInt ii = 0; ii < numTets(); ii++) {
		const emInt *verts = getTetConn(ii);
		addCellToPartitionData(verts, 4, ii, CGNS_ENUMV(TETRA_4), vecCPD, xmin,
				ymin, zmin, xmax, ymax, zmax);
	}
	for (emInt ii = 0; ii < numPyramids(); ii++) {
		const emInt *verts = getPyrConn(ii);
		addCellToPartitionData(verts, 5, ii, CGNS_ENUMV(PYRA_5), vecCPD, xmin,
				ymin, zmin, xmax, ymax, zmax);
	}
	for (emInt ii = 0; ii < numPrisms(); ii++) {
		const emInt *verts = getPrismConn(ii);
		addCellToPartitionData(verts, 6, ii, CGNS_ENUMV(PENTA_6), vecCPD, xmin,
				ymin, zmin, xmax, ymax, zmax);
	}
	for (emInt ii = 0; ii < numHexes(); ii++) {
		const emInt *verts = getHexConn(ii);
		addCellToPartitionData(verts, 8, ii, CGNS_ENUMV(HEXA_8), vecCPD, xmin,
				ymin, zmin, xmax, ymax, zmax);
	}
}

std::unique_ptr<ExaMesh> UMesh::extractCoarseMeshPseudoParallel(Part &P,
		std::vector<CellPartData> &vecCPD, const int numDivs,
		const std::unordered_set<TriFaceVerts> &tris,
		const std::unordered_set<QuadFaceVerts> &quads,
		const emInt partID) const {
	// Count the number of tris, quads, tets, pyrs, prisms and hexes.
	const emInt first = P.getFirst();
	const emInt last = P.getLast();

	exa_set<TriFaceVerts> partBdryTris;
	exa_set<QuadFaceVerts> partBdryQuads;

	emInt nTris(0), nQuads(0), nTets(0), nPyrs(0), nPrisms(0), nHexes(0);
	const emInt *conn;

	std::vector<bool> isVertUsed(numVerts(), false);
	std::vector<bool> isBdryVert(numVerts(), false);

	for (emInt ii = first; ii < last; ii++) {
		emInt type = vecCPD[ii].getCellType();
		emInt ind = vecCPD[ii].getIndex();
		switch (type) {
		default:
			// Panic! Should never get here.
			assert(0);
			break;
		case CGNS_ENUMV(TRI_3):
			break;
		case CGNS_ENUMV(QUAD_4):
			break;
		case CGNS_ENUMV(TETRA_4): {
			nTets++;
			conn = getTetConn(ind);
			TriFaceVerts TFV012(numDivs, conn[0], conn[1], conn[2], type, ind);
			TriFaceVerts TFV013(numDivs, conn[0], conn[1], conn[3], type, ind);
			TriFaceVerts TFV123(numDivs, conn[1], conn[2], conn[3], type, ind);
			TriFaceVerts TFV203(numDivs, conn[2], conn[0], conn[3], type, ind);
			addUniquely(partBdryTris, TFV012);
			addUniquely(partBdryTris, TFV013);
			addUniquely(partBdryTris, TFV123);
			addUniquely(partBdryTris, TFV203);
			isVertUsed[conn[0]] = true;
			isVertUsed[conn[1]] = true;
			isVertUsed[conn[2]] = true;
			isVertUsed[conn[3]] = true;
			break;
		}
		case CGNS_ENUMV(PYRA_5): {
			nPyrs++;
			conn = getPyrConn(ind);
			QuadFaceVerts QFV0123(numDivs, conn[0], conn[1], conn[2], conn[3],
					type, ind);
			TriFaceVerts TFV014(numDivs, conn[0], conn[1], conn[4], type, ind);
			TriFaceVerts TFV124(numDivs, conn[1], conn[2], conn[4], type, ind);
			TriFaceVerts TFV234(numDivs, conn[2], conn[3], conn[4], type, ind);
			TriFaceVerts TFV304(numDivs, conn[3], conn[0], conn[4], type, ind);
			addUniquely(partBdryQuads, QFV0123);
			addUniquely(partBdryTris, TFV014);
			addUniquely(partBdryTris, TFV124);
			addUniquely(partBdryTris, TFV234);
			addUniquely(partBdryTris, TFV304);
			isVertUsed[conn[0]] = true;
			isVertUsed[conn[1]] = true;
			isVertUsed[conn[2]] = true;
			isVertUsed[conn[3]] = true;
			isVertUsed[conn[4]] = true;
			break;
		}
		case CGNS_ENUMV(PENTA_6): {
			nPrisms++;
			conn = getPrismConn(ind);
			QuadFaceVerts QFV0143(numDivs, conn[0], conn[1], conn[4], conn[3],
					type, ind);
			QuadFaceVerts QFV1254(numDivs, conn[1], conn[2], conn[5], conn[4],
					type, ind);
			QuadFaceVerts QFV2035(numDivs, conn[2], conn[0], conn[3], conn[5],
					type, ind);
			TriFaceVerts TFV012(numDivs, conn[0], conn[1], conn[2], type, ind);
			TriFaceVerts TFV345(numDivs, conn[3], conn[4], conn[5], type, ind);
			addUniquely(partBdryQuads, QFV0143);
			addUniquely(partBdryQuads, QFV1254);
			addUniquely(partBdryQuads, QFV2035);
			addUniquely(partBdryTris, TFV012);
			addUniquely(partBdryTris, TFV345);
			isVertUsed[conn[0]] = true;
			isVertUsed[conn[1]] = true;
			isVertUsed[conn[2]] = true;
			isVertUsed[conn[3]] = true;
			isVertUsed[conn[4]] = true;
			isVertUsed[conn[5]] = true;
			break;
		}
		case CGNS_ENUMV(HEXA_8): {
			nHexes++;
			conn = getHexConn(ind);
			QuadFaceVerts QFV0154(numDivs, conn[0], conn[1], conn[5], conn[4],
					type, ind);
			QuadFaceVerts QFV1265(numDivs, conn[1], conn[2], conn[6], conn[5],
					type, ind);
			QuadFaceVerts QFV2376(numDivs, conn[2], conn[3], conn[7], conn[6],
					type, ind);
			QuadFaceVerts QFV3047(numDivs, conn[3], conn[0], conn[4], conn[7],
					type, ind);
			QuadFaceVerts QFV0123(numDivs, conn[0], conn[1], conn[2], conn[3],
					type, ind);
			QuadFaceVerts QFV4567(numDivs, conn[4], conn[5], conn[6], conn[7],
					type, ind);
			addUniquely(partBdryQuads, QFV0154);
			addUniquely(partBdryQuads, QFV1265);
			addUniquely(partBdryQuads, QFV2376);
			addUniquely(partBdryQuads, QFV3047);
			addUniquely(partBdryQuads, QFV0123);
			addUniquely(partBdryQuads, QFV4567);
			isVertUsed[conn[0]] = true;
			isVertUsed[conn[1]] = true;
			isVertUsed[conn[2]] = true;
			isVertUsed[conn[3]] = true;
			isVertUsed[conn[4]] = true;
			isVertUsed[conn[5]] = true;
			isVertUsed[conn[6]] = true;
			isVertUsed[conn[7]] = true;
			break;
		}
		} // end switch
	} // end loop to gather information

	// Now check to see which bdry entities are in this part.  That'll be the
	// ones whose verts are all marked as used.  Unfortunately, this requires
	// searching through -all- the bdry entities for each part.
	std::vector<emInt> realBdryTris;
	std::vector<emInt> realBdryQuads;
	for (emInt ii = 0; ii < numBdryTris(); ii++) {
		conn = getBdryTriConn(ii);
		if (isVertUsed[conn[0]] && isVertUsed[conn[1]] && isVertUsed[conn[2]]) {
			TriFaceVerts TFV(numDivs, conn[0], conn[1], conn[2]);
			auto iter = partBdryTris.find(TFV);
			// If this bdry tri is an unmatched tri from this part, match it, and
			// add the bdry tri to the list of things to copy to the part coarse
			// mesh.  Otherwise, do nothing.  This will keep the occasional wrong
			// bdry face from slipping through.
			if (iter != partBdryTris.end()) {
				partBdryTris.erase(iter);
				isBdryVert[conn[0]] = true;
				isBdryVert[conn[1]] = true;
				isBdryVert[conn[2]] = true;
				realBdryTris.push_back(ii);
				nTris++;
			}
		}
	}
	for (emInt ii = 0; ii < numBdryQuads(); ii++) {
		conn = getBdryQuadConn(ii);
		if (isVertUsed[conn[0]] && isVertUsed[conn[1]] && isVertUsed[conn[2]]
				&& isVertUsed[conn[3]]) {
			QuadFaceVerts QFV(numDivs, conn[0], conn[1], conn[2], conn[3]);
			auto iter = partBdryQuads.find(QFV);
			// If this bdry tri is an unmatched tri from this part, match it, and
			// add the bdry tri to the list of things to copy to the part coarse
			// mesh.  Otherwise, do nothing.  This will keep the occasional wrong
			// bdry face from slipping through.
			if (iter != partBdryQuads.end()) {
				partBdryQuads.erase(iter);
				isBdryVert[conn[0]] = true;
				isBdryVert[conn[1]] = true;
				isBdryVert[conn[2]] = true;
				isBdryVert[conn[3]] = true;
				realBdryQuads.push_back(ii);
				nQuads++;
			}
		}
	}

	emInt nPartBdryTris = partBdryTris.size();
	emInt nPartBdryQuads = partBdryQuads.size();

	for (auto tri : partBdryTris) {
		isBdryVert[tri.getCorner(0)] = true;
		isBdryVert[tri.getCorner(1)] = true;
		isBdryVert[tri.getCorner(2)] = true;
	}
	for (auto quad : partBdryQuads) {
		isBdryVert[quad.getCorner(0)] = true;
		isBdryVert[quad.getCorner(1)] = true;
		isBdryVert[quad.getCorner(2)] = true;
		isBdryVert[quad.getCorner(3)] = true;
	}
	emInt nBdryVerts = 0, nVerts = 0;
	for (emInt ii = 0; ii < numVerts(); ii++) {
		if (isVertUsed[ii])
			nVerts++;
		if (isBdryVert[ii])
			nBdryVerts++;
	}

	// Now set up the data structures for the new coarse UMesh
	auto extractedMesh = std::make_unique<UMesh>(nVerts, nBdryVerts,
			nTris + nPartBdryTris, nQuads + nPartBdryQuads, nTets, nPyrs,
			nPrisms, nHexes);

	// Store the vertices, while keeping a mapping from the full list of verts
	// to the restricted list so the connectivity can be copied properly.
	std::vector<emInt> newIndices(numVerts(), EMINT_MAX);
	for (emInt ii = 0; ii < numVerts(); ii++) {
		if (isVertUsed[ii]) {
			double coords[3];
			getCoords(ii, coords);
			newIndices[ii] = extractedMesh->addVert(coords);
			// Copy length scale for vertices from the parent; otherwise, there will be
			// mismatches in the refined meshes.
			extractedMesh->setLengthScale(newIndices[ii], getLengthScale(ii));
		}
	}

	// Now copy connectivity.
	emInt newConn[8];
	for (emInt ii = first; ii < last; ii++) {
		emInt type = vecCPD[ii].getCellType();
		emInt ind = vecCPD[ii].getIndex();
		switch (type) {
		default:
			// Panic! Should never get here.
			assert(0);
			break;
		case CGNS_ENUMV(TRI_3):
			break;
		case CGNS_ENUMV(QUAD_4):
			break;
		case CGNS_ENUMV(TETRA_4): {
			conn = getTetConn(ind);
			remapIndices(4, newIndices, conn, newConn);
			extractedMesh->addTet(newConn);
			break;
		}
		case CGNS_ENUMV(PYRA_5): {
			conn = getPyrConn(ind);
			remapIndices(5, newIndices, conn, newConn);
			extractedMesh->addPyramid(newConn);
			break;
		}
		case CGNS_ENUMV(PENTA_6): {
			conn = getPrismConn(ind);
			remapIndices(6, newIndices, conn, newConn);
			extractedMesh->addPrism(newConn);
			break;
		}
		case CGNS_ENUMV(HEXA_8): {
			conn = getHexConn(ind);
			remapIndices(8, newIndices, conn, newConn);
			extractedMesh->addHex(newConn);
			break;
		}
		} // end switch
	} // end loop to copy most connectivity

	for (std::size_t ii = 0; ii < realBdryTris.size(); ii++) {
		conn = getBdryTriConn(realBdryTris[ii]);
		remapIndices(3, newIndices, conn, newConn);
		extractedMesh->addBdryTri(newConn);
	}
	for (std::size_t ii = 0; ii < realBdryQuads.size(); ii++) {
		conn = getBdryQuadConn(realBdryQuads[ii]);
		remapIndices(4, newIndices, conn, newConn);
		extractedMesh->addBdryQuad(newConn);
	}

	// Now, finally, the part bdry connectivity.
	// TODO: Currently, there's nothing in the data structure that marks which
	// are part bdry faces.
	assert(partBdryTris.size() == tris.size());

	for (auto tri : partBdryTris) {
		emInt localConn[] = { newIndices[tri.getCorner(0)],
				newIndices[tri.getCorner(1)], newIndices[tri.getCorner(2)] };
		emInt global[3] =
				{ tri.getCorner(0), tri.getCorner(1), tri.getCorner(2) };
		TriFaceVerts TF(numDivs, global, partID, -1, true);
		auto itr = tris.find(TF);
		if (itr != tris.end()) {
			assert(
					itr->getGlobalCorner(0) == global[0]
							&& itr->getGlobalCorner(1) == global[1]
							&& itr->getGlobalCorner(2) == global[2]
							&& itr->getPartid() == partID);
			TriFaceVerts TFV(numDivs, localConn, global, partID,
					itr->getRemoteId(), 0, EMINT_MAX, false);
			// need to be corrected, I could not generate with correct bool value unless
			// I pass all arguments

			extractedMesh->addPartTritoSet(TFV);
		}
		extractedMesh->addBdryTri(localConn);
	}
	assert(extractedMesh->getSizePartTris() == tris.size());

	assert(partBdryQuads.size() == quads.size());
	for (auto quad : partBdryQuads) {
		emInt localConn[] = { newIndices[quad.getCorner(0)],
				newIndices[quad.getCorner(1)], newIndices[quad.getCorner(2)],
				newIndices[quad.getCorner(3)] };
		emInt global[4] = { quad.getCorner(0), quad.getCorner(1),
				quad.getCorner(2), quad.getCorner(3) };
		QuadFaceVerts QF(numDivs, global, partID, -1, true);
		auto itr = quads.find(QF);
		if (itr != quads.end()) {
			assert(
					itr->getGlobalCorner(0) == global[0]
							&& itr->getGlobalCorner(1) == global[1]
							&& itr->getGlobalCorner(2) == global[2]
							&& itr->getGlobalCorner(3) == global[3]
							&& itr->getPartid() == partID);
			QuadFaceVerts QFV(numDivs, localConn, global, partID,
					itr->getRemoteId(), 0, EMINT_MAX, false);
			// need to be corrected, I could not generate with correct bool value unless
			// I pass all arguments
			extractedMesh->addPartQuadtoSet(QFV);
		}
		extractedMesh->addBdryQuad(localConn);
	}
	assert(extractedMesh->getSizePartQuads() == quads.size());
	return extractedMesh;
}

void UMesh::calcMemoryRequirements(const UMesh &UMIn, const int nDivs) {
	setlocale(LC_ALL, "");
	// Assuming that input cell is never on the order of billions, so use emInt 

	// size_t totalInputCells = size_t(UMIn.m_nTets) 
	// + UMIn.m_nPyrs + UMIn.m_nPrisms + UMIn.m_nHexes;

	emInt totalInputCells = UMIn.m_nTets + UMIn.m_nPyrs + UMIn.m_nPrisms
			+ UMIn.m_nHexes;

	// fprintf(
	// 	stderr,
	// 	"Initial mesh has:\n %'15u verts,\n %'15u bdry tris,\n %'15u bdry quads,\n %'15u tets,\n %'15u pyramids,\n %'15u prisms,\n %'15u hexes,\n%'15u cells total\n",
	// 	UMIn.m_nVerts, UMIn.m_nTris, UMIn.m_nQuads, UMIn.m_nTets, UMIn.m_nPyrs,
	// 	UMIn.m_nPrisms, UMIn.m_nHexes, totalInputCells);

	MeshSize MSOut = UMIn.computeFineMeshSize(nDivs);
	// since we want to get error value as negative so ssize_t instead of size_t 

	ssize_t nVerts = MSOut.nVerts;
	ssize_t nBdryVerts = MSOut.nBdryVerts;
	ssize_t nBdryTris = MSOut.nBdryTris;
	ssize_t nBdryQuads = MSOut.nBdryQuads;
	ssize_t nTets = MSOut.nTets;
	ssize_t nPyramids = MSOut.nPyrs;
	ssize_t nPrisms = MSOut.nPrisms;
	ssize_t nHexes = MSOut.nHexes;

	// All sizes are computed in bytes.
	// Work out buffer size, including padding to ensure 8-byte alignment for the coordinates.
	size_t intSize = sizeof(emInt);
	size_t headerSize = 7 * intSize;
	// How many bytes to add to get eight-byte alignment for coordinates,
	// assuming eight-byte alignment for the buffer overall.
	size_t slack1Size = (intSize == 4) ? 4 : 0;
	size_t coordSize = 3 * sizeof(double) * nVerts;
	size_t connSize = (3 * nBdryTris + 4 * nBdryQuads + 4 * size_t(nTets)
			+ 5 * size_t(nPyramids) + 6 * size_t(nPrisms) + 8 * size_t(nHexes))
			* intSize;
	size_t BCSize = (nBdryTris + nBdryQuads) * intSize;

	// How many bytes to add to fill up the last eight-byte chunk?
	size_t slack2Size = ((((connSize + BCSize) / 8 + 1) * 8)
			- (connSize + BCSize)) % 8;
	size_t bufferBytes = headerSize + coordSize + connSize + BCSize + slack1Size
			+ slack2Size;
	assert((headerSize + slack1Size) % 8 == 0);
	assert((connSize + BCSize + slack2Size) % 8 == 0);
	assert(bufferBytes % 8 == 0);

	size_t fileImageSize = bufferBytes - slack1Size - slack2Size;
	std::cout << "Final Image size: "
			<< static_cast<double>(fileImageSize) / 1000000000 << " GB"
			<< std::endl;
	setlocale(LC_ALL, "");
	size_t numCells = nTets + nPyramids + nPrisms + nHexes;
	fprintf(
	stderr,
			"Final mesh has:\n %'15zd verts,\n %'15zd bdry tris,\n %'15zd bdry quads,\n %'15zd tets,\n %'15zd pyramids,\n %'15zd prisms,\n %'15zd hexes,\n%'15zd cells total\n",
			nVerts, nBdryTris, nBdryQuads, nTets, nPyramids, nPrisms, nHexes,
			numCells);
}

std::unique_ptr<ExaMesh> UMesh::extractCoarseMeshMPI(const emInt partID,
		const std::vector<emInt> &partcells, const int numDivs,
		const std::unordered_set<TriFaceVerts>& tris,
		const std::unordered_set<QuadFaceVerts>& quads) const {
	// Count the number of tris, quads, tets, pyrs, prisms and hexes.
	// const emInt first = P.getFirst();
	// const emInt last = P.getLast();

	exa_set<TriFaceVerts> partBdryTris;
	exa_set<QuadFaceVerts> partBdryQuads;

	emInt nTris(0), nQuads(0), nTets(0), nPyrs(0), nPrisms(0), nHexes(0);
	const emInt *conn;

	std::vector<bool> isBdryVert(numVerts(), false);
	std::vector<bool> isVertUsed(numVerts(), false);

	for (emInt ii = 0; ii < partcells.size(); ii++) {
		// emInt type = vecCPD[ii].getCellType();
		// emInt ind = vecCPD[ii].getIndex();
		emInt globalInd = partcells[ii];
		emInt ind = cellID2cellTypeLocalID[globalInd].second;
		emInt type = cellID2cellTypeLocalID[globalInd].first;

		switch (type) {
		default:
			// Panic! Should never get here.
			assert(0);
			break;
		case CGNS_ENUMV(TRI_3):
			break;
		case CGNS_ENUMV(QUAD_4):
			break;
		case CGNS_ENUMV(TETRA_4): {
			nTets++;
			conn = getTetConn(ind);
			TriFaceVerts TFV012(numDivs, conn[0], conn[1], conn[2]);
			TriFaceVerts TFV013(numDivs, conn[0], conn[1], conn[3]);
			TriFaceVerts TFV123(numDivs, conn[1], conn[2], conn[3]);
			TriFaceVerts TFV203(numDivs, conn[2], conn[0], conn[3]);
			addUniquely(partBdryTris, TFV012);
			addUniquely(partBdryTris, TFV013);
			addUniquely(partBdryTris, TFV123);
			addUniquely(partBdryTris, TFV203);
			isVertUsed[conn[0]] = true;
			isVertUsed[conn[1]] = true;
			isVertUsed[conn[2]] = true;
			isVertUsed[conn[3]] = true;
			break;
		}
		case CGNS_ENUMV(PYRA_5): {
			nPyrs++;
			conn = getPyrConn(ind);
			QuadFaceVerts QFV0123(numDivs, conn[0], conn[1], conn[2], conn[3]);
			TriFaceVerts TFV014(numDivs, conn[0], conn[1], conn[4]);
			TriFaceVerts TFV124(numDivs, conn[1], conn[2], conn[4]);
			TriFaceVerts TFV234(numDivs, conn[2], conn[3], conn[4]);
			TriFaceVerts TFV304(numDivs, conn[3], conn[0], conn[4]);
			addUniquely(partBdryQuads, QFV0123);
			addUniquely(partBdryTris, TFV014);
			addUniquely(partBdryTris, TFV124);
			addUniquely(partBdryTris, TFV234);
			addUniquely(partBdryTris, TFV304);
			isVertUsed[conn[0]] = true;
			isVertUsed[conn[1]] = true;
			isVertUsed[conn[2]] = true;
			isVertUsed[conn[3]] = true;
			isVertUsed[conn[4]] = true;
			break;
		}
		case CGNS_ENUMV(PENTA_6): {
			nPrisms++;
			conn = getPrismConn(ind);
			QuadFaceVerts QFV0143(numDivs, conn[0], conn[1], conn[4], conn[3]);
			QuadFaceVerts QFV1254(numDivs, conn[1], conn[2], conn[5], conn[4]);
			QuadFaceVerts QFV2035(numDivs, conn[2], conn[0], conn[3], conn[5]);
			TriFaceVerts TFV012(numDivs, conn[0], conn[1], conn[2]);
			TriFaceVerts TFV345(numDivs, conn[3], conn[4], conn[5]);
			addUniquely(partBdryQuads, QFV0143);
			addUniquely(partBdryQuads, QFV1254);
			addUniquely(partBdryQuads, QFV2035);
			addUniquely(partBdryTris, TFV012);
			addUniquely(partBdryTris, TFV345);
			isVertUsed[conn[0]] = true;
			isVertUsed[conn[1]] = true;
			isVertUsed[conn[2]] = true;
			isVertUsed[conn[3]] = true;
			isVertUsed[conn[4]] = true;
			isVertUsed[conn[5]] = true;
			break;
		}
		case CGNS_ENUMV(HEXA_8): {
			nHexes++;
			conn = getHexConn(ind);
			QuadFaceVerts QFV0154(numDivs, conn[0], conn[1], conn[5], conn[4]);
			QuadFaceVerts QFV1265(numDivs, conn[1], conn[2], conn[6], conn[5]);
			QuadFaceVerts QFV2376(numDivs, conn[2], conn[3], conn[7], conn[6]);
			QuadFaceVerts QFV3047(numDivs, conn[3], conn[0], conn[4], conn[7]);
			QuadFaceVerts QFV0123(numDivs, conn[0], conn[1], conn[2], conn[3]);
			QuadFaceVerts QFV4567(numDivs, conn[4], conn[5], conn[6], conn[7]);
			addUniquely(partBdryQuads, QFV0154);
			addUniquely(partBdryQuads, QFV1265);
			addUniquely(partBdryQuads, QFV2376);
			addUniquely(partBdryQuads, QFV3047);
			addUniquely(partBdryQuads, QFV0123);
			addUniquely(partBdryQuads, QFV4567);
			isVertUsed[conn[0]] = true;
			isVertUsed[conn[1]] = true;
			isVertUsed[conn[2]] = true;
			isVertUsed[conn[3]] = true;
			isVertUsed[conn[4]] = true;
			isVertUsed[conn[5]] = true;
			isVertUsed[conn[6]] = true;
			isVertUsed[conn[7]] = true;
			break;
		}
		} // end switch
	}	  // end loop to gather information

	// Now check to see which bdry entities are in this part.  That'll be the
	// ones whose verts are all marked as used.  Unfortunately, this requires
	// searching through -all- the bdry entities for each part.
	std::vector<emInt> realBdryTris;
	std::vector<emInt> realBdryQuads;
	for (emInt ii = 0; ii < numBdryTris(); ii++) {
		conn = getBdryTriConn(ii);
		if (isVertUsed[conn[0]] && isVertUsed[conn[1]] && isVertUsed[conn[2]]) {
			TriFaceVerts TFV(numDivs, conn[0], conn[1], conn[2]);
			auto iter = partBdryTris.find(TFV);
			// If this bdry tri is an unmatched tri from this part, match it, and
			// add the bdry tri to the list of things to copy to the part coarse
			// mesh.  Otherwise, do nothing.  This will keep the occasional wrong
			// bdry face from slipping through.
			if (iter != partBdryTris.end()) {
				partBdryTris.erase(iter);
				isBdryVert[conn[0]] = true;
				isBdryVert[conn[1]] = true;
				isBdryVert[conn[2]] = true;
				realBdryTris.push_back(ii);
				nTris++;
			}
		}
	}
	for (emInt ii = 0; ii < numBdryQuads(); ii++) {
		conn = getBdryQuadConn(ii);
		if (isVertUsed[conn[0]] && isVertUsed[conn[1]] && isVertUsed[conn[2]]
				&& isVertUsed[conn[3]]) {
			QuadFaceVerts QFV(numDivs, conn[0], conn[1], conn[2], conn[3]);
			auto iter = partBdryQuads.find(QFV);
			// If this bdry tri is an unmatched tri from this part, match it, and
			// add the bdry tri to the list of things to copy to the part coarse
			// mesh.  Otherwise, do nothing.  This will keep the occasional wrong
			// bdry face from slipping through.
			if (iter != partBdryQuads.end()) {
				partBdryQuads.erase(iter);
				isBdryVert[conn[0]] = true;
				isBdryVert[conn[1]] = true;
				isBdryVert[conn[2]] = true;
				isBdryVert[conn[3]] = true;
				realBdryQuads.push_back(ii);
				nQuads++;
			}
		}
	}

	emInt nPartBdryTris = partBdryTris.size();
	emInt nPartBdryQuads = partBdryQuads.size();

	for (auto tri : partBdryTris) {
		isBdryVert[tri.getCorner(0)] = true;
		isBdryVert[tri.getCorner(1)] = true;
		isBdryVert[tri.getCorner(2)] = true;
	}

	for (auto quad : partBdryQuads) {
		isBdryVert[quad.getCorner(0)] = true;
		isBdryVert[quad.getCorner(1)] = true;
		isBdryVert[quad.getCorner(2)] = true;
		isBdryVert[quad.getCorner(3)] = true;
	}
	emInt nBdryVerts = 0, nVerts = 0;
	for (emInt ii = 0; ii < numVerts(); ii++) {
		if (isBdryVert[ii])
			nBdryVerts++;
		if (isVertUsed[ii])
			nVerts++;
	}

	// Now set up the data structures for the new coarse UMesh
	auto UUM = std::make_unique<UMesh>(nVerts, nBdryVerts,
			nTris + nPartBdryTris, nQuads + nPartBdryQuads, nTets, nPyrs,
			nPrisms, nHexes);

	// Store the vertices, while keeping a mapping from the full list of verts
	// to the restricted list so the connectivity can be copied properly.
	std::vector<emInt> newIndices(numVerts(), EMINT_MAX);
	for (emInt ii = 0; ii < numVerts(); ii++) {
		if (isVertUsed[ii]) {
			double coords[3];
			getCoords(ii, coords);
			newIndices[ii] = UUM->addVert(coords);
			// Copy length scale for vertices from the parent; otherwise, there will be
			// mismatches in the refined meshes.
			UUM->setLengthScale(newIndices[ii], getLengthScale(ii));
		}
	}

	// Now copy connectivity.
	emInt newConn[8];
	for (emInt ii = 0; ii < partcells.size(); ii++) {
		//emInt type = vecCPD[ii].getCellType();
		//emInt ind = vecCPD[ii].getIndex();
		emInt globalInd = partcells[ii];
		emInt ind = cellID2cellTypeLocalID[globalInd].second;
		emInt type = cellID2cellTypeLocalID[globalInd].first;
		switch (type) {
		default:
			// Panic! Should never get here.
			assert(0);
			break;
		case CGNS_ENUMV(TRI_3):
			break;
		case CGNS_ENUMV(QUAD_4):
			break;
		case CGNS_ENUMV(TETRA_4): {
			conn = getTetConn(ind);
			remapIndices(4, newIndices, conn, newConn);
			UUM->addTet(newConn);
			break;
		}
		case CGNS_ENUMV(PYRA_5): {
			conn = getPyrConn(ind);
			remapIndices(5, newIndices, conn, newConn);
			UUM->addPyramid(newConn);
			break;
		}
		case CGNS_ENUMV(PENTA_6): {
			conn = getPrismConn(ind);
			remapIndices(6, newIndices, conn, newConn);
			UUM->addPrism(newConn);
			break;
		}
		case CGNS_ENUMV(HEXA_8): {
			conn = getHexConn(ind);
			remapIndices(8, newIndices, conn, newConn);
			UUM->addHex(newConn);
			break;
		}
		} // end switch
	}	  // end loop to copy most connectivity

	for (std::size_t ii = 0; ii < realBdryTris.size(); ii++) {
		conn = getBdryTriConn(realBdryTris[ii]);
		remapIndices(3, newIndices, conn, newConn);
		UUM->addBdryTri(newConn);
	}
	for (std::size_t ii = 0; ii < realBdryQuads.size(); ii++) {
		conn = getBdryQuadConn(realBdryQuads[ii]);
		remapIndices(4, newIndices, conn, newConn);
		UUM->addBdryQuad(newConn);
	}

	// Now, finally, the part bdry connectivity.
	// TODO: Currently, there's nothing in the data structure that marks which
	// are part bdry faces.

	assert(partBdryTris.size() == tris.size());
	for (auto tri : partBdryTris) {
		emInt conn[] = { newIndices[tri.getCorner(0)], newIndices[tri.getCorner(
				1)], newIndices[tri.getCorner(2)] };
		emInt global[3] =
				{ tri.getCorner(0), tri.getCorner(1), tri.getCorner(2) };
		TriFaceVerts TF(numDivs, global, partID, -1, true);
		auto itr = tris.find(TF);
		if (itr != tris.end()) {
			assert(
					itr->getGlobalCorner(0) == global[0]
							&& itr->getGlobalCorner(1) == global[1]
							&& itr->getGlobalCorner(2) == global[2]
							&& itr->getPartid() == partID);
			TriFaceVerts TFV(numDivs, conn, global, partID, itr->getRemoteId(),
					0, EMINT_MAX, false);
			//need to be corrected, I could not generate with correct bool value unless
			//I pass all arguments

			UUM->addPartTritoSet(TFV);
		}

		UUM->addBdryTri(conn);
	}

	assert(UUM->getSizePartTris() == tris.size());

	assert(partBdryQuads.size() == quads.size());
	for (auto quad : partBdryQuads) {
		emInt conn[] = { newIndices[quad.getCorner(0)],
				newIndices[quad.getCorner(1)], newIndices[quad.getCorner(2)],
				newIndices[quad.getCorner(3)] };
		emInt global[4] = { quad.getCorner(0), quad.getCorner(1),
				quad.getCorner(2), quad.getCorner(3) };
		QuadFaceVerts QF(numDivs, global, partID, -1, true);
		auto itr = quads.find(QF);
		if (itr != quads.end()) {
			assert(
					itr->getGlobalCorner(0) == global[0]
							&& itr->getGlobalCorner(1) == global[1]
							&& itr->getGlobalCorner(2) == global[2]
							&& itr->getGlobalCorner(3) == global[3]
							&& itr->getPartid() == partID);
			QuadFaceVerts QFV(numDivs, conn, global, partID, itr->getRemoteId(),
					0, EMINT_MAX, false);
			//need to be corrected, I could not generate with correct bool value unless
			//I pass all arguments
			UUM->addPartQuadtoSet(QFV);
		}

		UUM->addBdryQuad(conn);
	}
	assert(UUM->getSizePartQuads() == quads.size());
	return UUM;
}

