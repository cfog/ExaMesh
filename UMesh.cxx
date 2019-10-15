/*
 * UMesh.cxx
 *
 *  Created on: Jul. 2, 2019
 *      Author: cfog
 */

#include <iostream>
#include <cmath>

#include <algorithm>
#include <set>
#include <vector>

#include <string.h>
#include <zlib.h>

#include "examesh.h"
#include "UMesh.h"

#include "GMGW_unstr.hxx"
#include "GMGW_FileWrapper.hxx"

using std::cout;
using std::endl;

static bool memoryCheck(void* address, int nBytes) {
	char* checkPtr = reinterpret_cast<char*>(address);
	bool retVal = true;
	for (int ii = 0; ii < nBytes; ii++) {
		retVal = retVal && (checkPtr[ii] = 0xff);
	}
	return retVal;
}

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
											+ 5 * size_t(nPyramids) + 6 * size_t(nPrisms)
											+ 8 * size_t(nHexes))
										* intSize;
	// How many bytes to add to fill up the last eight-byte chunk?
	size_t slack2Size = ((connSize / 8 + 1) * 8) - connSize;
	size_t bufferBytes = headerSize + coordSize + connSize + slack1Size
											+ slack2Size;
	assert((headerSize + slack1Size) % 8 == 0);
	assert((connSize + slack2Size) % 8 == 0);
	assert(bufferBytes % 8 == 0);
	size_t bufferWords = bufferBytes / 8;
	// Use words instead of bytes to ensure 8-byte alignment.
	m_buffer = reinterpret_cast<char*>(calloc(bufferWords, 8));
	memset(m_buffer, 0xFF, bufferBytes);
	// The pointer arithmetic here is made more complicated because the pointers aren't
	// compatible with each other.
	m_header = reinterpret_cast<emInt*>(m_buffer + slack1Size);
	std::fill(m_header, m_header + 7, 0);
	m_coords =
			reinterpret_cast<double (*)[3]>(m_buffer + headerSize + slack1Size);
	m_TriConn = reinterpret_cast<emInt (*)[3]>(m_buffer + headerSize + slack1Size
																							+ coordSize);
	m_QuadConn = reinterpret_cast<emInt (*)[4]>(m_TriConn + nBdryTris);
	m_TetConn = m_QuadConn + nBdryQuads;
	m_PyrConn = reinterpret_cast<emInt (*)[5]>(m_TetConn + nTets);
	m_PrismConn = reinterpret_cast<emInt (*)[6]>(m_PyrConn + nPyramids);
	m_HexConn = reinterpret_cast<emInt (*)[8]>(m_PrismConn + nPrisms);
	m_fileImage = m_buffer + slack1Size;
	m_fileImageSize = bufferBytes - slack1Size - slack2Size;

	m_lenScale = new double[m_nVerts];
}

UMesh::UMesh(const emInt nVerts, const emInt nBdryVerts, const emInt nBdryTris,
		const emInt nBdryQuads, const emInt nTets, const emInt nPyramids,
		const emInt nPrisms, const emInt nHexes) :
		m_nVerts(nVerts), m_nBdryVerts(nBdryVerts), m_nTris(nBdryTris),
				m_nQuads(nBdryQuads), m_nTets(nTets), m_nPyrs(nPyramids),
				m_nPrisms(nPrisms), m_nHexes(nHexes), m_fileImageSize(0),
				m_header(nullptr), m_coords(nullptr), m_TriConn(nullptr),
				m_QuadConn(nullptr), m_TetConn(nullptr), m_PyrConn(nullptr),
				m_PrismConn(nullptr), m_HexConn(nullptr), m_buffer(nullptr),
				m_fileImage(nullptr) {

	// All sizes are computed in bytes.

	// Work out buffer size, including padding to ensure 8-byte alignment for the coordinates.
	init(nVerts, nBdryVerts, nBdryTris, nBdryQuads, nTets, nPyramids, nPrisms,
				nHexes);
}

emInt UMesh::addVert(const double newCoords[3]) {
	assert(memoryCheck(m_coords[m_header[eVert]], 24));
	m_coords[m_header[eVert]][0] = newCoords[0];
	m_coords[m_header[eVert]][1] = newCoords[1];
	m_coords[m_header[eVert]][2] = newCoords[2];
	return (m_header[eVert]++);
}

emInt UMesh::addBdryTri(const emInt verts[3]) {
	assert(memoryCheck(m_TriConn[m_header[eTri]], 3 * sizeof(emInt)));
	for (int ii = 0; ii < 3; ii++) {
		m_TriConn[m_header[eTri]][ii] = verts[ii];
	}
	return (m_header[eTri]++);
}

emInt UMesh::addBdryQuad(const emInt verts[4]) {
	assert(memoryCheck(m_QuadConn[m_header[eQuad]], 4 * sizeof(emInt)));
	for (int ii = 0; ii < 4; ii++) {
		m_QuadConn[m_header[eQuad]][ii] = verts[ii];
	}
	return (m_header[eQuad]++);
}

emInt UMesh::addTet(const emInt verts[4]) {
	assert(memoryCheck(m_TetConn[m_header[eTet]], 4 * sizeof(emInt)));
	for (int ii = 0; ii < 4; ii++) {
		m_TetConn[m_header[eTet]][ii] = verts[ii];
	}
	return (m_header[eTet]++);
}

emInt UMesh::addPyramid(const emInt verts[5]) {
	assert(memoryCheck(m_PyrConn[m_header[ePyr]], 5 * sizeof(emInt)));
	for (int ii = 0; ii < 5; ii++) {
		m_PyrConn[m_header[ePyr]][ii] = verts[ii];
	}
	return (m_header[ePyr]++);
}

emInt UMesh::addPrism(const emInt verts[6]) {
	assert(memoryCheck(m_PrismConn[m_header[ePrism]], 6 * sizeof(emInt)));
	for (int ii = 0; ii < 6; ii++) {
		m_PrismConn[m_header[ePrism]][ii] = verts[ii];
	}
	return (m_header[ePrism]++);
}

emInt UMesh::addHex(const emInt verts[8]) {
	assert(memoryCheck(m_HexConn[m_header[eHex]], 8 * sizeof(emInt)));
	for (int ii = 0; ii < 8; ii++) {
		m_HexConn[m_header[eHex]][ii] = verts[ii];
	}
	return (m_header[eHex]++);
}

UMesh::~UMesh() {
	free(m_buffer);
}

void checkConnectivitySize(const char cellType, const emInt nVerts) {
	emInt expected[] = { 0, 0, 0, 0, 0, 3, 0, 0, 0, 4, 4, 0, 8, 6, 5 };
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
	bool operator<(const vertTriple& that) const {
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
	bool operator<(const vertQuadruple& that) const {
		emInt thisTemp[4], thatTemp[4];
		sortVerts4(corners, thisTemp);
		sortVerts4(that.corners, thatTemp);
		return (thisTemp[0] < thatTemp[0]
				|| (thisTemp[0] == thatTemp[0] && thisTemp[1] < thatTemp[1])
				|| (thisTemp[0] == thatTemp[0] && thisTemp[1] == thatTemp[1]
						&& thisTemp[2] < thatTemp[2])
				|| (thisTemp[0] == thatTemp[0] && thisTemp[1] == thatTemp[1]
						&& thisTemp[2] == thatTemp[2] && thisTemp[3] < thatTemp[3]));
	}
	const emInt* getCorners() const {
		return corners;
	}
};

void updateTriSet(std::set<vertTriple>& triSet, const emInt v0, const emInt v1,
		const emInt v2) {
	vertTriple VT(v0, v1, v2);
	typename std::set<vertTriple>::iterator vertIter, VIend = triSet.end();

	vertIter = triSet.find(VT);
	if (vertIter == VIend) {
		triSet.insert(VT);
	}
	else {
		triSet.erase(vertIter);
	}
}

void updateQuadSet(std::set<vertQuadruple>& quadSet, const emInt v0,
		const emInt v1, const emInt v2, const emInt v3) {
	vertQuadruple VQ(v0, v1, v2, v3);
	typename std::set<vertQuadruple>::iterator vertIter, VQend = quadSet.end();

	vertIter = quadSet.find(VQ);
	if (vertIter == VQend) {
		quadSet.insert(VQ);
	}
	else {
		quadSet.erase(vertIter);
	}
}

UMesh::UMesh(const char baseFileName[], const char type[],
		const char ugridInfix[]) :
		m_nVerts(0), m_nBdryVerts(0), m_nTris(0), m_nQuads(0), m_nTets(0),
				m_nPyrs(0), m_nPrisms(0), m_nHexes(0), m_fileImageSize(0),
				m_header(nullptr), m_coords(nullptr), m_TriConn(nullptr),
				m_QuadConn(nullptr), m_TetConn(nullptr), m_PyrConn(nullptr),
				m_PrismConn(nullptr), m_HexConn(nullptr), m_buffer(nullptr),
				m_fileImage(nullptr) {
	// Use the same IO routines as the mesh analyzer code from GMGW.
	FileWrapper* reader = FileWrapper::factory(baseFileName, type,
																										ugridInfix);

	reader->scanFile();

	// Identify any bdry tris and quads that aren't in the file

	emInt numBdryTris = reader->getNumBdryTris();
	emInt numBdryQuads = reader->getNumBdryQuads();

	std::set<vertTriple> setTris;
	std::set<vertQuadruple> setQuads;

	reader->seekStartOfConnectivity();
	for (emInt ii = 0; ii < reader->getNumCells(); ii++) {
		char cellType = reader->getCellType(ii);
		emInt nConn, connect[8];
		reader->getNextCellConnectivity(nConn, connect);
		checkConnectivitySize(cellType, nConn);
		switch (cellType) {
			case BDRY_TRI:
				updateTriSet(setTris, connect[0], connect[1], connect[2]);
				break;
			case BDRY_QUAD:
				updateQuadSet(setQuads, connect[0], connect[1], connect[2], connect[3]);
				break;
			case TET:
				updateTriSet(setTris, connect[0], connect[1], connect[2]);
				updateTriSet(setTris, connect[0], connect[1], connect[3]);
				updateTriSet(setTris, connect[1], connect[2], connect[3]);
				updateTriSet(setTris, connect[2], connect[0], connect[3]);
				break;
			case PYRAMID:
				updateTriSet(setTris, connect[0], connect[1], connect[4]);
				updateTriSet(setTris, connect[1], connect[2], connect[4]);
				updateTriSet(setTris, connect[2], connect[3], connect[4]);
				updateTriSet(setTris, connect[3], connect[0], connect[4]);
				updateQuadSet(setQuads, connect[0], connect[1], connect[2], connect[3]);
				break;
			case PRISM:
				updateTriSet(setTris, connect[0], connect[1], connect[2]);
				updateTriSet(setTris, connect[3], connect[4], connect[5]);
				updateQuadSet(setQuads, connect[0], connect[1], connect[4], connect[3]);
				updateQuadSet(setQuads, connect[1], connect[2], connect[5], connect[4]);
				updateQuadSet(setQuads, connect[2], connect[0], connect[3], connect[5]);
				break;
			case HEX:
				updateQuadSet(setQuads, connect[0], connect[1], connect[2], connect[3]);
				updateQuadSet(setQuads, connect[4], connect[5], connect[6], connect[7]);
				updateQuadSet(setQuads, connect[0], connect[1], connect[5], connect[4]);
				updateQuadSet(setQuads, connect[1], connect[2], connect[6], connect[5]);
				updateQuadSet(setQuads, connect[2], connect[3], connect[7], connect[6]);
				updateQuadSet(setQuads, connect[3], connect[0], connect[4], connect[7]);
				break;
			default:
				assert(0);
		}
	}

	numBdryTris += setTris.size();
	numBdryQuads += setQuads.size();

	init(reader->getNumVerts(), reader->getNumBdryVerts(), numBdryTris,
				numBdryQuads, reader->getNumTets(), reader->getNumPyramids(),
				reader->getNumPrisms(), reader->getNumHexes());

	reader->seekStartOfCoords();
	for (emInt ii = 0; ii < m_nVerts; ii++) {
		double coords[3];
		reader->getNextVertexCoords(coords[0], coords[1], coords[2]);
		addVert(coords);
	}

	reader->seekStartOfConnectivity();
	for (emInt ii = 0; ii < reader->getNumCells(); ii++) {
		char cellType = reader->getCellType(ii);
		emInt nConn, connect[8];
		reader->getNextCellConnectivity(nConn, connect);
		checkConnectivitySize(cellType, nConn);
		switch (cellType) {
			case BDRY_TRI:
				addBdryTri(connect);
				break;
			case BDRY_QUAD:
				addBdryQuad(connect);
				break;
			case TET:
				addTet(connect);
				break;
			case PYRAMID:
				addPyramid(connect);
				break;
			case PRISM:
				addPrism(connect);
				break;
			case HEX:
				addHex(connect);
				break;
			default:
				assert(0);
		}
	}

	for (auto VT : setTris) {
		const emInt* const corners = VT.getCorners();
		addBdryTri(corners);
	}
	for (auto VQ : setQuads) {
		const emInt* const corners = VQ.getCorners();
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

	// If any of these fail, your file was invalid.
	assert(m_nVerts == m_header[eVert]);
	assert(m_nTris == m_header[eTri]);
	assert(m_nQuads == m_header[eQuad]);
	assert(m_nTets == m_header[eTet]);
	assert(m_nPyrs == m_header[ePyr]);
	assert(m_nPrisms == m_header[ePrism]);
	assert(m_nHexes == m_header[eHex]);

	delete reader;

	setupLengthScales();
}

UMesh::UMesh(const UMesh& UMIn, const int nDivs) :
		m_nVerts(0), m_nBdryVerts(0), m_nTris(0), m_nQuads(0), m_nTets(0),
				m_nPyrs(0), m_nPrisms(0), m_nHexes(0), m_fileImageSize(0),
				m_header(nullptr), m_coords(nullptr), m_TriConn(nullptr),
				m_QuadConn(nullptr), m_TetConn(nullptr), m_PyrConn(nullptr),
				m_PrismConn(nullptr), m_HexConn(nullptr), m_buffer(nullptr),
				m_fileImage(nullptr) {

	setlocale(LC_ALL, "");
	size_t totalCells = size_t(UMIn.m_nTets) + UMIn.m_nPyrs + UMIn.m_nPrisms
											+ UMIn.m_nHexes;
	fprintf(
			stderr,
			"Initial mesh has:\n %'15u verts,\n %'15u bdry tris,\n %'15u bdry quads,\n %'15u tets,\n %'15u pyramids,\n %'15u prisms,\n %'15u hexes,\n%'15lu cells total\n",
			UMIn.m_nVerts, UMIn.m_nTris, UMIn.m_nQuads, UMIn.m_nTets, UMIn.m_nPyrs,
			UMIn.m_nPrisms, UMIn.m_nHexes, totalCells);

	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = UMIn.m_nBdryVerts;
	MSIn.nVerts = UMIn.m_header[eVert];
	MSIn.nBdryTris = UMIn.m_header[eTri];
	MSIn.nBdryQuads = UMIn.m_header[eQuad];
	MSIn.nTets = UMIn.m_header[eTet];
	MSIn.nPyrs = UMIn.m_header[ePyr];
	MSIn.nPrisms = UMIn.m_header[ePrism];
	MSIn.nHexes = UMIn.m_header[eHex];
	bool sizesOK = computeMeshSize(MSIn, nDivs, MSOut);
	if (!sizesOK) exit(2);

	init(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris, MSOut.nBdryQuads,
				MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms, MSOut.nHexes);
	// Copy length scale data from the other mesh.
	for (emInt vv = 0; vv < UMIn.m_nVerts; vv++) {
		m_lenScale[vv] = UMIn.m_lenScale[vv];
	}

	double timeBefore = clock() / double(CLOCKS_PER_SEC);
	subdividePartMesh(&UMIn, this, nDivs);
	double timeAfter = clock() / double(CLOCKS_PER_SEC);
	double elapsed = timeAfter - timeBefore;
	setlocale(LC_ALL, "");
	totalCells = size_t(m_nTets) + m_nPyrs + m_nPrisms + m_nHexes;
	fprintf(
			stderr,
			"Final mesh has:\n %'15u verts,\n %'15u bdry tris,\n %'15u bdry quads,\n %'15u tets,\n %'15u pyramids,\n %'15u prisms,\n %'15u hexes,\n%'15lu cells total\n",
			m_nVerts, m_nTris, m_nQuads, m_nTets, m_nPyrs, m_nPrisms, m_nHexes,
			totalCells);
	fprintf(stderr, "CPU time for refinement = %5.2F seconds\n", elapsed);
	fprintf(stderr, "                          %5.2F million cells / minute\n",
					(totalCells / 1000000.) / (elapsed / 60));
}

UMesh::UMesh(const CubicMesh& CMIn, const int nDivs) :
		m_nVerts(0), m_nBdryVerts(0), m_nTris(0), m_nQuads(0), m_nTets(0),
				m_nPyrs(0), m_nPrisms(0), m_nHexes(0), m_fileImageSize(0),
				m_header(nullptr), m_coords(nullptr), m_TriConn(nullptr),
				m_QuadConn(nullptr), m_TetConn(nullptr), m_PyrConn(nullptr),
				m_PrismConn(nullptr), m_HexConn(nullptr), m_buffer(nullptr),
				m_fileImage(nullptr) {

	setlocale(LC_ALL, "");
	size_t totalCells = size_t(CMIn.numTets()) + CMIn.numPyramids()
											+ CMIn.numPrisms() + CMIn.numHexes();
	fprintf(
			stderr,
			"Initial mesh has:\n %'15u verts,\n %'15u bdry tris,\n %'15u bdry quads,\n %'15u tets,\n %'15u pyramids,\n %'15u prisms,\n %'15u hexes,\n%'15lu cells total\n",
			CMIn.numVerts(), CMIn.numBdryTris(), CMIn.numBdryQuads(), CMIn.numTets(),
			CMIn.numPyramids(), CMIn.numPrisms(), CMIn.numHexes(), totalCells);

	MeshSize MSIn, MSOut;
	MSIn.nBdryVerts = CMIn.numBdryVerts();
	MSIn.nVerts = CMIn.numVerts();
	MSIn.nBdryTris = CMIn.numBdryTris();
	MSIn.nBdryQuads = CMIn.numBdryQuads();
	MSIn.nTets = CMIn.numTets();
	MSIn.nPyrs = CMIn.numPyramids();
	MSIn.nPrisms = CMIn.numPrisms();
	MSIn.nHexes = CMIn.numHexes();
	bool sizesOK = computeMeshSize(MSIn, nDivs, MSOut);
	if (!sizesOK) exit(2);

	init(MSOut.nVerts, MSOut.nBdryVerts, MSOut.nBdryTris, MSOut.nBdryQuads,
				MSOut.nTets, MSOut.nPyrs, MSOut.nPrisms, MSOut.nHexes);
	// Copy length scale data from the other mesh.
	for (emInt vv = 0; vv < CMIn.numVerts(); vv++) {
		m_lenScale[vv] = CMIn.getLengthScale(vv);
	}

	double timeBefore = clock() / double(CLOCKS_PER_SEC);
	subdividePartMesh(&CMIn, this, nDivs);
	double timeAfter = clock() / double(CLOCKS_PER_SEC);
	double elapsed = timeAfter - timeBefore;
	setlocale(LC_ALL, "");
	totalCells = size_t(m_nTets) + m_nPyrs + m_nPrisms + m_nHexes;
	fprintf(
			stderr,
			"Final mesh has:\n %'15u verts,\n %'15u bdry tris,\n %'15u bdry quads,\n %'15u tets,\n %'15u pyramids,\n %'15u prisms,\n %'15u hexes,\n%'15lu cells total\n",
			m_nVerts, m_nTris, m_nQuads, m_nTets, m_nPyrs, m_nPrisms, m_nHexes,
			totalCells);
	fprintf(stderr, "CPU time for refinement = %5.2F seconds\n", elapsed);
	fprintf(stderr, "                          %5.2F million cells / minute\n",
					(totalCells / 1000000.) / (elapsed / 60));
}


bool UMesh::writeVTKFile(const char fileName[]) {
	double timeBefore = clock() / double(CLOCKS_PER_SEC);

	FILE* outFile = fopen(fileName, "w");
	if (!outFile) {
		fprintf(stderr, "Couldn't open file %s for writing.  Bummer!\n", fileName);
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
	for (emInt i = 0; i < nTris; i++) {
		const emInt *verts = getBdryTriConn(i);
		fprintf(outFile, "3 %d %d %d\n", verts[0], verts[1], verts[2]);
	}

	// Write all the bdry quads
	for (emInt i = 0; i < nQuads; i++) {
		const emInt *verts = getBdryQuadConn(i);
		fprintf(outFile, "4 %d %d %d %d\n", verts[0], verts[1], verts[2],
						verts[3]);
	}

	// Write all the tets
	for (emInt i = 0; i < nTets; i++) {
		const emInt *verts = getTetConn(i);
		fprintf(outFile, "4 %d %d %d %d\n", verts[0], verts[1], verts[2],
						verts[3]);
	}

	// Write all the pyramids
	for (emInt i = 0; i < nPyrs; i++) {
		const emInt *verts = getPyrConn(i);
		fprintf(outFile, "5 %d %d %d %d %d\n", verts[0], verts[1], verts[2],
						verts[3], verts[4]);
	}

	// Write all the prisms
	for (emInt i = 0; i < nPrisms; i++) {
		const emInt *verts = getPrismConn(i);
		fprintf(outFile, "6 %d %d %d %d %d %d\n", verts[0], verts[1], verts[2],
						verts[3], verts[4], verts[5]);
	}

	// Write all the hexes
	for (emInt i = 0; i < nHexes; i++) {
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
	for (emInt ct = 0; ct < nTris; ++ct)
		fprintf(outFile, "5\n");
	for (emInt ct = 0; ct < nQuads; ++ct)
		fprintf(outFile, "9\n");
	for (emInt ct = 0; ct < nTets; ++ct)
		fprintf(outFile, "10\n");
	for (emInt ct = 0; ct < nPyrs; ++ct)
		fprintf(outFile, "14\n");
	for (emInt ct = 0; ct < nPrisms; ++ct)
		fprintf(outFile, "13\n");
	for (emInt ct = 0; ct < nHexes; ++ct)
		fprintf(outFile, "12\n");

	fclose(outFile);
	double timeAfter = clock() / double(CLOCKS_PER_SEC);
	double elapsed = timeAfter - timeBefore;
	size_t totalCells = size_t(m_nTets) + m_nPyrs + m_nPrisms + m_nHexes;
	fprintf(stderr, "CPU time for VTK file write = %5.2F seconds\n", elapsed);
	fprintf(stderr, "                          %5.2F million cells / minute\n",
					(totalCells / 1000000.) / (elapsed / 60));

	return true;
}

bool UMesh::writeUGridFile(const char fileName[]) {
	double timeBefore = clock() / double(CLOCKS_PER_SEC);

	FILE* outFile = fopen(fileName, "w");
	if (!outFile) {
		fprintf(stderr, "Couldn't open file %s for writing.  Bummer!\n", fileName);
		return false;
	}

	fwrite(m_fileImage, m_fileImageSize, 1, outFile);
	fclose(outFile);

	double timeAfter = clock() / double(CLOCKS_PER_SEC);
	double elapsed = timeAfter - timeBefore;
	size_t totalCells = size_t(m_nTets) + m_nPyrs + m_nPrisms + m_nHexes;
	fprintf(stderr, "CPU time for UGRID file write = %5.2F seconds\n", elapsed);
	fprintf(stderr, "                          %5.2F million cells / minute\n",
					(totalCells / 1000000.) / (elapsed / 60));

	return true;
}

//bool UMesh::writeCompressedUGridFile(const char fileName[]) {
//	double timeBefore = clock() / double(CLOCKS_PER_SEC);
//
//	gzFile outFile = gzopen(fileName, "wb");
//	if (!outFile) {
//		fprintf(stderr, "Couldn't open file %s for writing.  Bummer!\n", fileName);
//		return false;
//	}
//
//	int bytesWritten = gzwrite(outFile, reinterpret_cast<void*>(m_fileImage),
//															m_fileImageSize);
//	gzclose(outFile);
//	fprintf(stderr, "Wrote %d bytes into a compressed file as %d bytes\n",
//					m_fileImageSize, bytesWritten);
//
//	double timeAfter = clock() / double(CLOCKS_PER_SEC);
//	double elapsed = timeAfter - timeBefore;
//	size_t totalCells = size_t(m_nTets) + m_nPyrs + m_nPrisms + m_nHexes;
//	fprintf(stderr, "CPU time for compressed UGRID file write = %5.2F seconds\n",
//					elapsed);
//	fprintf(stderr, "                          %5.2F million cells / minute\n",
//					(totalCells / 1000000.) / (elapsed / 60));
//
//	return true;
//}

