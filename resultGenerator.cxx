#include <iostream>
#include <cstdio>
#include <inttypes.h>
#include "exa-defs.h"

void inline writeAllTimeResults(FILE *file, int nP, double Tpartition,
		double TpartFaceMatching, double Tserial, double Textraction,
		double Trefinement, double TfaceExchange, double syncTime,
		double TtriMatch, double TquadMatch, double Ttotal, size_t triSize,
		size_t quadSize, size_t nCells) {
	// setlocale(LC_ALL, "");
	fseek(file, 0, SEEK_END);
	long size = ftell(file);
	if (size == 0) {
		fprintf(file,
				"%-5s %-12s %-18s %-12s %-12s %-12s %-18s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n",
				"nP", "Tpartition", "TpartFaceMatch", "Tserial", "Textract",
				"Trefine", "TfaceExchange", "Tsync", "TtriMatch", "TquadMatch",
				"Ttoal", "triSize", "quadSize", "nCells");
	}

	// fprintf(file, "%-5u %-12f %-18f %-12f %-12f %-12f %-12f %-12f %-12" PRId64"\n",
	// nP,Tpartition,TpartFaceMatching,Textraction,Trefinement,
	// TtriMatch,TquadMatch,Ttotal,nCells);
	// if(sizeof(nCells)==sizeof(int32_t))
	// {

	fprintf(file,
			"%-5u %-12f %-18f %-12f %-12f %-12f %-18f %-12f %-12f %-12f %-12f  %-12lu %-12lu %'zd\n",
			nP, Tpartition, TpartFaceMatching, Tserial, Textraction,
			Trefinement, TfaceExchange, syncTime, TtriMatch, TquadMatch, Ttotal,
			triSize, quadSize, nCells);
	//}
	// if(sizeof(nCells)==sizeof(int64_t))
	// {
	//     fprintf(file, "%-5u %-12f %-18f %-12f %-12f %-12f %-12f %-12f %-12" PRId64"\n",
	//     nP,Tpartition,TpartFaceMatching,Textraction,Trefinement,
	//     TtriMatch,TquadMatch,Ttotal,nCells);

	// }
}

void inline writeMeshStatics(FILE *file, int ndivs, emInt nCells) {
	fseek(file, 0, SEEK_END);
	long size = ftell(file);
	if (size == 0) {
		// if(sizeof(nCells)==sizeof(int64_t))
		//     fprintf(file, "%-12s %-12" PRId64"\n","nDivs", "nCells");
		if (sizeof(nCells) == sizeof(int32_t))
			// fprintf(file, "%-12s %-12" PRId32"\n","nDivs", "nCells");
			fprintf(file, "%-12s %-12s \n", "nDivs", "nCells");
	}
	fprintf(file, "%-12u %-12u\n", ndivs, nCells);
}

void inline writeEachRankMeshStatics(int rank, int nCells, int nDivs, int nProc,
		std::string fileName) {
	fileName = fileName + "-nDivs-" + std::to_string(nDivs) + "-nProc"
			+ std::to_string(nProc) + ".txt";

	FILE *file = fopen(fileName.c_str(), "a");

	if (file == NULL) {
		fprintf(stderr, "Error opening the file!\n");
	}

	fseek(file, 0, SEEK_END);
	long size = ftell(file);
	if (size == 0) {
		fprintf(file, "%-12s %-12s\n", "Rank", "nRefinedCells");
	}

	fprintf(file, "%-12u %-12u\n", rank, nCells);
}
void inline printTimeEachRank(FILE *file, emInt rank, double partitioning,
		double partfaceMatching, double serialTime, double extraction,
		double refinment, double faceExchange, double totalSyncTime,
		double triMatch, double quadMatch, double total) {
	fseek(file, 0, SEEK_END);
	long size = ftell(file);
	if (size == 0) {
		fprintf(file,
				"%-5s %-12s %-18s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s \n",
				"Rank", "Tpartition", "TpartFaceMatch", "Serial Time",
				"Textract", "Trefine", "TfaceExchange", "TotalSyncTime",
				"TtriMatch", "TquadMatch", "Ttoal");
	}

	fprintf(file,
			"%-5u %-12f %-18f %-12f %-12f %-12f %-12f %-12f %-12f %-12f %-12f\n",
			rank, partitioning, partfaceMatching, serialTime, extraction,
			refinment, faceExchange, totalSyncTime, triMatch, quadMatch, total);
}

void inline writeEachRankMeshStatics(FILE *file, emInt rank, int tet, int pyrm,
		int prism, int hex, int nBdryTris, int nBdryQuads, int total) {

	setlocale(LC_NUMERIC, "");
	fseek(file, 0, SEEK_END);
	long size = ftell(file);
	if (size == 0) {
		fprintf(file, "%-10s %-12s %-12s %-12s %-8s %-40s %-40s %-40s\n",
				"Rank", "tet", "pyrm", "prism", "hex", "nprtBdryTris",
				"nprtBdryQuads", "total");
	}
	setlocale(LC_NUMERIC, "");
	fprintf(file, "%-10u %-12u %-12u %-12u %-8u %-40u %-40u %-40u\n", rank, tet,
			pyrm, prism, hex, nBdryTris, nBdryQuads, total);
}
