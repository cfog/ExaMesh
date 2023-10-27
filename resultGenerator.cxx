#include <iostream>
#include <cstdio>
#include <inttypes.h>
#include "exa-defs.h"

void 
inline writeAllTimeResults (FILE* file, int nP, double Tpartition, 
double TpartFaceMatching, double Textraction, double Trefinement, 
double TtriMatch, double TquadMatch, double Ttotal ,emInt nCells)
{
    //setlocale(LC_ALL, "");
    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    if (size == 0) 
	{
	 	fprintf(file, "%-5s %-12s %-18s %-12s %-12s %-12s %-12s %-12s %-12s\n",
          "nP", "Tpartition", "TpartFaceMatch", 
         "Textract","Trefine","TtriMatch","TquadMatch","Ttoal","nCells");
   	}
	
    // fprintf(file, "%-5u %-12f %-18f %-12f %-12f %-12f %-12f %-12f %-12" PRId64"\n",
    // nP,Tpartition,TpartFaceMatching,Textraction,Trefinement, 
    // TtriMatch,TquadMatch,Ttotal,nCells);
    if(sizeof(nCells)==sizeof(int32_t))
    {
        fprintf(file, "%-5u %-12f %-18f %-12f %-12f %-12f %-12f %-12f %-12" PRId32"\n",
        nP,Tpartition,TpartFaceMatching,Textraction,Trefinement, 
        TtriMatch,TquadMatch,Ttotal,nCells);
    }
    // if(sizeof(nCells)==sizeof(int64_t))
    // {
    //     fprintf(file, "%-5u %-12f %-18f %-12f %-12f %-12f %-12f %-12f %-12" PRId64"\n",
    //     nP,Tpartition,TpartFaceMatching,Textraction,Trefinement, 
    //     TtriMatch,TquadMatch,Ttotal,nCells);

    // }

}

void 
inline writeMeshStatics(FILE* file, int ndivs ,emInt nCells)
{
    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    if (size == 0) 
	{
        //if(sizeof(nCells)==sizeof(int64_t))
        //    fprintf(file, "%-12s %-12" PRId64"\n","nDivs", "nCells");
        if(sizeof(nCells)==sizeof(int32_t))
            //fprintf(file, "%-12s %-12" PRId32"\n","nDivs", "nCells");
            fprintf(file, "%-12s %-12s \n","nDivs", "nCells");
    }
	fprintf(file, "%-12u %-12u\n",
    ndivs,nCells);
}

void
inline  
writeEachRankMeshStatics (int rank, int nCells, int nDivs, int nProc ,std::string fileName)
{   
    fileName = fileName+"-nDivs-"+std::to_string(nDivs)+"-nProc"+std::to_string(nProc)+".txt";

    FILE *file = fopen(fileName.c_str(), "a");

    if (file == NULL) 
	{
        fprintf(stderr, "Error opening the file!\n");
    }

    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    if (size == 0) 
	{
	 	fprintf(file, "%-12s %-12s\n",
          "Rank", "nRefinedCells");
   	}
	
    fprintf(file, "%-12u %-12u\n",
    rank,nCells);

}
