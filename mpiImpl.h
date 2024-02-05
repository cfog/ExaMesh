#include <unistd.h>
#include <cstdio>
#include "ExaMesh.h"
#include "CubicMesh.h"
#include "PARMETIS.h"
#include "UMesh.h"


std::unique_ptr<UMesh>  ReadMesh ( const char  baseFileName[] , const char type[], 
                                     const char  ugridInfix[]   , const char CGNSFileName[], 
                                     const char MeshType); 


void refineForMPI ( const char  baseFileName[] , const char type[], 
                    const char  ugridInfix[]   , const char CGNSFileName[],
                    const int   numDivs        , const char MeshType, 
                    std::string mshName, FILE* eachRank);



