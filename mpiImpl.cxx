
#include <unistd.h>
#include <cstdio>
#include <boost/mpi.hpp>
#include "ExaMesh.h"
#include "CubicMesh.h"
#include "PARMETIS.h"
#include "UMesh.h"
#include <boost/serialization/unique_ptr.hpp>

std::unique_ptr<UMesh>  ReadMesh ( const char  baseFileName[] , const char type[], 
                                     const char  ugridInfix[]   , const char CGNSFileName[], 
                                     const char MeshType)
{

    // if (MeshType == 'C')
    // {
    //     auto    pInitialMsh = std::make_unique<CubicMesh>(CGNSFileName);
    //     return  pInitialMsh;
    // }
   // if (MeshType == 'U')
   // {
        auto    pInitialMsh = std::make_unique<UMesh>(baseFileName, type, ugridInfix);
        return  pInitialMsh;
   // }
   //return nullptr; 

}



void refineForMPI ( const char  baseFileName[] , const char type[], 
                    const char  ugridInfix[]   , const char CGNSFileName[],
                    const int   numDivs        , const char MeshType, 
                    std::string mshName        , FILE* fileAllTimes)
{
    // boost::mpi::environment   env; 
	//boost::mpi::communicator  world;
    //if(world.rank()==MASTER)
   // {
    std::unique_ptr<UMesh> pEM = ReadMesh(baseFileName,type,ugridInfix,CGNSFileName,MeshType);
    std::cout<<"Done with Reading the mesh "<<std::endl; 
    assert(pEM!=nullptr); 
    emInt nParts = 4; 
    auto part2cell= partitionMetis(pEM,nParts); 


	std::cout<<"Done with Metis Partitioning"<<std::endl; 
   // }
    for (auto ipart=0 ; ipart<nParts ; ipart++)
    {
	    auto UM = pEM->Extract(part2cell[ipart],numDivs); 
        std::cout<<"Done with extraction of part ID of "<<ipart<<std::endl; 
        char filename[100];
        sprintf(filename, "Trash/Testsubmesh%03d.vtk", ipart);
        UM->writeVTKFile(filename);
    } 

    // std::unique_ptr<UMesh> pEM = ReadMesh(baseFileName,type,ugridInfix,CGNSFileName,MeshType); 
    // std::vector<Part> parts;
	// std::vector<CellPartData> vecCPD;
    // emInt nParts = 10;  
	// partitionCells(pEM.get(), nParts, parts, vecCPD);

    // for (auto ipart=0 ; ipart<nParts ; ipart++)
    // {
    //     pEM->extractCoarseMesh(parts[ipart],vecCPD,numDivs); 
    // }

}



