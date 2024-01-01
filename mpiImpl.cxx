
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

    emInt nParts = 4; 

    double start = exaTime(); 
    std::unique_ptr<UMesh> pEM = ReadMesh(baseFileName,type,ugridInfix,CGNSFileName,MeshType);
    double time = exaTime()-start; 
    std::cout<<"Reading mesh took: "<<time<<std::endl; 
    assert(pEM!=nullptr); 
    
    start = exaTime(); 
    auto part2cell= partitionMetis(pEM,nParts); 
    time = exaTime()-start; 
    std::cout<<"Metis Partitioning took: "<<time<<std::endl; 
    vecHashTri  tris; 
	vecHashQuad quads;
    std::size_t triSize; 
    std::size_t quadSize; 
    start = exaTime(); 
    pEM->partFaceMatching(part2cell,tris,quads,triSize,quadSize); 
    time = exaTime()-start; 
    std::cout<<"Part face Matching took: "<<time<<std::endl; 


	// for (auto ipart=0 ; ipart<nParts ; ipart++)
    // {
	//     auto UM = pEM->Extract(ipart, part2cell[ipart],numDivs,tris[ipart],quads[ipart]); 
    //     std::cout<<"Done with extraction of part ID of "<<ipart<<std::endl; 
    //     char filename[100];
    //     sprintf(filename, "Trash/Testsubmesh%03d.vtk", ipart);
    //     UM->writeVTKFile(filename);
    // } 

    // double start = exaTime(); 
    // std::unique_ptr<UMesh> pEM = ReadMesh(baseFileName,type,ugridInfix,CGNSFileName,MeshType); 
    // double time = exaTime()- start; 
    // std::cout<<"Reading mesh took: "<<time<<std::endl; 
    // std::vector<Part> parts;
	// std::vector<CellPartData> vecCPD;
    // start = exaTime(); 
	// partitionCells(pEM.get(), nParts, parts, vecCPD);
    // time  = exaTime()-start; 
    // std::cout<<"Partitioning took: "<<time<<std::endl; 


}



