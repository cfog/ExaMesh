
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

void printCellIDs(std::vector<std::pair<emInt,emInt>> &vec)
{
    for(auto i=0 ; i<vec.size(); i++)
    {
        std::cout<<"Cell ID: "<<i <<" First: "<<vec[i].first<<" Second: "<<vec[i].second<<std::endl;
    }
}
void printTetConn(const emInt (*tetconn)[4], const emInt nTets)
{
    for(auto iTet=0 ; iTet<nTets ; iTet++)
    {
        for(auto j=0; j<4; j++)
        {
            std::cout<<tetconn[iTet][j]<<std::endl; 
        }  
    }
}

void sendUMesh(boost::mpi::communicator  world, UMesh* pEM)
{

    emInt nCells       =   pEM->cellID2cellTypeLocalID.size();
    auto lengthscale   =   pEM->getAllLenghtScale(); 

    std::vector<double> vLengthScale(lengthscale,lengthscale+ pEM->m_header[0]); 
    std::vector<emInt>  vheader(pEM->m_header, pEM->m_header+7); 

    for(auto irank=1 ; irank<world.size();irank++)
    {
        world.send(irank,0,vheader);
        world.send(irank,0,nCells); 
        world.send(irank,0, pEM-> cellID2cellTypeLocalID);
        world.send(irank,0, pEM->vTetConns);
        world.send(irank,0, pEM->vTriConns);
        world.send(irank,0, vLengthScale); 
        world.send(irank,0, pEM->vQuadConns); 
        world.send(irank,0, pEM->vPyrmConns); 
        world.send(irank,0, pEM->vPrsimConns); 
        world.send(irank,0, pEM->vHexConns); 
    }
     

}

std::unique_ptr<UMesh>
RecvUMesh(boost::mpi::communicator  world)
{
   
    emInt nCells; 
 

    std::vector<std::pair<emInt,emInt>> cellId2type;
    std::vector<emInt>  vheader(7);  
    

    world.recv(MASTER,0,vheader); 

    std::vector<std::vector<emInt>> vTriConns  (vheader[1],std::vector<emInt>(3)); 
    std::vector<std::vector<emInt>> vQuadConn  (vheader[2],std::vector<emInt>(4));
    std::vector<std::vector<emInt>> vTetConns  (vheader[3],std::vector<emInt>(4)); 
    std::vector<std::vector<emInt>> vPyrmConn  (vheader[4],std::vector<emInt>(5)); 
    std::vector<std::vector<emInt>> vPrsimConn (vheader[5],std::vector<emInt>(6)); 
    std::vector<std::vector<emInt>> vHexConn   (vheader[6],std::vector<emInt>(8)); 

 


    world.recv(MASTER,0,nCells);

    cellId2type.resize(nCells); 
    world.recv(MASTER,0,cellId2type);

    world.recv(MASTER,0,vTetConns);
    

  
    world.recv(MASTER,0,vTriConns); 

    std::vector<double> vLengthScale(vheader[0]);
    world.recv(MASTER,0,vLengthScale); 
    world.recv(MASTER,0,vQuadConn); 
    world.recv(MASTER,0,vPyrmConn);
    world.recv(MASTER,0,vPrsimConn); 
    world.recv(MASTER,0,vHexConn); 



    auto    iniMesh = std::make_unique<UMesh>(vheader[0],0,vheader[1],vheader[2],
    vheader[3],vheader[4],vheader[5],vheader[6],vTetConns,vheader,cellId2type,vTriConns,vLengthScale,
    vQuadConn,vPyrmConn,vPrsimConn,vHexConn); 

  

    return iniMesh; 


}
void sendTrisQuads(boost::mpi::communicator  world, 
    const  vecHashTri  &tris,  
	const  vecHashQuad &quads, 
    hashTri  &trisS, 
	hashQuad &quadsS)
{
    // Sending Tri and Quad Data
    vecVecTri   VecTriVec; 
	vecVecQuad  vecQuadVec; 
    int tag=0; 


    for(size_t  itri=0 ; itri<tris.size(); itri++)
	{
		vecTri TriVec; 
		SetToVector(tris[itri],TriVec); 
		VecTriVec.emplace_back(TriVec); 
	}
	for(size_t iquad=0 ; iquad<quads.size(); iquad++)
	{
		vecQuad QuadVec; 
		SetToVector(quads[iquad],QuadVec); 
		vecQuadVec.emplace_back(QuadVec); 
	}

	trisS = tris[0]; // For MASTER 
	quadsS= quads[0];
		
	for(auto irank=1 ; irank<world.size();irank++)
	{

		world.send(irank,tag,VecTriVec[irank]); 
		world.send(irank,tag,vecQuadVec[irank]); 

	}
}



void refineForMPI ( const char  baseFileName[] , const char type[], 
                    const char  ugridInfix[]   , const char CGNSFileName[],
                    const int   numDivs        , const char MeshType, 
                    std::string mshName        , FILE* fileAllTimes)
{
    //boost::mpi::environment   env; 
	//boost::mpi::communicator  world;
    //double appTimeStart =exaTime(); 

    emInt nParts = world.size(); 

   // std::unique_ptr<UMesh> pEM(nullptr); 

    std::vector<emInt> partCells; 

    hashTri  trisS; 
	hashQuad quadsS;

    vecTri   triV;
	vecQuad  quadV; 

    int tag=0; 
   

    //if(world.rank()==MASTER)
   // {
        double start = exaTime();
        UMesh inimesh(baseFileName, type, ugridInfix); 
        double time = exaTime()-start;
        

        std::cout<<"Reading mesh & building connectivities in total on MASTER took: "<<time<<std::endl; 
        //auto celltypes= inimesh.getCellID2CellType2LocalID();
        //emInt celltypesSize =celltypes.size(); 
        //start=exaTime();
       // for(auto irank=1 ; irank<world.size(); irank++)
       // {
            //world.send(irank,0,inimesh); 
            //world.send(irank,0,celltypesSize);
           // world.send(irank,0,celltypes);
        //}
        //time=exaTime()-start; 
        //std::cout<<"Sending Umesh took: "<<time<<std::endl; 


        start = exaTime(); 
        std::vector<emInt> vaicelltopart;
        auto part2cell= partitionMetis(inimesh,nParts,vaicelltopart); 
       // partCells = part2cell[0];
        time = exaTime()-start; 
        std::cout<<"Metis Partitioning on MASTER took: "<<time<<std::endl;

        //for(auto irank=1; irank<world.size();irank++)
        //{
         //   world.send(irank,0,part2cell[irank]); 
        //}
        vecHashTri  hashtris; 
	    vecHashQuad hashquads;

        vecVecTri  tris;
        vecVecQuad quads;

        std::size_t triSize; 
        std::size_t quadSize; 
 
        start = exaTime(); 
        inimesh.partFaceMatching(part2cell,hashtris,hashquads,triSize,quadSize); 
        time = exaTime()-start; 
        std::cout<<"Old Part face Matching on MASTER took: "<<time<<std::endl; 

        start = exaTime();
        inimesh.FastpartFaceMatching(nParts, part2cell,vaicelltopart,tris,quads); 
        time = exaTime()-start;
        std::cout<<"New Part face Matching on MASTER took: "<<time<<std::endl;
        TestPartFaceMatching(nParts,hashtris,hashquads,tris,quads);



        //start=exaTime(); 
        //sendTrisQuads(world,tris,quads,trisS,quadsS);
        //time=exaTime()-start; 
        //std::cout<<"Sending Boundary Faces took: "<<time<<std::endl; 

    //}
    //if(world.rank()!=MASTER)
   // {
    //     std::vector<std::pair<emInt,emInt>> cellId2type;
    //     emInt celltypesSize; 
       // UMesh recvMesh; 
        //world.recv(MASTER,0,recvMesh);
    //     //world.recv(MASTER,0,celltypesSize); 
    //    // cellId2type.resize(celltypesSize);
    //    // world.recv(MASTER,0,cellId2type); 
    //     world.recv(MASTER,0,partCells); 


    //     recvMesh.setCellId2CellTypeLocal(cellId2type); 

    //     world.recv(MASTER,tag,triV); 
	//  	world.recv(MASTER,tag,quadV);

    //     vectorToSet(triV,trisS);
	//  	vectorToSet(quadV,quadsS); 

    //     double start = exaTime(); 
    //     recvMesh.convertToUmeshFormat();
    //     double time = exaTime()-start; 
    //     std::cout<<"Converting Umesh pack: "<<time<<std::endl;
       // recvMesh.Extract(world.rank(),partCells,numDivs,trisS,quadsS); 
        
        
        
   // }

    //world.barrier(); 
    //double appTime = exaTime()-appTimeStart; 
    //std::cout<<"My rank: "<<world.rank()<<" My total time: "<<appTime<<std::endl; 

  


    // if(world.rank()==MASTER)
    // {
    //     double start = exaTime(); 
    //     pEM = ReadMesh(baseFileName,type,ugridInfix,CGNSFileName,MeshType);
        
    //     double time = exaTime()-start; 
    //     std::cout<<"Reading mesh & building connectivities in total on MASTER took: "<<time<<std::endl; 
    //     assert(pEM!=nullptr); 


    //     start=exaTime(); 
    //     sendUMesh(world,pEM.get()); 
    //     time= exaTime()-start; 
    //     std::cout<<"Packing UMesh Data and sending it to others took: "<<time<<std::endl; 
        
    //     start = exaTime(); 
    //     auto part2cell= partitionMetis(pEM,nParts); 
    //     partCells = part2cell[0];
    //     time = exaTime()-start; 
    //     std::cout<<"Metis Partitioning on MASTER took: "<<time<<std::endl; 

    //     vecHashTri  tris; 
	//     vecHashQuad quads;

    //     std::size_t triSize; 
    //     std::size_t quadSize; 
    //     for(auto irank=1; irank<world.size();irank++)
    //     {
    //         world.send(irank,0,part2cell[irank]); 
    //     }

    //     start = exaTime(); 
    //     pEM->partFaceMatching(part2cell,tris,quads,triSize,quadSize); 
    //     time = exaTime()-start; 
    //     std::cout<<"Part face Matching on MASTER took: "<<time<<std::endl; 


    //     start=exaTime(); 
    //     sendTrisQuads(world,tris,quads,trisS,quadsS);
    //     time=exaTime()-start; 
    //     std::cout<<"Sending Boundary Faces took: "<<time<<std::endl; 

    //     pEM->Extract(world.rank(),partCells,numDivs,trisS,quadsS); 
    //     int a=3 ; 
    //     world.isend(1,0,a);
        
        

    // }

    // if(world.rank()!=MASTER)
    // {
 
    //     pEM= RecvUMesh(world);
    //     assert(pEM!=nullptr);
    //     world.recv(MASTER,0,partCells); 

    //     world.recv(MASTER,tag,triV); 
	// 	world.recv(MASTER,tag,quadV);


    //     hashTri trisHash;

	// 	vectorToSet(triV,trisS);
	// 	vectorToSet(quadV,quadsS); 


    //    pEM->Extract(world.rank(),partCells,numDivs,trisS,quadsS);   
   
 
    // }

 
    // world.barrier(); 
    // double appTime = exaTime()-appTimeStart; 
    // std::cout<<"My rank: "<<world.rank()<<" My total time: "<<appTime<<std::endl; 

}



