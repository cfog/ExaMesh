
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

    // emInt nCells       =   pEM->cellID2cellTypeLocalID.size();
    // auto lengthscale   =   pEM->getAllLenghtScale(); 

    // std::vector<double> vLengthScale(lengthscale,lengthscale+ pEM->m_header[0]); 
    // std::vector<emInt>  vheader(pEM->m_header, pEM->m_header+7); 

    // for(auto irank=1 ; irank<world.size();irank++)
    // {
    //     world.send(irank,0,vheader);
    //     world.send(irank,0,nCells); 
    //     world.send(irank,0, pEM-> cellID2cellTypeLocalID);
    //     world.send(irank,0, pEM->vTetConns);
    //     world.send(irank,0, pEM->vTriConns);
    //     world.send(irank,0, vLengthScale); 
    //     world.send(irank,0, pEM->vQuadConns); 
    //     world.send(irank,0, pEM->vPyrmConns); 
    //     world.send(irank,0, pEM->vPrsimConns); 
    //     world.send(irank,0, pEM->vHexConns); 
    // }
     

}

std::unique_ptr<UMesh>
RecvUMesh(boost::mpi::communicator  world)
{
   
    // emInt nCells; 
 

    // std::vector<std::pair<emInt,emInt>> cellId2type;
    // std::vector<emInt>  vheader(7);  
    

    // world.recv(MASTER,0,vheader); 

    // std::vector<std::vector<emInt>> vTriConns  (vheader[1],std::vector<emInt>(3)); 
    // std::vector<std::vector<emInt>> vQuadConn  (vheader[2],std::vector<emInt>(4));
    // std::vector<std::vector<emInt>> vTetConns  (vheader[3],std::vector<emInt>(4)); 
    // std::vector<std::vector<emInt>> vPyrmConn  (vheader[4],std::vector<emInt>(5)); 
    // std::vector<std::vector<emInt>> vPrsimConn (vheader[5],std::vector<emInt>(6)); 
    // std::vector<std::vector<emInt>> vHexConn   (vheader[6],std::vector<emInt>(8)); 

 


    // world.recv(MASTER,0,nCells);

    // cellId2type.resize(nCells); 
    // world.recv(MASTER,0,cellId2type);

    // world.recv(MASTER,0,vTetConns);
    

  
    // world.recv(MASTER,0,vTriConns); 

    // std::vector<double> vLengthScale(vheader[0]);
    // world.recv(MASTER,0,vLengthScale); 
    // world.recv(MASTER,0,vQuadConn); 
    // world.recv(MASTER,0,vPyrmConn);
    // world.recv(MASTER,0,vPrsimConn); 
    // world.recv(MASTER,0,vHexConn); 



    // auto    iniMesh = std::make_unique<UMesh>(vheader[0],0,vheader[1],vheader[2],
    // vheader[3],vheader[4],vheader[5],vheader[6],vTetConns,vheader,cellId2type,vTriConns,vLengthScale,
    // vQuadConn,vPyrmConn,vPrsimConn,vHexConn); 

  

    // return iniMesh; 


}
void sendTrisQuads(boost::mpi::communicator  world, 
    const  vecHashTri  &tris,  
	const  vecHashQuad &quads, 
    const  vecVecTri   &vtris, 
    const  vecVecQuad  &vquads,
    hashTri  &trisS, 
	hashQuad &quadsS)
{
    // Sending Tri and Quad Data
    vecVecTri   VecTriVec; 
	vecVecQuad  vecQuadVec; 
    int tag=0; 


    // for(size_t  itri=0 ; itri<tris.size(); itri++)
	// {
	// 	vecTri TriVec; 
	// 	SetToVector(tris[itri],TriVec); 
	// 	VecTriVec.emplace_back(TriVec); 
	// }
	// for(size_t iquad=0 ; iquad<quads.size(); iquad++)
	// {
	// 	vecQuad QuadVec; 
	// 	SetToVector(quads[iquad],QuadVec); 
	// 	vecQuadVec.emplace_back(QuadVec); 
	// }

	//trisS = tris[0]; // For MASTER 
	//quadsS= quads[0];

    for(auto itri=0; itri<vtris[0].size(); itri++)
    {
        trisS.insert(vtris[0][itri]);
    }
    for(auto iquad=0; iquad<vquads[0].size(); iquad++)
    {
        quadsS.insert(vquads[0][iquad]);
    }
		
	for(auto irank=1 ; irank<world.size();irank++)
	{

		world.send(irank,tag,vtris[irank]); 
		world.send(irank,tag,vquads[irank]); 

	}
}

void printHashTris(const hashTri& hashTris) 
{
    for (const auto& hashTri : hashTris)
    {
        std::cout << "Global ID: " << hashTri.getGlobalCorner(0)<<" "<< hashTri.getGlobalCorner(1)<<" "<< hashTri.getGlobalCorner(2) << std::endl;
        std::cout << "Sorted Global ID: " << hashTri.getSortedGlobal(0)<<" "<< hashTri.getSortedGlobal(1)<<" "<< hashTri.getSortedGlobal(2) << std::endl;
        std::cout << "Part ID: " <<hashTri.getPartid() << std::endl;
        std::cout << "Remote Part ID: " << hashTri.getRemoteId()<< std::endl;
    }
}
void printHashQuads(const hashQuad& hashQuads) 
{
    for (const auto& hashQuad : hashQuads)
    {
        std::cout << "Global ID: " << hashQuad.getGlobalCorner(0)<<" "<< hashQuad.getGlobalCorner(1)<<" "<< hashQuad.getGlobalCorner(2)<<" "<< hashQuad.getGlobalCorner(3) << std::endl;
        std::cout << "Sorted Global ID: " << hashQuad.getSortedGlobal(0)<<" "<< hashQuad.getSortedGlobal(1)<<" "<< hashQuad.getSortedGlobal(2)<<" "<< hashQuad.getSortedGlobal(3) << std::endl;
        std::cout << "Part ID: " <<hashQuad.getPartid() << std::endl;
        std::cout << "Remote Part ID: " << hashQuad.getRemoteId()<< std::endl;
    }
}

void refineForMPI ( const char  baseFileName[] , const char type[], 
                    const char  ugridInfix[]   , const char CGNSFileName[],
                    const int   numDivs        , const char MeshType, 
                    std::string mshName        , FILE* fileAllTimes)
{



    double appTimeStart =exaTime();
    boost::mpi::environment   env; 
	boost::mpi::communicator  world;
   

    emInt nParts = world.size(); 
    std::vector<boost::mpi::request> reqs;
    std::vector<boost::mpi::request> sizes; 

    std::vector<boost::mpi::request> reqforPartCellSizes;
    std::vector<boost::mpi::request> temp;


    std::vector<std::size_t>  partCellsizes(nParts);
    std::vector<std::size_t>  trisizes(nParts); 
    std::vector<std::size_t>  quadsize(nParts);
    std::vector<std::vector<emInt>> partCells(nParts); 
    


    hashTri  hashTris; 
	hashQuad hashQuads;

    vecVecTri   triV(nParts);
	vecVecQuad  quadV(nParts); 

    double starttimereading=exaTime();
   
    UMesh inimesh(baseFileName, type, ugridInfix); 
    double timereading=exaTime()-starttimereading;
   
    std::cout<<"Reading mesh on rank: "<<world.rank()<<" took: "<<timereading<<std::endl;
   
    
    
    if(world.rank()==MASTER)
    {
        
        std::vector<emInt> vaicelltopart;
        auto part2cell= partitionMetis(inimesh,nParts,vaicelltopart); 
        partCells[world.rank()] = part2cell[0];

 
        for(auto irank=1; irank<world.size();irank++)
        {
            boost::mpi::request rq0=world.isend(irank,0,part2cell[irank].size());
            reqforPartCellSizes.emplace_back(rq0);

            boost::mpi::request rq1 =world.isend(irank,1,part2cell[irank]); 
            reqs.emplace_back(rq1);
        }
        vecHashTri  hashtris; 
	    vecHashQuad hashquads;

        vecVecTri  tris;
        vecVecQuad quads;

        inimesh.FastpartFaceMatching(nParts, part2cell,vaicelltopart,tris,quads); 
        for(auto irank=1 ; irank<world.size();irank++)
	    {


            boost::mpi::request rq1= world.isend(irank,2,tris[irank].size());
            sizes.emplace_back(rq1);

            boost::mpi::request rq2= world.isend(irank,3,quads[irank].size());
            sizes.emplace_back(rq2);

            boost::mpi::request rq3= world.isend(irank,4,tris[irank]); 
            temp.emplace_back(rq3);

            boost::mpi::request rq4= world.isend(irank,5,quads[irank]);
            temp.emplace_back(rq4);

        }
        vectorToSet(tris[MASTER],hashTris);
        vectorToSet(quads[MASTER],hashQuads);
        inimesh.Extract(world.rank(),partCells[MASTER],numDivs,hashTris,hashQuads);

      
    }
    if(world.rank()!=MASTER)
    {
       
        boost::mpi::request rq0=world.irecv(MASTER,0,partCellsizes[world.rank()]);
        reqforPartCellSizes.emplace_back(rq0);

        boost::mpi::wait_all(reqforPartCellSizes.begin(),reqforPartCellSizes.end());
        partCells[world.rank()].resize(partCellsizes[world.rank()]);


        boost::mpi::request req=world.irecv(MASTER,1,partCells[world.rank()]); 
        reqs.emplace_back(req);


        

        boost::mpi::request rq1= world.irecv(MASTER,2,trisizes[world.rank()]);
        sizes.emplace_back(rq1);
        boost::mpi::request rq3= world.irecv(MASTER,3,quadsize[world.rank()]);
        sizes.emplace_back(rq3);
        boost::mpi::wait_all(sizes.begin(),sizes.end());
        


        triV[world.rank()].resize(trisizes[world.rank()]);
        quadV[world.rank()].resize(quadsize[world.rank()]);

        boost::mpi::request rq2= world.irecv(MASTER,4,triV[world.rank()]);
        temp.emplace_back(rq2);

        boost::mpi::request rq4= world.irecv(MASTER,5,quadV[world.rank()]);
        temp.emplace_back(rq4);

    }
    boost::mpi::wait_all(reqs.begin(),reqs.end());
    boost::mpi::wait_all(temp.begin(),temp.end());
  
  
   
   
    if(world.rank()!=MASTER)
    {
       vectorToSet(triV[world.rank()],hashTris);
       vectorToSet(quadV[world.rank()],hashQuads); 
       inimesh.Extract(world.rank(),partCells[world.rank()],numDivs,hashTris,hashQuads);
    }
    double time=exaTime()-appTimeStart;
    std::cout<<"My rank: "<<world.rank()<<" My total time: "<<time<<std::endl;


  
   

}



