
#include <unistd.h>
#include <cstdio>
#include <boost/mpi.hpp>
#include "ExaMesh.h"
#include "CubicMesh.h"
#include "PARMETIS.h"
#include "UMesh.h"
#include <boost/serialization/unique_ptr.hpp>

void inline NewWriteTimes(
    FILE *file, int nP, const timeResults &times,
    size_t nCells)
{
    setlocale(LC_ALL, "");
    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    if (size == 0)
    {
        fprintf(file, "%-5s %-10s %-14s %-15s %-12s %-12s %-14s %-14s %-14s %-16s %-10s %-14s\n",
                "nP", "Read", "Partition", "PartFaceMatch", "Serial",
                "Extract", "Refine", "FaceExchange",
                "TotalSync", "TotalMatchTime", "Total",
                "nCells");
    }

    fprintf(file, "%-5u %-10.2f %-14.2f %-15.2f %-12.2f %-12.2f %-14.2f %-14.2f %-14.2f %-16.2f %-10.2f %'-zd\n",
            nP, times.read, times.partition, times.partfacematching, times.serial,
            times.extract, times.refine, times.faceExchange,
            times.syncTri+times.syncQuad, times.matchtris+times.matchquads, times.total, nCells);



}
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

void printPartCells (const std::vector<emInt>& partCells, const emInt myRank)
{
    std::cout<<"My rank: "<<myRank<<std::endl;
    for(auto irank=0 ; irank<10; irank++)
    {
        std::cout<<partCells[irank]<<std::endl;
    }
}

void refineForMPI ( const char  baseFileName[] , const char type[], 
                    const char  ugridInfix[]   , const char CGNSFileName[],
                    const int   numDivs        , const char MeshType, std::string mshName 
                    )
{

    size_t lastSlashPos        = std::string(baseFileName).find_last_of('/');
	mshName                    = std::string(baseFileName).substr(lastSlashPos + 1)+"_WeakScaleTimes.txt";
    std::string totalTimeFile  = std::string(baseFileName).substr(lastSlashPos + 1)+"_WeakOnlyTotalTime.txt";

    timeResults times, maxTimes;
    times.total=exaTime();
    times.serial=exaTime();
   
    boost::mpi::environment   env; 
	boost::mpi::communicator  world;
    emInt nParts = world.size(); 

    std::vector<boost::mpi::request> reqForPartCells,reqForCoarseFaces;
    std::vector<boost::mpi::request> triReqs, quadReqs;
    std::vector<std::vector<emInt>>   partCells; 

   
  

    hashTri  hashTris; 
	hashQuad hashQuads;

    vecVecTri   triV(nParts);
	vecVecQuad  quadV(nParts);

    intToVecTri  remoteTovecTris;
	intToVecQuad remoteTovecQuads; 

	std::set<int> triNeighbrs;  
	std::set<int> quadNeighbrs; 

    std::vector<vecTri>  trisTobeRcvd; 
	std::vector<vecQuad> quadsTobeRcvd;

    hashTri  recvdTris;
	hashQuad recvdQuads; 

	TableTri2TableIndex2Index   matchedTris; 
	TableQuad2TableIndex2Index  matchedQuads; 

    std::vector<std::size_t>  trisizes(nParts); 
    std::vector<std::size_t>  quadsize(nParts);
    
    times.read=exaTime(); 
	std::unique_ptr<UMesh> inimesh = 
    std::make_unique<UMesh>(baseFileName, type, ugridInfix);
    times.read=exaTime()-times.read;


    if(world.rank()==MASTER)
    {
        times.partition=exaTime();
        std::vector<emInt> vaicelltopart;
        auto part2cell= partitionMetis(inimesh.get(),nParts,vaicelltopart); 
        partCells = part2cell;
        times.partition=exaTime()-times.partition;

        vecHashTri  hashtris; 
	    vecHashQuad hashquads;

        vecVecTri  tris;
        vecVecQuad quads;

        times.partfacematching=exaTime();
        inimesh->FastpartFaceMatching(nParts, part2cell,vaicelltopart,tris,quads);
        times.partfacematching=exaTime()-times.partfacematching; 
       

        for(auto irank=1 ; irank<world.size();irank++)
	    {
            boost::mpi::request rq0 =world.isend(irank,1,part2cell); 
            reqForPartCells.emplace_back(rq0);

            boost::mpi::request rq1= world.isend(irank,2,tris[irank]); 
            reqForCoarseFaces.emplace_back(rq1);

            boost::mpi::request rq2= world.isend(irank,3,quads[irank]);
            reqForCoarseFaces.emplace_back(rq2);

        }
        vectorToSet(tris[MASTER],hashTris);
        vectorToSet(quads[MASTER],hashQuads);

    }
    if(world.rank()!=MASTER)
    {
        boost::mpi::request rq0=world.irecv(MASTER,1,partCells);
        reqForPartCells.emplace_back(rq0);

        boost::mpi::request rq1= world.irecv(MASTER,2,triV[world.rank()]);
        reqForCoarseFaces.emplace_back(rq1);

        boost::mpi::request rq2= world.irecv(MASTER,3,quadV[world.rank()]);
        reqForCoarseFaces.emplace_back(rq2);
    }

    boost::mpi::wait_all(reqForPartCells.begin(),reqForPartCells.end());
    boost::mpi::wait_all(reqForCoarseFaces.begin(),reqForCoarseFaces.end());

    if(world.rank()!=MASTER)
    {
         vectorToSet(triV[world.rank()],hashTris);
         vectorToSet(quadV[world.rank()],hashQuads); 
    }
    times.serial=exaTime()-times.serial;
    times.extract=exaTime();
    std::unique_ptr<UMesh>  extractedMsh =
    inimesh->Extract(world.rank(),partCells[world.rank()],
    numDivs,hashTris,hashQuads);
    times.extract=exaTime()-times.extract;

    times.refine= exaTime();
    auto refinedMsh = std::make_unique<UMesh>(*((extractedMsh.get())),numDivs,world.rank());
    times.refine= exaTime()-times.refine;
    
    times.faceExchange=exaTime();
    auto refinedTris = refinedMsh->getRefinedPartTris();

    buildTrisMap (refinedTris ,remoteTovecTris ,triNeighbrs);

    for(auto tri:remoteTovecTris)
	{
		boost::mpi::request req= world.isend(tri.first,6,tri.second);
		triReqs.push_back(req);
	}
    
    auto refinedQuads= refinedMsh->getRefinedPartQuads();

    buildQuadsMap(refinedQuads,remoteTovecQuads,quadNeighbrs); 
    for(auto quad:remoteTovecQuads)
	{
		boost::mpi::request req = world.isend(quad.first,7,quad.second);
		quadReqs.push_back(req);  
	}
    
	quadsTobeRcvd.resize(quadNeighbrs.size());

    int Quadjj=0; 

	for(auto isource:quadNeighbrs)
	{
		auto findSource = remoteTovecQuads.find(isource); 
		// same amount of message will be sent and received from other prcoeszsor
		// Be cuatious if this assumption later on will be break
		quadsTobeRcvd[Quadjj].resize(findSource->second.size()); 
		Quadjj++; 
	}
    int Quadkk=0 ; 
	for(auto source: quadNeighbrs)
	{
		boost::mpi::request req = world.irecv(source,7,quadsTobeRcvd[Quadkk]);
		quadReqs.push_back(req); 
		Quadkk++; 
	}
    
    trisTobeRcvd.resize (triNeighbrs.size()); 

    int Trijj=0 ; 
	for(auto isource:triNeighbrs)
	{
		auto findSource = remoteTovecTris.find(isource); 
		// same amount of message will be sent and received from other prcoeszsor
		// Be cuatious if this assumption later on will be break
		trisTobeRcvd[Trijj].resize(findSource->second.size()); 
		Trijj++; 
	}

    int Trikk=0 ; 
	for(auto source: triNeighbrs)
	{
		
		boost::mpi::request req= world.irecv(source,6,trisTobeRcvd[Trikk]);
		triReqs.push_back(req);
		Trikk++; 
	}
    times.faceExchange=exaTime()-times.faceExchange;


    times.syncTri=exaTime();
    boost::mpi::wait_all(triReqs.begin(),triReqs.end()); 
    times.syncTri=exaTime()-times.syncTri;
   
    times.matchtris=exaTime();
    for(const auto& tri: trisTobeRcvd)
	{
		// I'm collecting the whole data into a set 
		// I should have received the whole data from it 
		
		recvdTris.insert(tri.begin(),tri.end()); 
	}
    for(auto it=refinedTris.begin(); it!=refinedTris.end(); it++)
	{
		std::unordered_map<emInt, emInt> localRemote;
		//int rotation = getTriRotation(*it,recvdTris,numDivs);
		//matchTri(*it,rotation,numDivs,recvdTris,localRemote); 
		findRotationAndMatchTris(*it,numDivs,recvdTris,localRemote); 
	}
    times.matchtris=exaTime()-times.matchtris;


    times.syncQuad=exaTime();
    boost::mpi::wait_all(quadReqs.begin(),quadReqs.end()); 
    times.syncQuad=exaTime()-times.syncQuad;
   
    times.matchquads=exaTime();
    for(const auto& quad:quadsTobeRcvd)
	{
		recvdQuads.insert(quad.begin(),quad.end()); 
	}

	for(auto it=refinedQuads.begin(); it!=refinedQuads.end();it++)
	{
		std::unordered_map<emInt,emInt> localRemote; 
		//int rotation = getQuadRotation(*it,recvdQuads,numDivs); 
		//matchQuad(*it,rotation,numDivs,recvdQuads,localRemote);
		findRotationAndMatchQuads(*it,numDivs,recvdQuads,localRemote); 

	}
    times.matchquads=exaTime()-times.matchquads;

    times.total=exaTime()-times.total;
   
    size_t         nCells = refinedMsh->numPyramids()+
    refinedMsh->numPrisms()+refinedMsh->numHexes()+refinedMsh->numTets()+ 
    refinedMsh->numBdryTris()+refinedMsh->numBdryQuads();

	size_t         totalCells;

    boost::mpi::reduce(world, times.read            , maxTimes.read              ,boost:: mpi::maximum<double>(), MASTER);

	boost::mpi::reduce(world, times.serial          , maxTimes.serial            ,boost:: mpi::maximum<double>(), MASTER);

	boost::mpi::reduce(world, times.extract         , maxTimes.extract           ,boost:: mpi::maximum<double>(), MASTER);

	boost::mpi::reduce(world,times.refine           , maxTimes.refine            ,boost:: mpi::maximum<double>(), MASTER);

	boost::mpi::reduce(world, times.faceExchange    , maxTimes.faceExchange      ,boost:: mpi::maximum<double>(), MASTER);

    boost::mpi::reduce(world, times.syncTri         , maxTimes.syncTri           ,boost:: mpi::maximum<double>(), MASTER);

    boost::mpi::reduce(world, times.syncQuad        , maxTimes.syncQuad          ,boost:: mpi::maximum<double>(), MASTER);

    boost::mpi::reduce(world, times.matchtris       , maxTimes.matchtris         ,boost:: mpi::maximum<double>(), MASTER);

    boost::mpi::reduce(world, times.matchquads      , maxTimes.matchquads        ,boost:: mpi::maximum<double>(), MASTER);

    boost::mpi::reduce(world, times.total           , maxTimes.total             ,boost:: mpi::maximum<double>(), MASTER);

    boost::mpi::reduce(world, size_t(nCells)        , totalCells                 ,std::plus<size_t>()          , MASTER);


    if(world.rank()==MASTER)
    {
        FILE *fileWeakScaliblity = fopen(mshName.c_str(), "a");
    	if (fileWeakScaliblity == NULL) 
		{
        	fprintf(stderr, "Error opening the file!\n");
   		}
        maxTimes.partfacematching = times.partfacematching;
        maxTimes.partition = times.partition;
        NewWriteTimes(fileWeakScaliblity,nParts,maxTimes,totalCells);

        FILE *fileTotalTime = fopen(totalTimeFile.c_str(), "a");
        if (fileTotalTime == NULL)
        {
            fprintf(stderr, "Error opening the file!\n");
        }
           
        fseek(fileTotalTime, 0, SEEK_END);
        long size = ftell(fileTotalTime);
        if (size == 0)
        {
            fprintf(fileTotalTime, "%-5s %-20s %-12s\n",
                "nP", "nCells", "Total");
        }

        fprintf(fileTotalTime, "%-5u %'-20zd %-12.2f \n",
            nParts,totalCells, maxTimes.total) ;

      

    }
   

}
