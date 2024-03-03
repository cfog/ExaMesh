
#include <unistd.h>
#include <cstdio>
#include <boost/mpi.hpp>
#include "ExaMesh.h"
#include "CubicMesh.h"
#include "PARMETIS.h"
#include "UMesh.h"
#include <boost/serialization/unique_ptr.hpp>
#include <cstdio>

void Reduce (boost::mpi::communicator world, 
            const timeResults times, 
            timeResults& maxTimes, 
            const size_t nCells, 
            size_t& totalCells)
{

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

}

void writeEachRankData 
( boost::mpi::communicator  world, FILE* file, timeResults times, size_t nCells)
{
    setlocale(LC_ALL, "");
    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    if(world.rank()!=MASTER)
    {
        times.partfacematching=0;
        times.partition=0;
    }
    if(world.rank()==MASTER)
    {
        if (size == 0) 
        {
        fprintf(file, "%-5s %-10s %-14s %-15s %-12s %-12s %-12s %-12s %-12s %-14s %-14s %-14s %-14s %-16s %-16s\n",
                "nP", "Read", "Partition", "PartFaceMatch", "Sync1", "Sync2", "Serial",
                "Extract", "Refine", "FaceExchange",
                "TriSync", "QuadSync", "MatchTri", "MatchQuad", "Total");
        }
    }


   fprintf(file, "%-5u %-10.2f %-14.2f %-15.2f %-12.2f %-12.2f %-12.2f %-12.2f %-12.2f %-14.2f %-14.2f %-14.2f %-14.2f %-16.2f %-16.2f\n",
          world.rank(), times.read, times.partition, times.partfacematching, times.sync1, times.sync2, times.serial,
           times.extract, times.refine, times.faceExchange,
           times.syncTri, times.syncQuad, times.matchtris, times.matchquads, times.total);
}
void Wait( vecReqs& reqs)
{
    boost::mpi::wait_all(reqs.begin(),reqs.end());
}


void refineForMPI ( const char  baseFileName[] , const char type[], 
                    const char  ugridInfix[]   , const char CGNSFileName[],
                    const int   numDivs        , const char MeshType, std::string mshName, FILE* eachRank
                    )
{

    size_t lastSlashPos        = std::string(baseFileName).find_last_of('/');
	//mshName                    = std::string(baseFileName).substr(lastSlashPos + 1)+"_WeakScaleTimes.txt";
    std::string totalTimeFile  = std::string(baseFileName).substr(lastSlashPos + 1)+"_WeakOnlyTotalTime.txt";

    

    timeResults times, maxTimes;
    times.total=exaTime();
    times.serial=exaTime();
   
    boost::mpi::environment   env; 
	boost::mpi::communicator  world;
    emInt nParts = world.size(); 


    std::string fileName = mshName+"-nDivs-"+std::to_string(numDivs)+
    "-nCPUS-"+std::to_string(world.size())+ "AllTimes.txt";
    FILE *outFileAllTimes = fopen(fileName.c_str(), "a");
    if (outFileAllTimes  == NULL) 
	{
        fprintf(stderr, "Error opening the file!\n");
        exit(1);
    }

    vecVecReqs reqForPartCells(nParts);
    vecVecReqs reqForCoarseFaces(nParts);
    vecVecReqs triReqs  (nParts); 
    vecVecReqs quadReqs (nParts);

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
            Request rq0 =world.isend(irank,1,part2cell); 
            reqForPartCells[world.rank()].emplace_back(rq0);

            Request rq1= world.isend(irank,2,tris[irank]); 
            reqForCoarseFaces[world.rank()].emplace_back(rq1);

            Request rq2= world.isend(irank,3,quads[irank]);
            reqForCoarseFaces[world.rank()].emplace_back(rq2);

        }
        vectorToSet(tris[MASTER],hashTris);
        vectorToSet(quads[MASTER],hashQuads);

    }
    if(world.rank()!=MASTER)
    {
        Request rq0=world.irecv(MASTER,1,partCells);
        reqForPartCells[world.rank()].emplace_back(rq0);

        Request rq1= world.irecv(MASTER,2,triV[world.rank()]);
        reqForCoarseFaces[world.rank()].emplace_back(rq1);

        Request rq2= world.irecv(MASTER,3,quadV[world.rank()]);
        reqForCoarseFaces[world.rank()].emplace_back(rq2);
    }


    times.sync1=exaTime();
    Wait(reqForPartCells[world.rank()]);
    times.sync1=exaTime()-times.sync1;

    times.sync2=exaTime();
    Wait(reqForCoarseFaces[world.rank()]);
    times.sync2=exaTime()-times.sync2;

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


    times.refine=exaTime();
    auto refinedMsh = std::make_unique<UMesh>(*((extractedMsh.get())),numDivs,world.rank());
    times.refine=exaTime()-times.refine;

    times.faceExchange=exaTime();

    auto refinedTris = refinedMsh->getRefinedPartTris();

    buildTrisMap (refinedTris ,remoteTovecTris ,triNeighbrs);

    for(auto tri:remoteTovecTris)
	{
		Request req= world.isend(tri.first,6,tri.second);
		triReqs[world.rank()].push_back(req);
	}
    
    auto refinedQuads= refinedMsh->getRefinedPartQuads();

    buildQuadsMap(refinedQuads,remoteTovecQuads,quadNeighbrs); 
    for(auto quad:remoteTovecQuads)
	{
		Request req = world.isend(quad.first,7,quad.second);
		quadReqs[world.rank()].push_back(req);  
	}
    
	quadsTobeRcvd.resize(quadNeighbrs.size());

    int Quadjj=0; 

	for(auto isource:quadNeighbrs)
	{
		auto findSource = remoteTovecQuads.find(isource); 
		quadsTobeRcvd[Quadjj].resize(findSource->second.size()); 
		Quadjj++; 
	}
    int Quadkk=0 ; 
	for(auto source: quadNeighbrs)
	{
		Request req = world.irecv(source,7,quadsTobeRcvd[Quadkk]);
		quadReqs[world.rank()].push_back(req); 
		Quadkk++; 
	}
    
    trisTobeRcvd.resize (triNeighbrs.size()); 

    int Trijj=0 ; 
	for(auto isource:triNeighbrs)
	{
		auto findSource = remoteTovecTris.find(isource); 
		trisTobeRcvd[Trijj].resize(findSource->second.size()); 
		Trijj++; 
	}

    int Trikk=0 ; 
	for(auto source: triNeighbrs)
	{
	    Request req= world.irecv(source,6,trisTobeRcvd[Trikk]);
		triReqs[world.rank()].push_back(req);
		Trikk++; 
	}

    times.faceExchange=exaTime()-times.faceExchange;

    times.syncTri=exaTime();
    Wait(triReqs[world.rank()]); 
    times.syncTri=exaTime()-times.syncTri;


    times.matchtris=exaTime();

    for(const auto& tri: trisTobeRcvd)
	{
        recvdTris.insert(tri.begin(),tri.end()); 
	}

    for(auto it=refinedTris.begin(); it!=refinedTris.end(); it++)
	{
		std::unordered_map<emInt, emInt> localRemote;
		findRotationAndMatchTris(*it,numDivs,recvdTris,localRemote); 
	}
    times.matchtris=exaTime()-times.matchtris;

    times.syncQuad=exaTime();
    
    Wait(quadReqs[world.rank()]);

    times.syncQuad=exaTime()-times.syncQuad;

    times.matchquads=exaTime();

    for(const auto& quad:quadsTobeRcvd)
	{
		recvdQuads.insert(quad.begin(),quad.end()); 
	}
    for(auto it=refinedQuads.begin(); it!=refinedQuads.end();it++)
	{
		std::unordered_map<emInt,emInt> localRemote; 
		findRotationAndMatchQuads(*it,numDivs,recvdQuads,localRemote); 
    }
    times.matchquads=exaTime()-times.matchquads;

    times.total=exaTime()-times.total;

    size_t         nCells = refinedMsh->numPyramids()+
    refinedMsh->numPrisms()+refinedMsh->numHexes()+refinedMsh->numTets()+ 
    refinedMsh->numBdryTris()+refinedMsh->numBdryQuads();

    size_t         totalCells;

    world.barrier();    

    boost::mpi::reduce(world, size_t(nCells), totalCells ,     std::plus<size_t>(),            MASTER);
    boost::mpi::reduce(world, times.total,    maxTimes.total,  boost:: mpi::maximum<double>(), MASTER);
    boost::mpi::reduce(world, times.serial,   maxTimes.serial, boost:: mpi::maximum<double>(), MASTER);


    if(world.rank()==MASTER)
    {
        FILE *file = fopen(totalTimeFile.c_str(), "a");
        if (file == NULL) 
            fprintf(stderr, "Error opening the file!\n");
        fseek(file, 0, SEEK_END);
        long size = ftell(file);
        if (size == 0)  
            fprintf(file, "%-5s %-20s %-12s %-12s\n", "nP", "nRefinedCells", "MaxSerial" ,"MaxTotal");
        fprintf(file, "%-5u %'-20zd %-12.2f %-12.2f \n",
            nParts,totalCells, maxTimes.serial,maxTimes.total);
    }
  
    writeEachRankData(world,outFileAllTimes,times,partCells[world.rank()].size());


}
