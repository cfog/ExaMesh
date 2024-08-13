#include <unistd.h>
#include <cstdio>
#include <boost/mpi.hpp>
#include "ExaMesh.h"
#include "CubicMesh.h"
#include "PARMETIS.h"
#include "UMesh.h"
#include <boost/serialization/unique_ptr.hpp>
#include <boost/mpi/collectives.hpp>
#include <cstdio>
#include <mpi.h>

void writeForMaster(
		// std::string maxFile,
		// std::string avgFile,        
		// std::string totalFile, 
		const size_t nParts, const size_t totalCells,
		const timeResults maxTimes, timeResults avgTimes,
		const timeResults times, const size_t nDivs,
		const size_t nAvgRefinedCells, const size_t nAvgTriFaces,
		const size_t nAvgQuadFaces, const size_t nAvgInitalcells,
		const std::string baseFileName) {

	size_t lastSlashPos = baseFileName.find_last_of('/');
	std::string sTrailingName = baseFileName.substr(lastSlashPos + 1);
	std::string sMaxTimes = sTrailingName + "_MAXTimes.xlsx";
	std::string sAvgTimes = sTrailingName + "_AvgTimes.xlsx";
	std::string sAvgPartData = sTrailingName + "_AvgPartData.csv";
	std::string sAvgExtract = sTrailingName + "_AvgExtract.xlsx";
	std::string sAvgRefine = sTrailingName + "_AvgRefine.xlsx";
	std::string sAvgTotal = sTrailingName + "_AvgTotal.xlsx";

	auto fMaxTimes = openFile(sMaxTimes);
	auto fAvgTimes = openFile(sAvgTimes);
	auto fAvgPart = openFile(sAvgPartData);
	auto fAvgExtract = openFile(sAvgExtract);
	auto fAvgRefine = openFile(sAvgRefine);
	auto fAvgTotal = openFile(sAvgTotal);

	// avgTimes.serial= (times.preProcessing)*nParts+(times.partition)*nParts+
	//(times.InitialFaceMatching)*nParts+avgTimes.broadcasting; 

	// avgTimes.total= avgTimes.serial+avgTimes.extract+avgTimes.refine+
	//                avgTimes.faceExchange+avgTimes.totalFacesWait+
	//                avgTimes.totalMatch;  

	double maxRate = ((double) totalCells / (double) 1000000000
			/ (maxTimes.total / 60)) / (double) nParts;
	double avgRate = ((double) totalCells / (double) 1000000000
			/ (avgTimes.total / 60)) / (double) nParts;

	if (getSeekFile(fMaxTimes) == 0) {
		fprintf(fMaxTimes,
				"%-5s %-10s %-14s %-15s %-12s %-12s %-12s %-12s %-14s %-14s %-16s %-12s %-20s %-12s\n",
				"Rank", "PreProcessing", "Partitioning", "IniFaceMatching",
				"Broadcasting", "Serial", "Extracting", "Refining",
				"FaceExchanging", "Wait", "FaceMatching", "Total", "TotalCells",
				"Rate");
	}

	fprintf(fMaxTimes,
			"%-5lu %-10.2f %-14.2f %-15.2f %-12.2f %-12.2f %-12.2f %-12.2f %-14.2f %-14.2f %-16.2f %-12.2f %'-20zd %-12.2f\n",
			nParts, times.preProcessing, times.partition,
			times.InitialFaceMatching, maxTimes.broadcasting, maxTimes.serial,
			maxTimes.extract, maxTimes.refine, maxTimes.faceExchange,
			maxTimes.totalFacesWait, maxTimes.totalMatch, maxTimes.total,
			totalCells, maxRate);

	if (getSeekFile(fAvgTimes) == 0) {
		fprintf(fAvgTimes,
				"%-5s %-10s %-14s %-15s %-12s %-12s %-12s %-12s %-14s %-14s %-16s %-12s %-20s %-12s\n",
				"nP", "Read", "Partitioning", "IniFaceMatching", "Broadcasting",
				"Serial", "Extracting", "Refining", "FaceExchanging", "Wait",
				"FaceMatching", "Total", "TotalCells", "Rate");
	}

	double dParts = nParts;
	fprintf(fAvgTimes,
			"%-5lu %-10.2f %-14.2f %-15.2f %-12.2f %-12.2f %-12.2f %-12.2f %-14.2f %-14.2f %-16.2f %-12.2f %'-20zd %-12.2f\n",
			nParts, times.preProcessing, times.partition,
			times.InitialFaceMatching, avgTimes.broadcasting, avgTimes.serial,
			avgTimes.extract, avgTimes.refine, avgTimes.faceExchange,
			avgTimes.totalFacesWait, avgTimes.totalMatch, avgTimes.total,
			totalCells, avgRate);

	if (getSeekFile(fAvgPart) == 0)
		fprintf(fAvgPart, "%-5s %-5s %-14s %-14s %-14s\n", "nP", "nDivs",
				"nInitialCells", "nRefinedCells", "nCoarseFaces");

	fprintf(fAvgPart, "%-5lu %-5lu %-14zu %-14zu %-14zu\n", nParts, nDivs,
			nAvgInitalcells, nAvgRefinedCells, nAvgTriFaces + nAvgQuadFaces);

	if (getSeekFile(fAvgExtract) == 0)
		fprintf(fAvgExtract, "%-5s %-5s \n", "nP", "Extracting");

	fprintf(fAvgExtract, "%-5lu %-5.2f \n", nParts, avgTimes.extract);

	if (getSeekFile(fAvgRefine) == 0)
		fprintf(fAvgRefine, "%-5s %-5s \n", "nP", "Refining");

	fprintf(fAvgRefine, "%-5lu %-5.2f \n", nParts, avgTimes.refine);

	if (getSeekFile(fAvgTotal) == 0)
		fprintf(fAvgTotal, "%-5s %-12s\n", "nP", "Total");
	fprintf(fAvgTotal, "%-5lu %-12.2f\n", nParts, avgTimes.total);
}

void writePartData(boost::mpi::communicator world,
		const std::string baseFileName,
		std::vector<std::vector<emInt>> &partCells, size_t numDivs,
		size_t nRefinedCells, size_t nTriFaces, size_t nQuadFaces) {
	size_t lastSlashPos = std::string(baseFileName).find_last_of('/');

	std::string fileName = baseFileName.substr(lastSlashPos + 1)
			+ "-nDivs-" + std::to_string(numDivs) + "-nCPUS-"
			+ std::to_string(world.size()) + "_PartData.txt";

	FILE *file = fopen(fileName.c_str(), "a");

	if (file == NULL)
		fprintf(stderr, "Error opening the file!\n");
	fseek(file, 0, SEEK_END);
	long size = ftell(file);

	if (world.rank() == MASTER) {
		if (size == 0)
			fprintf(file, "%-5s %-14s %-14s %-14s\n", "Rank", "nInitialCells",
					"nRefinedCells", "nCoarseFaces");
	}
	fprintf(file, "%-5u %-14zu %-14zu %-14zu\n", world.rank(),
			partCells[world.rank()].size(), nRefinedCells,
			nTriFaces + nQuadFaces);

}
void writeEachRankData(boost::mpi::communicator world, FILE *file,
		timeResults times, size_t nCells) {
	setlocale(LC_ALL, "");
	fseek(file, 0, SEEK_END);
	long size = ftell(file);

	if (world.rank() == MASTER) {
		if (size == 0) {
			fprintf(file,
					"%-5s %-10s %-14s %-15s %-12s %-12s %-12s %-12s %-14s %-14s %-16s %-16s\n",
					"Rank", "PreProcessing", "Partitioning", "IniFaceMatching",
					"Broadcasting", "Serial", "Extracting", "Refining",
					"FaceExchanging", "Wait", "FaceMatching", "Total");
		}
	}
	fprintf(file,
			"%-5u %-10.2f %-14.2f %-15.2f %-12.2f %-12.2f %-12.2f %-12.2f %-14.2f %-14.2f %-16.2f %-16.2f\n",
			world.rank(), times.preProcessing, times.partition,
			times.InitialFaceMatching, times.broadcasting, times.serial,
			times.extract, times.refine, times.faceExchange,
			times.totalFacesWait, times.totalMatch, times.total);
}
void Wait(vecReqs &reqs) {
	boost::mpi::wait_all(reqs.begin(), reqs.end());
}

void ReduceToMax(boost::mpi::communicator world, timeResults times,
		timeResults &maxTimes) {

	boost::mpi::reduce(world, times.total, maxTimes.total,
			boost::mpi::maximum<double>(), MASTER);
	boost::mpi::reduce(world, times.serial, maxTimes.serial,
			boost::mpi::maximum<double>(), MASTER);
	boost::mpi::reduce(world, times.preProcessing, maxTimes.preProcessing,
			boost::mpi::maximum<double>(), MASTER);
	boost::mpi::reduce(world, times.extract, maxTimes.extract,
			boost::mpi::maximum<double>(), MASTER);
	boost::mpi::reduce(world, times.refine, maxTimes.refine,
			boost::mpi::maximum<double>(), MASTER);
	boost::mpi::reduce(world, times.faceExchange, maxTimes.faceExchange,
			boost::mpi::maximum<double>(), MASTER);
	boost::mpi::reduce(world, times.broadcasting, maxTimes.broadcasting,
			boost::mpi::maximum<double>(), MASTER);
	boost::mpi::reduce(world, times.totalFacesWait, maxTimes.totalFacesWait,
			boost::mpi::maximum<double>(), MASTER);
	boost::mpi::reduce(world, times.totalMatch, maxTimes.totalMatch,
			boost::mpi::maximum<double>(), MASTER);
	// boost::mpi::reduce(world, times.vecToSetTime,   maxTimes.vecToSetTime,   boost:: mpi::maximum<double>(),   MASTER);

}

void ReduceToSUM(boost::mpi::communicator world, timeResults times,
		timeResults &avgTimes) {

	MPI_Reduce(&times.total, &avgTimes.total, 1, MPI_DOUBLE, MPI_SUM, MASTER,
			world);
	// MPI_Reduce( &times.vecToSetTime,   &avgTimes.vecToSetTime,  1, MPI_DOUBLE, MPI_SUM, MASTER , world);
	MPI_Reduce(&times.serial, &avgTimes.serial, 1, MPI_DOUBLE, MPI_SUM, MASTER,
			world);
	//MPI_Reduce( &times.preProcessing,  &avgTimes.preProcessing, 1, MPI_DOUBLE, MPI_SUM, MASTER , world);
	MPI_Reduce(&times.extract, &avgTimes.extract, 1, MPI_DOUBLE, MPI_SUM,
			MASTER, world);
	MPI_Reduce(&times.refine, &avgTimes.refine, 1, MPI_DOUBLE, MPI_SUM, MASTER,
			world);
	MPI_Reduce(&times.faceExchange, &avgTimes.faceExchange, 1, MPI_DOUBLE,
			MPI_SUM, MASTER, world);
	MPI_Reduce(&times.broadcasting, &avgTimes.broadcasting, 1, MPI_DOUBLE,
			MPI_SUM, MASTER, world);
	MPI_Reduce(&times.totalFacesWait, &avgTimes.totalFacesWait, 1, MPI_DOUBLE,
			MPI_SUM, MASTER, world);
	MPI_Reduce(&times.totalMatch, &avgTimes.totalMatch, 1, MPI_DOUBLE, MPI_SUM,
			MASTER, world);

	avgTimes.total = avgTimes.total / double(world.size());
	// avgTimes.preProcessing  =  avgTimes.preProcessing/double(world.size());
	avgTimes.broadcasting = avgTimes.broadcasting / double(world.size());
	avgTimes.serial = avgTimes.serial / double(world.size());
	avgTimes.extract = avgTimes.extract / double(world.size());
	avgTimes.refine = avgTimes.refine / double(world.size());
	avgTimes.faceExchange = avgTimes.faceExchange / double(world.size());
	avgTimes.totalFacesWait = avgTimes.totalFacesWait / double(world.size());
	avgTimes.totalMatch = avgTimes.totalMatch / double(world.size());
	avgTimes.total = avgTimes.total / double(world.size());

}

void ReducePartDataToAvg(boost::mpi::communicator world,
		const size_t nRefinedCells, const size_t nTriFaces,
		const size_t nQuadFaces, size_t &nAvgRefinedCells, size_t &nAvgTriFaces,
		size_t &nAvgQuadFaces, size_t &nAvgInitalcells, emInt nParts,
		std::vector<std::vector<emInt>> &partCells) {

	MPI_Reduce(&nRefinedCells, &nAvgRefinedCells, 1, MPI_UNSIGNED_LONG, MPI_SUM,
			MASTER, world);
	MPI_Reduce(&nTriFaces, &nAvgTriFaces, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER,
			world);
	MPI_Reduce(&nQuadFaces, &nAvgQuadFaces, 1, MPI_UNSIGNED_LONG, MPI_SUM,
			MASTER, world);

	nAvgRefinedCells /= nParts;
	nAvgTriFaces /= nParts;
	nAvgQuadFaces /= nParts;

	size_t sum = 0;
	for (const auto &cells : partCells) {
		sum += cells.size();
	}
	nAvgInitalcells = sum / nParts;

}

void refineForMPI(const std::string  baseFileName,
		const std::string fileSuffix,
		const std::string ugridInfix,
		const emInt   numDivs,
		std::string mshName) {

	size_t lastSlashPos = std::string(baseFileName).find_last_of('/');
	//mshName                    = std::string(baseFileName).substr(lastSlashPos + 1)+"_WeakScaleTimes.txt";
	//std::string maxFile       = std::string(baseFileName).substr(lastSlashPos + 1)+"_MAXTimes.xlsx";
	//std::string avgFile       = std::string(baseFileName).substr(lastSlashPos + 1)+"_AvgTimes.xlsx";
	//std::string totalFile     = std::string(baseFileName).substr(lastSlashPos + 1)+"_TotalTime.xlsx";

	timeResults times, maxTimes, avgTimes;
	//MeshStatics mshData, allMshData;

	boost::mpi::environment env;
	boost::mpi::communicator world;

	times.serial = exaTime();
	times.total = exaTime();
	emInt nParts = world.size();


	emInt sizeCellparts;
	size_t nCoarseTris;
	size_t nCoarseQuads;

	size_t nAvgTriFaces, nAvgQuadFaces, nAvgInitalcells, nAvgRefinedCells;
	std::string fileName = mshName + "-nDivs-" + std::to_string(numDivs)
			+ "-nCPUS-" + std::to_string(world.size())
			+ "TimesForAllRanks.xlsx";
	FILE *fTimeForEachRank = openFile(fileName);

	vecVecReqs triReqs(nParts);
	vecVecReqs quadReqs(nParts);

	std::vector<std::vector<emInt>> partCells;

	// All processors read the mesh 
	times.preProcessing = exaTime();
	std::unique_ptr<ExaMesh> inimesh = ExaMesh::readMeshFromFile(
			baseFileName, fileSuffix, ugridInfix);
	times.preProcessing = exaTime() - times.preProcessing;

	hashTri hashTris;
	hashQuad hashQuads;

	{
		vecVecTri tris;
		vecVecQuad quads;

		if (world.rank() == MASTER) {
			times.partition = exaTime();
			std::vector<emInt> vaicelltopart;
			// Paritions the mesh
			auto part2cell = partitionMetis(inimesh.get(), nParts,
					vaicelltopart);
			partCells = part2cell;
			times.partition = exaTime() - times.partition;

			// Find coarse part boundary faces and update info of part ID and remote part ID
			times.InitialFaceMatching = exaTime();
			sizeCellparts = inimesh->FastpartFaceMatching(nParts, partCells,
					vaicelltopart, tris, quads);
			times.InitialFaceMatching = exaTime() - times.InitialFaceMatching;

		}

		times.broadcasting = exaTime();
		// Broadcasting all data
		boost::mpi::broadcast(world, partCells, MASTER);
		boost::mpi::broadcast(world, tris, MASTER);
		boost::mpi::broadcast(world, quads, MASTER);
		vectorToSet(tris[world.rank()], hashTris);
		vectorToSet(quads[world.rank()], hashQuads);
		times.broadcasting = exaTime() - times.broadcasting;
		times.serial = exaTime() - times.serial;

		// tris and quads can go out of scope now
	}

	//Extract
	times.extract = exaTime();
	std::unique_ptr<ExaMesh> extractedMsh = inimesh->extractCoarseMeshMPI(
			world.rank(), partCells[world.rank()], numDivs, hashTris,
			hashQuads);
	times.extract = exaTime() - times.extract;

	//Refine 
	times.refine = exaTime();
	auto refinedMsh = extractedMsh->subdivideMesh(numDivs, world.rank());

	times.refine = exaTime() - times.refine;

	// Beginning of face exchange
	times.faceExchange = exaTime();

	// Send and post receives for tris
	std::vector<vecTri> trisTobeRcvd;
	auto refinedTris = refinedMsh->getRefinedPartTris();
	{
		std::set<int> triNeighbrs;
		intToVecTri remoteTovecTris;

		// Assemble and send tri data
		{
			// Establishing the list of neighbors whom processor share tris
			buildTrisMap(refinedTris, remoteTovecTris, triNeighbrs);

			for (auto tri : remoteTovecTris) {
				Request req = world.isend(tri.first, 6, tri.second);
				triReqs[world.rank()].push_back(req);
			}
		}
		trisTobeRcvd.resize(triNeighbrs.size());

		// Post receives for tris
		{
			int Trijj = 0;
			for (auto isource : triNeighbrs) {
				auto findSource = remoteTovecTris.find(isource);
				trisTobeRcvd[Trijj].resize(findSource->second.size());
				Trijj++;
			}

			int Trikk = 0;
			for (auto source : triNeighbrs) {
				Request req = world.irecv(source, 6, trisTobeRcvd[Trikk]);
				triReqs[world.rank()].push_back(req);
				Trikk++;
			}
		}
	}

	// Assemble and send quad data
	std::vector<vecQuad> quadsTobeRcvd;
	auto refinedQuads = refinedMsh->getRefinedPartQuads();
	{
		std::set<int> quadNeighbrs;
		intToVecQuad remoteTovecQuads;

		{

			// Establishing the list of neighbors whom processor share quads
			buildQuadsMap(refinedQuads, remoteTovecQuads, quadNeighbrs);
			for (auto quad : remoteTovecQuads) {
				Request req = world.isend(quad.first, 7, quad.second);
				quadReqs[world.rank()].push_back(req);
			}
		}

		// Post receives for quads
		quadsTobeRcvd.resize(quadNeighbrs.size());
		{
			int Quadjj = 0;

			for (auto isource : quadNeighbrs) {
				auto findSource = remoteTovecQuads.find(isource);
				quadsTobeRcvd[Quadjj].resize(findSource->second.size());
				Quadjj++;
			}
			int Quadkk = 0;
			for (auto source : quadNeighbrs) {
				Request req = world.irecv(source, 7, quadsTobeRcvd[Quadkk]);
				quadReqs[world.rank()].push_back(req);
				Quadkk++;
			}
		}
	}

	times.faceExchange = exaTime() - times.faceExchange;

	// Wait for tri data to arrive, and store it
	{
		times.waitTri = exaTime();
		Wait(triReqs[world.rank()]);
		times.waitTri = exaTime() - times.waitTri;

		times.matchtris = exaTime();
		hashTri recvdTris;

		for (const auto &tri : trisTobeRcvd) {
			recvdTris.insert(tri.begin(), tri.end());
		}

		for (auto it = refinedTris.begin(); it != refinedTris.end(); it++) {
			std::unordered_map<emInt, emInt> localRemote;
			findRotationAndMatchTris(*it, numDivs, recvdTris, localRemote);
		}
		times.matchtris = exaTime() - times.matchtris;
	}

	// Wait for quad data to arrive, and store it
	{
		times.waitQuad = exaTime();

		Wait(quadReqs[world.rank()]);

		times.waitQuad = exaTime() - times.waitQuad;

		times.matchquads = exaTime();
		hashQuad recvdQuads;

		for (const auto &quad : quadsTobeRcvd) {
			recvdQuads.insert(quad.begin(), quad.end());
		}
		for (auto it = refinedQuads.begin(); it != refinedQuads.end(); it++) {
			std::unordered_map<emInt, emInt> localRemote;
			findRotationAndMatchQuads(*it, numDivs, recvdQuads, localRemote);
		}
		times.matchquads = exaTime() - times.matchquads;
	}
	times.total = exaTime() - times.total;

	// Post Processing 

	times.totalFacesWait = times.waitTri + times.waitQuad;
	times.totalMatch = times.matchtris + times.matchquads;

	size_t nCells = refinedMsh->numPyramids() + refinedMsh->numPrisms()
			+ refinedMsh->numHexes() + refinedMsh->numTets()
			+ refinedMsh->numBdryTris() + refinedMsh->numBdryQuads();

	size_t totalCells;

	world.barrier();

	if (world.rank() != MASTER) {
		times.InitialFaceMatching = 0;
		times.partition = 0;
		times.preProcessing = 0;
	}
	writeEachRankData(world, fTimeForEachRank, times,
			partCells[world.rank()].size());

	boost::mpi::reduce(world, size_t(nCells), totalCells, std::plus<size_t>(),
			MASTER);
	boost::mpi::reduce(world, times.total, maxTimes.total,
			boost::mpi::maximum<double>(), MASTER);
	boost::mpi::reduce(world, times.serial, maxTimes.serial,
			boost::mpi::maximum<double>(), MASTER);

	ReduceToSUM(world, times, avgTimes);
	ReduceToMax(world, times, maxTimes);
	ReducePartDataToAvg(world, nCells, refinedTris.size(), refinedQuads.size(),
			nAvgRefinedCells, nAvgTriFaces, nAvgQuadFaces, nAvgInitalcells,
			nParts, partCells);

	if (world.rank() == MASTER) {
		writeForMaster(nParts, totalCells, maxTimes, avgTimes, times, numDivs,
				nAvgRefinedCells, nAvgTriFaces, nAvgQuadFaces, nAvgInitalcells,
				baseFileName);

		auto file = openFile("CellPartData.txt");
		fseek(file, 0, SEEK_END);
		long size = ftell(file);
		if (size == 0)
			fprintf(file, "%-5s %-10s\n", "nP", "CellDataSize");

		fprintf(file, "%-5u %-10u \n", nParts, sizeCellparts);

	}

	writePartData(world, baseFileName, partCells, numDivs, nCells,
			refinedTris.size(), refinedQuads.size());

	std::string outFile = "outputfile" + std::to_string(world.rank()) + ".vtk";
	printf("Writing %s\n", outFile.c_str());
	refinedMsh->writeVTKFile(outFile.c_str());
}
