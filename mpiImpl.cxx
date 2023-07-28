 #include "mpiImpl.h"

void
mpiImpl:: 
refineMPI(const int numDivs)
{
    // boost::mpi::environment   env; 
	// boost::mpi::communicator  world;

	// std::vector<boost::mpi::request> Trireqs;  

	// vecPart         parts; 
	// vecCellPartData vecCPD; 
	
	// std::size_t    vecCPDSize; 
	// std::size_t    trisSize; 
	// std::size_t    quadsSize; 
	// std::size_t nParts = world.size();

	// hashTri  trisS; 
	// hashQuad quadsS;

	// vecTri   triV;
	// vecQuad  quadV;  

	// int MASTER = 0; 
	// int tag    = 0; 

	// intToVecTri  remoteTovecTris;
	// intToVecQuad remoteTovecQuads; 

	// std::set<int> triNeighbrs;  
	// std::set<int> quadNeighbrs; 

	
	// std::vector<vecTri>  trisTobeRcvd; 
	// std::vector<vecQuad> quadsTobeRcvd;

	// hashTri  recvdTris;
	// hashQuad recvdQuads;  

	// if(world.rank()==MASTER)
	// {
	// 	partitionCells(mCoarseMesh, nParts, parts,vecCPD); 

	// 	vecCPDSize = vecCPD.size(); 

	// 	assert(vecCPDSize>0); 

	// 	vecHashTri  VectrisHash; 
	// 	vecHashQuad VecquadsHash;

	// 	vecVecTri   VecTriVec; 
	// 	vecVecQuad  vecQuadVec; 

	// 	mCoarseMesh->partFaceMatching(parts,vecCPD,VectrisHash,VecquadsHash); 

	// 	for(auto  itri=0 ; itri<VectrisHash.size(); itri++)
	// 	{
	// 		vecTri TriVec; 
	// 		SetToVector(VectrisHash[itri],TriVec); 
	// 		VecTriVec.emplace_back(TriVec); 
	// 	}
	// 	for(auto iquad=0 ; iquad<VecquadsHash.size(); iquad++)
	// 	{
	// 		vecQuad QuadVec; 
	// 		SetToVector(VecquadsHash[iquad],QuadVec); 
	// 		vecQuadVec.emplace_back(QuadVec); 
	// 	}

	// 	trisS = VectrisHash[0]; // For MASTER 
	// 	quadsS= VecquadsHash[0];
		
	// 	for(auto irank=1 ; irank<world.size();irank++)
	// 	{
	// 		world.send(irank,tag,parts[irank]); 
	// 		world.send(irank,tag,VecTriVec[irank]); 
	// 		world.send(irank,tag,vecQuadVec[irank]); 
	// 	}
	// }
	// else
	// {
	// 	parts.resize(world.size()); 
	// 	world.recv(MASTER,tag,parts[world.rank()]);
	// 	world.recv(MASTER,tag,triV); 
	// 	world.recv(MASTER,tag,quadV);
	// 	vectorToSet(triV,trisS);
	// 	vectorToSet(quadV,quadsS); 
	// }

	// boost::mpi::broadcast(world,vecCPDSize,MASTER);

	// if(world.rank()==MASTER)
	// {
	// 	for(auto irank=1 ; irank<world.size();irank++)
	// 	{
	// 		world.send(irank,tag,vecCPD); 
	// 	}
	// }
	// if(world.rank()!= MASTER)
	// {
	// 	vecCPD.resize(vecCPDSize); 
	// 	world.recv(MASTER,tag,vecCPD); 
	// }

	// auto coarse= mCoarseMesh->APIextractCoarseMesh(parts[world.rank()],vecCPD,numDivs, 
	//  trisS,quadsS,world.rank()); 


	// auto refinedMesh = std::make_shared<UMesh>(
	//  		*(coarse.get()), numDivs, world.rank());

	// auto tris  = refinedMesh->getRefinedPartTris();

	// auto quads = refinedMesh->getRefinedPartQuads();
	
	// buildTrisMap(tris,remoteTovecTris,triNeighbrs);

	// trisTobeRcvd.resize(triNeighbrs.size());

	// int jj=0 ; 
	// for(auto isource:triNeighbrs)
	// {
	// 	auto findSource = remoteTovecTris.find(isource); 
	// 	// same amount of message will be sent and received from other prcoeszsor
	// 	// Be cuatious if this assumption later on will be break
	// 	trisTobeRcvd[jj].resize(findSource->second.size()); 
	// 	jj++; 

	// }

	// for(auto tri:remoteTovecTris)
	// {

	// 	boost::mpi::request req= world.isend(tri.first,tag,tri.second);
	// 	Trireqs.push_back(req);
	
	// }

	// int kk=0 ; 
	// for(auto source: triNeighbrs)
	// {
		
	// 	boost::mpi::request req= world.irecv(source,tag,trisTobeRcvd[kk]);
	// 	//req.test();
	// 	Trireqs.push_back(req);
	// 	kk++; 
	// }

	// boost::mpi::wait_all(Trireqs.begin(),Trireqs.end()); 

	// //boost::mpi::test_all(reqs.begin(),reqs.end());

	// for(const auto& tri: trisTobeRcvd)
	// {

	// 	recvdTris.insert(tri.begin(),tri.end()); 
	// }

	// for(auto it=tris.begin(); it!=tris.end(); it++)
	// {
	// 	std::unordered_map<emInt, emInt> localRemote;
	// 	int rotation = getTriRotation(*it,recvdTris,numDivs);
	// 	matchTri(*it,rotation,numDivs,recvdTris,localRemote); 
		

	// }
	

}