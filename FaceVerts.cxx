/*
 * FaceVerts.cxx
 *
 *  Created on: Jul 23, 2024
 *      Author: cfog
 */

#include "FaceVerts.h"

int getTriRotation(const TriFaceVerts &localTri,
		const exa_set<TriFaceVerts> &remote, emInt nDivs) {

	//emInt Lvert0 = localTri.getRemoteIndices(0);
	//emInt Lvert1 = localTri.getRemoteIndices(1);
	//emInt Lvert2 = localTri.getRemoteIndices(2);

	emInt global[3] = { localTri.getGlobalCorner(0), localTri.getGlobalCorner(
			1), localTri.getGlobalCorner(2) };

	TriFaceVerts TF(nDivs, global);
	TF.setCompare(true);
	//TriFaceVerts TF(nDivs,Lvert0,Lvert1,Lvert2);

	auto iterTris = remote.find(TF);
	assert(iterTris!=remote.end());
	emInt vert0 = localTri.getGlobalCorner(0);
	emInt vert1 = localTri.getGlobalCorner(1);
	emInt vert2 = localTri.getGlobalCorner(2);
	// emInt globalIterTris [3]={
	// 	iterTris->getGlobalCorner(0),
	// 	iterTris->getGlobalCorner(1),
	// 	iterTris->getGlobalCorner(2)

	// };
	int rotCase = 0;
	for (int cc = 0; cc < 3; cc++) {
		if (vert0 == iterTris->getGlobalCorner(cc)) {
			if (vert1 == iterTris->getGlobalCorner((cc + 1) % 3)) {
				assert(vert2 == iterTris->getGlobalCorner((cc+2)%3));
				rotCase = cc + 1;
			} else {
				assert(vert1 == iterTris->getGlobalCorner((cc+2)%3));
				assert(vert2 == iterTris->getGlobalCorner((cc+1)%3));
				rotCase = -(cc + 1);
			}
		}
	}
	assert(rotCase != 0);
	return rotCase;

}

int getQuadRotation(const QuadFaceVerts &localQuad,
		const exa_set<QuadFaceVerts> &remote, emInt nDivs) {

	// emInt Lvert0 = localQuad.getRemoteIndices(0);
	// emInt Lvert1 = localQuad.getRemoteIndices(1);
	// emInt Lvert2 = localQuad.getRemoteIndices(2);
	// emInt Lvert3 = localQuad.getRemoteIndices(3);
	emInt global[4] = { localQuad.getGlobalCorner(0), localQuad.getGlobalCorner(
			1), localQuad.getGlobalCorner(2), localQuad.getGlobalCorner(3) };

	QuadFaceVerts QF(nDivs, global);
	QF.setCompare(true);
	//printQuads(remote,nDivs);

	//QuadFaceVerts QF(nDivs,Lvert0,Lvert1,Lvert2,Lvert3);

	auto iterQuads = remote.find(QF);

	assert(iterQuads!=remote.end());

	emInt vert0 = localQuad.getGlobalCorner(0);
	emInt vert1 = localQuad.getGlobalCorner(1);
	emInt vert2 = localQuad.getGlobalCorner(2);
	emInt vert3 = localQuad.getGlobalCorner(3);
	// emInt globalIterTris [4]={
	// 	iterQuads->getGlobalCorner(0),
	// 	iterQuads->getGlobalCorner(1),
	// 	iterQuads->getGlobalCorner(2),
	// 	iterQuads->getGlobalCorner(3)

	// };
	int rotCase = 0;

	for (int cc = 0; cc < 4; cc++) {
		if (vert0 == iterQuads->getGlobalCorner(cc)) {
			if (vert1 == iterQuads->getGlobalCorner((cc + 1) % 4)) {
				// Oriented forward; bdry quad
				assert(vert2 == iterQuads->getGlobalCorner((cc+2)%4));
				assert(vert3 == iterQuads->getGlobalCorner((cc+3)%4));
				rotCase = cc + 1;
			} else {
				assert(vert1 == iterQuads->getGlobalCorner((cc+3)%4));
				assert(vert2 == iterQuads->getGlobalCorner((cc+2)%4));
				assert(vert3 == iterQuads->getGlobalCorner((cc+1)%4));
				rotCase = -(cc + 1);
			}
		}
	}
	assert(rotCase != 0);

	return rotCase;

}

void matchTri(const TriFaceVerts &localTri, const emInt rotation,
		const emInt nDivs, const exa_set<TriFaceVerts> &remoteTriSet,
		std::unordered_map<emInt, emInt> &localRemote) {

	// emInt vert0 = localTri.getRemoteIndices(0);
	// emInt vert1 = localTri.getRemoteIndices(1);
	// emInt vert2 = localTri.getRemoteIndices(2);

	emInt global[3] = { localTri.getGlobalCorner(0), localTri.getGlobalCorner(
			1), localTri.getGlobalCorner(2) };

	TriFaceVerts TF(nDivs, global);
	TF.setCompare(true);

	//TriFaceVerts TF(nDivs,vert0,vert1,vert2);

	auto itRemote = remoteTriSet.find(TF);
	assert(itRemote!=remoteTriSet.end());
	assert(localTri.getPartid()==itRemote->getRemoteId());
	assert(localTri.getRemoteId()==itRemote->getPartid());
	// for(auto i=0; i<3; i++){
	// 	assert(localTri.getCorner(i)==
	// 	itRemote->getRemoteIndices(i));
	// 	assert(itRemote->getCorner(i)==localTri.getRemoteIndices(i));
	// }

	for (emInt ii = 0; ii <= nDivs; ii++) {
		for (emInt jj = 0; jj <= nDivs - ii; jj++) {

			emInt trueI;
			emInt trueJ;
			itRemote->getTrueIJ(ii, jj, trueI, trueJ, rotation);

			emInt vertLocal = localTri.getIntVertInd(trueI, trueJ);
			emInt vertRemote = itRemote->getIntVertInd(ii, jj);
			localRemote.insert( { vertLocal, vertRemote });

			// std::cout<<"ii: "<<ii<<" jj: "<<jj<<"--- ( "<<
			// vertLocal<<" "<<
			// vertRemote<<" )"
			// <<"-> ("<<local->getX(vertLocal)<<", "<<
			// local->getY(vertLocal)<<", "<<local->getZ(vertLocal)<<" ) "<<" ( "<<
			// remote->getX(vertRemote)<<", "<<remote->getY(vertRemote)<<", "<<
			// remote->getZ(vertRemote)<<" ) "
			//   <<std::endl;

		}
	}
	//std::cout<<std::endl;

}

void matchQuad(const QuadFaceVerts &localQuad, const emInt rotation,
		const emInt nDivs, const exa_set<QuadFaceVerts> &remoteQuadSet,
		std::unordered_map<emInt, emInt> &localRemote) {
	//exa_set<TriFaceVerts> remoteTriSet= remote->getRefinedPartTris();

	// emInt vert0 = localQuad.getRemoteIndices(0);
	// emInt vert1 = localQuad.getRemoteIndices(1);
	// emInt vert2 = localQuad.getRemoteIndices(2);
	// emInt vert3 = localQuad.getRemoteIndices(3);

	emInt global[4] = { localQuad.getGlobalCorner(0), localQuad.getGlobalCorner(
			1), localQuad.getGlobalCorner(2), localQuad.getGlobalCorner(3) };
	QuadFaceVerts QF(nDivs, global);
	QF.setCompare(true);

	//QuadFaceVerts QF(nDivs,vert0,vert1,vert2,vert3);

	auto itRemote = remoteQuadSet.find(QF);
	assert(itRemote!=remoteQuadSet.end());
	assert(localQuad.getPartid()==itRemote->getRemoteId());
	assert(localQuad.getRemoteId()==itRemote->getPartid());
	// for(auto i=0; i<4; i++){
	// 	assert(localQuad.getCorner(i)==
	// 	itRemote->getRemoteIndices(i));
	// 	assert(itRemote->getCorner(i)==localQuad.getRemoteIndices(i));
	// }
	for (emInt ii = 0; ii <= nDivs; ii++) {
		for (emInt jj = 0; jj <= nDivs; jj++) {

			emInt trueI;
			emInt trueJ;
			itRemote->getTrueIJ(ii, jj, trueI, trueJ, rotation);
			emInt vertLocal = localQuad.getIntVertInd(ii, jj);
			emInt vertRemote = itRemote->getIntVertInd(trueI, trueJ);
			//emInt vertLocal  = localQuad.getIntVertInd(trueI,trueJ);
			//emInt vertRemote = itRemote->getIntVertInd(ii,jj);
			//std::cout<<vertLocal<<" "<<vertRemote<<std::endl;
			localRemote.insert( { vertLocal, vertRemote });

			// std::cout<<"ii: "<<ii<<" jj: "<<jj<<"--- ( "<<
			// vertLocal<<" "<<
			// vertRemote<<" )"
			// <<"-> ("<<local->getX(vertLocal)<<", "<<
			// local->getY(vertLocal)<<", "<<local->getZ(vertLocal)<<" ) "<<" ( "<<
			// remote->getX(vertRemote)<<", "<<remote->getY(vertRemote)<<", "<<
			// remote->getZ(vertRemote)<<" ) "
			//   <<std::endl;

		}
	}
}

void findRotationAndMatchTris(const TriFaceVerts &localTri, const emInt nDivs,
		const exa_set<TriFaceVerts> &remoteTriSet,
		std::unordered_map<emInt, emInt> &localRemote) {

	emInt global[3] = { localTri.getGlobalCorner(0), localTri.getGlobalCorner(
			1), localTri.getGlobalCorner(2) };

	TriFaceVerts TF(nDivs, global);
	TF.setCompare(true);

	auto iterTris = remoteTriSet.find(TF);
	assert(iterTris!=remoteTriSet.end());
	assert(localTri.getPartid() == iterTris->getRemoteId());
	assert(localTri.getRemoteId()== iterTris->getPartid());

	emInt vert0 = localTri.getGlobalCorner(0);
	emInt vert1 = localTri.getGlobalCorner(1);
	emInt vert2 = localTri.getGlobalCorner(2);

	int rotCase = 0;
	for (int cc = 0; cc < 3; cc++) {
		if (vert0 == iterTris->getGlobalCorner(cc)) {
			if (vert1 == iterTris->getGlobalCorner((cc + 1) % 3)) {
				assert(vert2 == iterTris->getGlobalCorner((cc+2)%3));
				rotCase = cc + 1;
			} else {
				assert(vert1 == iterTris->getGlobalCorner((cc+2)%3));
				assert(vert2 == iterTris->getGlobalCorner((cc+1)%3));
				rotCase = -(cc + 1);
			}
		}
	}
	assert(rotCase != 0);

	// Once found totation, go a haed to match tris

	for (emInt ii = 0; ii <= nDivs; ii++) {
		for (emInt jj = 0; jj <= nDivs - ii; jj++) {

			emInt trueI;
			emInt trueJ;
			iterTris->getTrueIJ(ii, jj, trueI, trueJ, rotCase);

			emInt vertLocal = localTri.getIntVertInd(trueI, trueJ);
			emInt vertRemote = iterTris->getIntVertInd(ii, jj);
			localRemote.insert( { vertLocal, vertRemote });

		}
	}

}

void findRotationAndMatchQuads(const QuadFaceVerts &localQuad,
		const emInt nDivs, const exa_set<QuadFaceVerts> &remoteQuadSet,
		std::unordered_map<emInt, emInt> &localRemote) {

	emInt global[4] = { localQuad.getGlobalCorner(0), localQuad.getGlobalCorner(
			1), localQuad.getGlobalCorner(2), localQuad.getGlobalCorner(3) };

	QuadFaceVerts QF(nDivs, global);
	QF.setCompare(true);

	auto iterQuads = remoteQuadSet.find(QF);

	assert(iterQuads!=remoteQuadSet.end());
	assert(localQuad.getPartid() == iterQuads->getRemoteId());
	assert(localQuad.getRemoteId()== iterQuads->getPartid());

	emInt vert0 = localQuad.getGlobalCorner(0);
	emInt vert1 = localQuad.getGlobalCorner(1);
	emInt vert2 = localQuad.getGlobalCorner(2);
	emInt vert3 = localQuad.getGlobalCorner(3);

	int rotCase = 0;

	for (int cc = 0; cc < 4; cc++) {
		if (vert0 == iterQuads->getGlobalCorner(cc)) {
			if (vert1 == iterQuads->getGlobalCorner((cc + 1) % 4)) {
				// Oriented forward; bdry quad
				assert(vert2 == iterQuads->getGlobalCorner((cc+2)%4));
				assert(vert3 == iterQuads->getGlobalCorner((cc+3)%4));
				rotCase = cc + 1;
			} else {
				assert(vert1 == iterQuads->getGlobalCorner((cc+3)%4));
				assert(vert2 == iterQuads->getGlobalCorner((cc+2)%4));
				assert(vert3 == iterQuads->getGlobalCorner((cc+1)%4));
				rotCase = -(cc + 1);
			}
		}
	}
	assert(rotCase != 0);

	for (emInt ii = 0; ii <= nDivs; ii++) {
		for (emInt jj = 0; jj <= nDivs; jj++) {

			emInt trueI;
			emInt trueJ;
			iterQuads->getTrueIJ(ii, jj, trueI, trueJ, rotCase);
			emInt vertLocal = localQuad.getIntVertInd(ii, jj);
			emInt vertRemote = iterQuads->getIntVertInd(trueI, trueJ);

			localRemote.insert( { vertLocal, vertRemote });

		}
	}
}

void buildTrisMap(hashTri const &tris, std::map<int, vecTri> &remoteTotris,
		std::set<int> &neighbors) {

	for (auto it = tris.begin(); it != tris.end(); it++) {
		int remoteId = it->getRemoteId();
		neighbors.insert(remoteId);
		remoteTotris[remoteId].push_back(*it);

	}
}
void buildQuadsMap(hashQuad const &quads, std::map<int, vecQuad> &remoteToquads,
		std::set<int> &neighbors) {

	for (auto it = quads.begin(); it != quads.end(); it++) {
		int remoteId = it->getRemoteId();
		neighbors.insert(remoteId);
		remoteToquads[remoteId].push_back(*it);

	}
}

void printMultiMap(
		const std::multimap<std::set<emInt>, std::pair<emInt, emInt>> &map) {
	for (const auto &entry : map) {
		const auto &key = entry.first;
		const auto &value = entry.second;

		// Print the key (set of emInt)
		std::cout << "Key: { ";
		for (const auto &elem : key) {
			std::cout << elem << " ";
		}
		std::cout << "} ";

		// Print the value (pair of emInt and char*)
		std::cout << "Value: {" << value.first << ", " << value.second << "}\n";
	}
}

void printTriFaceVerts(const std::vector<TriFaceVerts> &tris) {
	for (const auto &tri : tris) {

		std::cout << "Global corners: [" << tri.getGlobalCorner(0) << ", "
				<< tri.getGlobalCorner(1) << ", " << tri.getGlobalCorner(2)
				<< "]" << std::endl;
		std::cout << "Part ID: " << tri.getPartid() << std::endl;
		std::cout << "Remote Part ID: " << tri.getRemoteId() << std::endl;
		std::cout << "Is True: " << std::boolalpha << tri.getGlobalCompare()
				<< std::endl;
		std::cout << std::endl;
	}
}

void preMatchingPartBdryTris(const emInt numDivs, const setTri &partBdryTris,
		vecVecTri &tris) {
	std::size_t k = 0;
	for (auto itr = partBdryTris.begin(); itr != partBdryTris.end(); itr++) {
		k++;
		auto next = std::next(itr, 1);
		if (k != partBdryTris.size()) {
			if (itr->getSortedGlobal(0) == next->getSortedGlobal(0) &&

			itr->getSortedGlobal(1) == next->getSortedGlobal(1)
					&& itr->getSortedGlobal(2) == next->getSortedGlobal(2)) {

				emInt global[3] = { itr->getGlobalCorner(0),
						itr->getGlobalCorner(1), itr->getGlobalCorner(2) };
				emInt globalNext[3] = { next->getGlobalCorner(0),
						next->getGlobalCorner(1), next->getGlobalCorner(2) };

				TriFaceVerts tripart(numDivs, global, itr->getPartid(),
						next->getPartid(), true);

				TriFaceVerts tripartNext(numDivs, globalNext, next->getPartid(),
						itr->getPartid(), true);

				tris[itr->getPartid()].push_back(tripart);
				tris[next->getPartid()].push_back(tripartNext);
			}
		}
	}
}

void preMatchingPartBdryQuads(const emInt numDivs, const setQuad &partBdryQuads,
		vecVecQuad &quads) {
	std::size_t kquad = 0;
	for (auto itr = partBdryQuads.begin(); itr != partBdryQuads.end(); itr++) {

		auto next = std::next(itr, 1);
		kquad++;
		if (kquad != partBdryQuads.size()) {
			emInt v0Global = next->getGlobalCorner(0);
			emInt v1Global = next->getGlobalCorner(1);
			emInt v2Global = next->getGlobalCorner(2);
			emInt v3Global = next->getGlobalCorner(3);

			emInt partid = next->getPartid();

			emInt v0SortedGlobal = next->getSortedGlobal(0);
			emInt v1SortedGlobal = next->getSortedGlobal(1);
			emInt v2SortedGlobal = next->getSortedGlobal(2);
			emInt v3SortedGlobal = next->getSortedGlobal(3);

			emInt v0Global_ = itr->getGlobalCorner(0);
			emInt v1Global_ = itr->getGlobalCorner(1);
			emInt v2Global_ = itr->getGlobalCorner(2);
			emInt v3Global_ = itr->getGlobalCorner(3);

			emInt partid_ = itr->getPartid();

			emInt v0SortedGlobal_ = itr->getSortedGlobal(0);
			emInt v1SortedGlobal_ = itr->getSortedGlobal(1);
			emInt v2SortedGlobal_ = itr->getSortedGlobal(2);
			emInt v3SortedGlobal_ = itr->getSortedGlobal(3);

			if (v0SortedGlobal_ == v0SortedGlobal
					&& v1SortedGlobal == v1SortedGlobal_
					&& v2SortedGlobal == v2SortedGlobal_
					&& v3SortedGlobal == v3SortedGlobal_) {
				emInt global[4] = { v0Global, v1Global, v2Global, v3Global };
				emInt global_[4] =
						{ v0Global_, v1Global_, v2Global_, v3Global_ };

				QuadFaceVerts quadpart(numDivs, global, partid, partid_, true);
				QuadFaceVerts quadpart_(numDivs, global_, partid_, partid,
						true);
				quads[partid].push_back(quadpart);
				quads[partid_].push_back(quadpart_);
			}
		}
	}
}

void TestPartFaceMatching(const std::size_t nPart, const vecHashTri &HashTri,
		const vecHashQuad &HashQuad, const vecVecTri &vecTris,
		const vecVecQuad &vecquads) {
	assert(HashTri.size()==nPart);
	assert(vecTris.size()==nPart);
	assert(HashQuad.size()==nPart);
	assert(vecquads.size()==nPart);
	for (std::size_t i = 0; i < nPart; i++) {
		assert(vecTris[i].size()==HashTri[i].size());
		assert(vecquads[i].size()==HashQuad[i].size());
		for (auto itri : HashTri[i]) {
			auto find = std::find(vecTris[i].begin(), vecTris[i].end(), itri);
			assert(find!=vecTris[i].end());

		}
		for (auto iquad : HashQuad[i]) {
			auto find = std::find(vecquads[i].begin(), vecquads[i].end(),
					iquad);
			if (find == vecquads[i].end()) {
				assert(find!=vecquads[i].end());
			}
		}

	}
	std::cout << "Successful Testing Part Face Matching" << std::endl;
}

