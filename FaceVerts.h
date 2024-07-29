/*
 * FaceVerts.h
 *
 *  Created on: Jul 23, 2024
 *      Author: cfog
 */

#ifndef FACEVERTS_H_
#define FACEVERTS_H_

#include "exa-defs.h"

class FaceVerts {
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int /*version*/) {
		ar & m_remoteId;
		ar & m_global;
		ar & m_sortedGlobal;
		//ar &m_remote;
		//ar &m_sortedRemote;
		//ar &m_corners;
		//ar &m_sorted;
		//ar &m_cornerUVW;
		ar & m_nCorners;
		ar & m_nDivs;
		ar & m_intVerts;
		//ar &m_param_st;
		//ar &m_param_uvw;
		//ar &m_volElem;
		//ar &m_volElemType;
		//ar &m_bothSidesDone;
		ar & m_partId;
		ar & m_globalComparison;

	}
protected:
	// CUATION: ANY CHANGE IN MEMBER SIZE OR DELETION OR ADDING SHOULD BE REFLECTED IN ITS MPI TYPE

	emInt m_global[4], m_sortedGlobal[4];
	emInt m_remote[4], m_sortedRemote[4];
	emInt m_corners[4], m_sorted[4];
	double m_cornerUVW[4][3];
	emInt m_nCorners, m_nDivs;
	std::vector<std::vector<emInt>> m_intVerts;
	//emInt m_intVerts [MAX_DIVS + 1][MAX_DIVS + 1];
	//double m_param_st[MAX_DIVS + 1][MAX_DIVS + 1][2];
	std::vector<std::vector<std::vector<double>>> m_param_st;
	//double m_param_uvw[MAX_DIVS + 1][MAX_DIVS + 1][3];
	std::vector<std::vector<std::vector<double>>> m_param_uvw;
	emInt m_volElem, m_volElemType;
	bool m_bothSidesDone;
	emInt m_partId;
	emInt m_remoteId;
	bool m_globalComparison;
public:
	FaceVerts() :
			m_nCorners(0), m_nDivs(0), m_volElem(EMINT_MAX), m_volElemType(
					TETRA_4), m_bothSidesDone(false), m_partId(EMINT_MAX), m_remoteId(
			EMINT_MAX), m_globalComparison(false) {
	}
	;
	FaceVerts(const int nDivs, const emInt NC = 0) :
			m_nCorners(NC), m_nDivs(nDivs), m_intVerts(nDivs + 1,
					std::vector<emInt>(nDivs + 1)), m_param_st(nDivs + 1,
					std::vector<std::vector<double>>(nDivs + 1,
							std::vector<double>(2))), m_param_uvw(nDivs + 1,
					std::vector<std::vector<double>>(nDivs + 1,
							std::vector<double>(3))), m_volElem(EMINT_MAX), m_volElemType(
					0), m_bothSidesDone(false), m_partId(EMINT_MAX), m_remoteId(
			EMINT_MAX), m_globalComparison(false)

	{
		assert(NC == 3 || NC == 4);
	}
	//virtual ~FaceVerts() {};
	bool isValidIJ(const emInt ii, const emInt jj) const {
		bool retVal = true;
		if (m_nCorners == 3) {
			retVal = retVal && ((ii + jj) <= m_nDivs);
		} else {
			assert(m_nCorners == 4);
			retVal = retVal && (ii <= m_nDivs) && (jj <= m_nDivs);
		}
		return retVal;
	}
	bool isValidParam(const double param) const {
		// This isn't comprehensive, in that not all faces
		// can have this full parameter range.  But the more
		// accurate test requires significantly more information.
		return (param >= 0 && param <= 1);
	}
	//virtual void setupSorted() = 0;
	void setupSorted() {
	}
	;
	emInt getSorted(const int ii) const {
		return m_sorted[ii];
	}
	void setIntVertInd(const emInt ii, const emInt jj, const emInt vert) {
		assert(isValidIJ(ii, jj));
		m_intVerts[ii][jj] = vert;
	}
	emInt getIntVertInd(const emInt ii, const emInt jj) const {
		assert(isValidIJ(ii, jj));
		return m_intVerts[ii][jj];
	}
	void setVertSTParams(const emInt ii, const emInt jj, const double st[2]) {
		assert(isValidIJ(ii, jj));
		assert(isValidParam(st[0]));
		assert(isValidParam(st[1]));
		m_param_st[ii][jj][0] = st[0];
		m_param_st[ii][jj][1] = st[1];
	}
	//virtual void getVertAndST(const int ii, const int jj, emInt &vert,
	//					  double st[2], const int rotCase = 0) const = 0;
	void setVertUVWParams(const emInt ii, const emInt jj, const double uvw[3]) {
		assert(isValidIJ(ii, jj));
		assert(isValidParam(uvw[0]));
		assert(isValidParam(uvw[1]));
		assert(isValidParam(uvw[2]));
		m_param_uvw[ii][jj][0] = uvw[0];
		m_param_uvw[ii][jj][1] = uvw[1];
		m_param_uvw[ii][jj][2] = uvw[2];
	}
	void getVertUVWParams(const emInt ii, const emInt jj, double uvw[3]) const {
		assert(isValidIJ(ii, jj));
		uvw[0] = m_param_uvw[ii][jj][0];
		uvw[1] = m_param_uvw[ii][jj][1];
		uvw[2] = m_param_uvw[ii][jj][2];
		assert(isValidParam(uvw[0]));
		assert(isValidParam(uvw[1]));
		assert(isValidParam(uvw[2]));
	}
	// void setCorners(const emInt cA, const emInt cB, const emInt cC,
	// 				const emInt cD = EMINT_MAX)
	// {
	// 	m_corners[0] = cA;
	// 	m_corners[1] = cB;
	// 	m_corners[2] = cC;
	// 	m_corners[3] = cD;
	// 	setupSorted();
	// }
	void setGlobalCorners(const emInt cA, const emInt cB, const emInt cC,
			const emInt cD = EMINT_MAX) {
		m_global[0] = cA;
		m_global[1] = cB;
		m_global[2] = cC;
		m_global[3] = cD;
		// setupSorted();
	}
	emInt getCorner(const emInt ii) const {
		assert(ii < m_nCorners);
		return m_corners[ii];
	}
	emInt getGlobalCorner(const emInt ii) const {
		assert(ii < m_nCorners);
		return m_global[ii];
	}
	emInt getSortedGlobal(const emInt ii) const {
		assert(ii < m_nCorners);
		return m_sortedGlobal[ii];
	}
	emInt getRemoteIndices(const emInt ii) const {
		assert(ii < m_nCorners);
		return m_remote[ii];
	}
	emInt getSortedRemoteIndices(const emInt ii) const {
		assert(ii < m_nCorners);
		return m_sortedRemote[ii];
	}
	//virtual void computeParaCoords(const int ii, const int jj,
	//   double st[2]) const = 0;
	emInt getVolElement() const {
		return m_volElem;
	}

	emInt getVolElementType() const {
		return m_volElemType;
	}
	emInt getPartid() const {
		return m_partId;
	}
	void setRemotePartID(const emInt remotePartid_) {
		m_remoteId = remotePartid_;
	}
	void setPartID(const emInt partID_) {
		m_partId = partID_;
	}
	// void setRemoteIndices(const emInt remote[3]){
	// 	for(auto i=0; i<3; i++){
	// 		remoteIndices[i]=remote[i];
	// 	}

	// }
	void setRemoteIndices(const emInt *remote) {
		for (auto i = 0; i < 4; i++) {
			m_remote[i] = *(remote + i);
		}

	}
	emInt getRemoteId() const {
		return m_remoteId;
	}

	bool getGlobalCompare() const {
		return m_globalComparison;
	}
	emInt getNumDivs() const {
		return m_nDivs;
	}
	void setCompare(const bool compare) {
		m_globalComparison = compare;
	}
	//friend MPI_Datatype register_mpi_type(FaceVerts const &);

};
//Notes from Documentation :
// It may not be immediately obvious how this one template serves
//for both saving data to an archive as well as loading data from the archive
// The key is that the & operator is defined as << for output archives and as >> input archives.
class TriFaceVerts: public FaceVerts {
private:

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int /*version*/) {
		ar & boost::serialization::base_object<FaceVerts>(*this);

	}

public:
	TriFaceVerts() {
	}
	; // Dangerous! CHANGE IT LATER ON
	TriFaceVerts(const emInt nDivs, const emInt partID = -1,
			bool globalComparison = false) :
			FaceVerts(nDivs, 3) {
		m_partId = partID;
		m_globalComparison = globalComparison;
	}
	TriFaceVerts(const emInt nDivs, emInt v0, const emInt v1, const emInt v2,
			const emInt type = 0, const emInt elemInd = EMINT_MAX,
			const emInt partID = -1, const emInt remoteId = -1,
			bool globalComparison = false);

	TriFaceVerts(const emInt nDivs, const emInt local[3], const emInt global[3],
			const emInt partid_ = -1, const emInt remoteID = -1,
			const emInt type = 0, const emInt elemInd = EMINT_MAX,
			bool globalComparison = false);

	// TriFaceVerts(const int nDivs,const emInt global[3],
	// const emInt partid_=-1, const emInt remoteID=-1 ,const emInt type=0 ,
	// const emInt elemInd=EMINT_MAX,const bool globalComparison=false);

	TriFaceVerts(const emInt nDivs, const emInt global[3], const emInt partid_ =
			-1, const emInt remoteID = -1, const bool globalComparison = false,
			const emInt type = 0, const emInt elemInd = EMINT_MAX);

	TriFaceVerts(const emInt nDivs, const emInt local[3], const emInt global[3],
			const emInt remoteIndices_[3], const emInt partid_,
			const emInt remoteID, const emInt type = 0, const emInt elemInd =
			EMINT_MAX, bool globalComparison = false);

	//virtual ~TriFaceVerts() {};
	//	void allocVertMemory() {
	//		m_intVerts = new emInt[MAX_DIVS - 2][MAX_DIVS - 2];
	//	}
	//	void freeVertMemory() const {
	//		if (m_intVerts) delete[] m_intVerts;
	//	}
	void computeParaCoords(const emInt ii, const emInt jj, double st[2]) const;
	// virtual void computeParaCoords(const int ii, const int jj,
	// 							   double st[2]) const;
	//virtual void setupSorted();
	void setupSorted();
	void getVertAndST(const emInt ii, const emInt jj, emInt &vert, double st[2],
			const int rotCase = 0) const;
	void getTrueIJ(const emInt ii, const emInt jj, emInt &trueI, emInt &trueJ,
			const int rotCase = 0) const;
	friend bool operator<(const TriFaceVerts &a, const TriFaceVerts &b);
	friend bool operator==(const TriFaceVerts &a, const TriFaceVerts &b);
	friend inline bool compare(const TriFaceVerts &a, const TriFaceVerts &b);
	friend inline MPI_Datatype register_mpi_type(TriFaceVerts const&);
	void setCorners(const emInt cA, const emInt cB, const emInt cC,
			const emInt cD = EMINT_MAX) {
		m_corners[0] = cA;
		m_corners[1] = cB;
		m_corners[2] = cC;
		m_corners[3] = cD;
		setupSorted();
	}
};

struct QuadFaceVerts: public FaceVerts {
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int /*version*/) {
		ar & boost::serialization::base_object<FaceVerts>(*this);

	}

public:
	QuadFaceVerts() {
	}
	; // Dangerous! CHANGE IT LATER ON
	QuadFaceVerts(const emInt nDivs, const emInt partID = -1,
			const emInt _remotePartid = -1, bool globalCompare = false) :
			FaceVerts(nDivs, 4) {
		m_partId = partID;
		m_remoteId = _remotePartid;
		m_globalComparison = globalCompare;
	}
	QuadFaceVerts(const emInt nDivs, const emInt v0, const emInt v1,
			const emInt v2, const emInt v3,

			const emInt type = 0, const emInt elemInd = EMINT_MAX,
			const emInt partID = -1, const emInt remoteID = -1,
			bool globalCompare = false);

	QuadFaceVerts(const emInt nDivs, const emInt local[4],
			const emInt global[4], const emInt partid_ = -1,
			const emInt remoteID = -1, const emInt type = 0,
			const emInt elemInd = EMINT_MAX, bool globalCompare = false);
	// QuadFaceVerts(const int nDivs,const emInt global[4],const emInt partid_=-1, const emInt remoteID=-1
	// ,const emInt type=0 ,const emInt elemInd=EMINT_MAX,
	// bool globalCompare=false);
	QuadFaceVerts(const emInt nDivs, const emInt global[4],
			const emInt partid_ = -1, const emInt remoteID = -1,
			bool globalCompare = false, const emInt type = 0,
			const emInt elemInd = EMINT_MAX);
	QuadFaceVerts(const emInt nDivs, const emInt local[4],
			const emInt global[4], const emInt remotelocal[4],
			const emInt partid_ = -1, const emInt remoteID = -1,
			const emInt type = 0, const emInt elemInd = EMINT_MAX,
			bool globalCompare = false);

	//virtual ~QuadFaceVerts() {}
	void computeParaCoords(const emInt ii, const emInt jj, double st[2]) const;
	// virtual void computeParaCoords(const int ii, const int jj,
	// 							   double st[2]) const;
	//virtual void setupSorted();
	void setupSorted();
	void getVertAndST(const emInt ii, const emInt jj, emInt &vert, double st[2],
			const int rotCase = 0) const;
	void getTrueIJ(const emInt ii, const emInt jj, emInt &trueI, emInt &trueJ,
			const int rotCase = 0) const;
	friend bool operator<(const QuadFaceVerts &a, const QuadFaceVerts &b);
	friend bool operator==(const QuadFaceVerts &a, const QuadFaceVerts &b);
	friend inline MPI_Datatype register_mpi_type(QuadFaceVerts const&);
	void setCorners(const emInt cA, const emInt cB, const emInt cC,
			const emInt cD = EMINT_MAX) {
		m_corners[0] = cA;
		m_corners[1] = cB;
		m_corners[2] = cC;
		m_corners[3] = cD;
		setupSorted();
	}
};

#ifndef USE_ORDERED
namespace std {
template<> struct hash<TriFaceVerts> {
	typedef TriFaceVerts argument_type;
	typedef std::size_t result_type;
	result_type operator()(const argument_type &TFV) const noexcept {
		result_type h0, h1, h2;

		if (TFV.getGlobalCompare() == false) {
			// const result_type h0 = TFV.getSorted(0);
			// const result_type h1 = TFV.getSorted(1);
			// const result_type h2 = TFV.getSorted(2);
			h0 = TFV.getSorted(0);
			h1 = TFV.getSorted(1);
			h2 = TFV.getSorted(2);
			//return (h0 ^ (h1 << 1)) ^ (h2 << 2);
		}
		if (TFV.getGlobalCompare() == true) {
			h0 = TFV.getSortedGlobal(0);
			h1 = TFV.getSortedGlobal(1);
			h2 = TFV.getSortedGlobal(2);

			// const result_type h0 = TFV.getSortedGlobal(0);
			// const result_type h1 = TFV.getSortedGlobal(1);
			// const result_type h2 = TFV.getSortedGlobal(2);
			//return (h0 ^ (h1 << 1)) ^ (h2 << 2);
		}
		return (h0 ^ (h1 << 1)) ^ (h2 << 2);

	}
};

template<> struct hash<QuadFaceVerts> {
	typedef QuadFaceVerts argument_type;
	typedef std::size_t result_type;
	result_type operator()(const argument_type &QFV) const noexcept {
		if (QFV.getGlobalCompare() == false) {
			const result_type h0 = QFV.getSorted(0);
			const result_type h1 = QFV.getSorted(1);
			const result_type h2 = QFV.getSorted(2);
			const result_type h3 = QFV.getSorted(3);
			return h0 ^ (h1 << 1) ^ (h2 << 2) ^ (h3 << 3);
		} else {
			const result_type h0 = QFV.getSortedGlobal(0);
			const result_type h1 = QFV.getSortedGlobal(1);
			const result_type h2 = QFV.getSortedGlobal(2);
			const result_type h3 = QFV.getSortedGlobal(3);
			return h0 ^ (h1 << 1) ^ (h2 << 2) ^ (h3 << 3);

		}

	}
};

template<> struct hash<Edge> {
	typedef Edge argument_type;
	typedef std::size_t result_type;
	result_type operator()(const argument_type &E) const noexcept {
		const result_type h0 = E.getV0();
		const result_type h1 = E.getV1();
		return (h0 ^ (h1 << 1));
	}
};
}
#endif

bool operator==(const TriFaceVerts &a, const TriFaceVerts &b);
bool operator<(const TriFaceVerts &a, const TriFaceVerts &b);
bool operator==(const QuadFaceVerts &a, const QuadFaceVerts &b);
bool operator<(const QuadFaceVerts &a, const QuadFaceVerts &b);

int getTriRotation(const TriFaceVerts &localTri,
		const exa_set<TriFaceVerts> &remote, emInt nDivs);
int getQuadRotation(const QuadFaceVerts &localQuad,
		const exa_set<QuadFaceVerts> &remote, emInt nDivs);
void matchTri(const TriFaceVerts &localTri, const emInt rotation,
		const emInt nDivs, const exa_set<TriFaceVerts> &remoteTriSet,
		std::unordered_map<emInt, emInt> &localRemote);
void matchQuad(const QuadFaceVerts &localQuad, const emInt rotation,
		const emInt nDivs, const exa_set<QuadFaceVerts> &remoteQuadSet,
		std::unordered_map<emInt, emInt> &localRemote);

void findRotationAndMatchTris(const TriFaceVerts &localTri, const emInt nDivs,
		const exa_set<TriFaceVerts> &remoteTriSet,
		std::unordered_map<emInt, emInt> &localRemote);

void findRotationAndMatchQuads(const QuadFaceVerts &localQuad,
		const emInt nDivs, const exa_set<QuadFaceVerts> &remoteQuadSet,
		std::unordered_map<emInt, emInt> &localRemote);

using setTri = std::set<TriFaceVerts>;
using setQuad = std::set<QuadFaceVerts>;
using hashTri = std::unordered_set<TriFaceVerts>;
using hashQuad = std::unordered_set<QuadFaceVerts>;
using vecHashTri = std::vector<hashTri>;
using vecHashQuad = std::vector<hashQuad>;
using vecTri = std::vector<TriFaceVerts>;
using vecQuad = std::vector<QuadFaceVerts>;
using vecVecTri = std::vector<vecTri>;
using vecVecQuad = std::vector<vecQuad>;
using intToVecTri = std::map<int,vecTri>;
using intToVecQuad = std::map<int,vecQuad>;
using TableTri2TableIndex2Index = std::unordered_map< TriFaceVerts , std::unordered_map<emInt,emInt>>;
using TableQuad2TableIndex2Index = std::unordered_map< QuadFaceVerts, std::unordered_map<emInt,emInt>>;

void buildTrisMap(hashTri const &tris, std::map<int, vecTri> &remoteTotris,
		std::set<int> &neighbors);

void buildQuadsMap(hashQuad const &quads, std::map<int, vecQuad> &remoteToquads,
		std::set<int> &neighbors);

BOOST_SERIALIZATION_ASSUME_ABSTRACT(FaceVerts)

void printMultiMap(
		const std::multimap<std::set<emInt>, std::pair<emInt, emInt>> &map);

void printTriFaceVerts(const std::vector<TriFaceVerts> &tris);
void preMatchingPartBdryTris(const emInt numDivs, const setTri &partBdryTris,
		vecVecTri &tris);
void preMatchingPartBdryQuads(const emInt numDivs, const setQuad &partBdryQuads,
		vecVecQuad &quads);
void TestPartFaceMatching(const std::size_t nPart, const vecHashTri &HashTri,
		const vecHashQuad &HashQuad, const vecVecTri &vecTris,
		const vecVecQuad &vecquads);
#endif /* FACEVERTS_H_ */
