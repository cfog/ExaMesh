#include "exa-defs.h" 
#include "Part.h"
using vecPart               =  std::vector<Part>; 
using vecCellPartData       =  std::vector<CellPartData>; 
using vecMatchedTris        =  std::vector<TableTri2TableIndex2Index>;
using vecMatchedQuads       =  std::vector<TableQuad2TableIndex2Index>;

class ParallelTester {
private: 
    vecPart         m_vecPart; 
    vecCellPartData m_vecPartData; 
    vecHashTri      m_tris; 
    vecHashQuad     m_quads; 
    vecMatchedTris  m_matchedTrisAllParts; 
	vecMatchedQuads m_matchedQuadsAllParts;



public: 
    ParallelTester(){}; 
    ~ParallelTester( ); 
    void setVecpart          (const vecPart& inVecPart); 
    void setvecCellPartData  (const vecCellPartData& inVecCellPartData); 
    void setInputTri         (const vecHashTri& inTris); 
    void setMatchedTris      (const vecMatchedTris& inMatchedTris); 
    void setMatchedQuads     (const vecMatchedQuads& inMatchedQuads);

    void testVecPart         (const int rank,const Part& inPart) const; 
    void testVecCellPartData (const vecCellPartData& inVecCellPartData)const; 
    void testInputTri        (const int rank, const hashTri& inTri) const; 

    void testMatchedTris     (const TableTri2TableIndex2Index& matchedTris  , const int rank) const;

    void testMatchedQuads    (const TableQuad2TableIndex2Index& matchedQuads, const int rank) const;


}; 