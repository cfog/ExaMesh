#include "exa-defs.h" 
#include "Part.h"
using vecPart               = std::vector<Part>; 
using vecCellPartData       = std::vector<CellPartData>; 

class ParallelTester {
private: 
    vecPart         m_vecPart; 
    vecCellPartData m_vecPartData; 
    vecTriHash      m_tris; 
    vecQuadHash     m_quads; 



public: 
    ParallelTester(){}; 
    ~ParallelTester( ); 
    void setVecpart(const vecPart& inVecPart); 
    void setvecCellPartData(const vecCellPartData& inVecCellPartData); 
    void setInputTri (const vecTriHash& inTris); 

    void testVecPart(const int rank,const Part& inPart) const; 
    void testVecCellPartData(const vecCellPartData& inVecCellPartData)const; 
    void testInputTri (const int rank, const triHash& inTri) const; 


}; 