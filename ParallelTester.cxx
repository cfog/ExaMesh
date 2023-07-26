#include "ParallelTester.h"
#include <algorithm>

void 
ParallelTester:: setVecpart(const vecPart& inVecPart)
{
    m_vecPart.assign(inVecPart.begin(), inVecPart.end());
}
void 
ParallelTester::setvecCellPartData(const vecCellPartData& inVecCellPartData)
{   
    m_vecPartData.assign(inVecCellPartData.begin(),inVecCellPartData.end());
}
void 
ParallelTester::testVecPart(const int rank, const Part& inPart) const
{
    assert(m_vecPart[rank]==inPart); 
}
void
ParallelTester::testVecCellPartData(const vecCellPartData& inVecCellPartData) const
{
    assert(m_vecPartData.size()==inVecCellPartData.size()); 
    // since the data is broadcasted no need for rank  
    assert(std::equal(m_vecPartData.begin(),m_vecPartData.end(),inVecCellPartData.begin())); 
    
};
void 
ParallelTester::setInputTri (const vecHashTri& inTris)
{
    for(auto itri=0 ; itri<inTris.size(); itri++)
    {
        hashTri Settri; 
        for(auto&iset:inTris[itri])
        {
            int nDivs    =  iset.getNumDivs(); 
            int partId   =  iset.getPartid(); 
            int remoteId =  iset.getRemoteId(); 

            int global0  =   iset.getGlobalCorner(0); 
            int global1  =   iset.getGlobalCorner(1); 
            int global2  =   iset.getGlobalCorner(2); 
            int global[3]=   {global0,global1,global2}; 
            TriFaceVerts tri(nDivs,global,partId,remoteId,true); 
            Settri.insert(tri); 
        }
        m_tris.push_back(Settri); 


    
    }
    assert(m_tris.size()==inTris.size()); 
    for(auto itri=0 ; itri<inTris.size(); itri++)
    {
  
       assert(inTris[itri].size()==m_tris[itri].size()); 
    }
    
} 
void 
ParallelTester::testInputTri (const int rank, const hashTri& inTri) const{
    assert(inTri.size()==m_tris[rank].size()); 
   // for(const auto& itri:inTri){
        // How to check whether two unordered set are equal 

   // }
}
ParallelTester::~ParallelTester(){
    // delete m_vecPart; 
    // delete m_vecPartData; 
}

void 
ParallelTester::setMatchedTris (const vecMatchedTris& inMatchedTris)
{
    m_matchedTrisAllParts= inMatchedTris; 
}
void
ParallelTester::setMatchedQuads (const vecMatchedQuads& inMatchedQuads)
{
    m_matchedQuadsAllParts= inMatchedQuads;
}

void 
ParallelTester::testMatchedTris     
(const TableTri2TableIndex2Index& matchedTris, const int rank) 
const
{

    auto targetTris = m_matchedTrisAllParts[rank]; 

    for(auto tri:matchedTris)
    {
        
        auto findtri = targetTris.find(tri.first); 
        assert(findtri!=targetTris.end()); 

        auto targetMap  = findtri->second; 
        auto thisTriMap = tri.second; 
        assert(thisTriMap==targetMap);


    }

}

void
ParallelTester:: testMatchedQuads    
(const TableQuad2TableIndex2Index& matchedQuads, const int rank) 
const
{

    auto targetQuads = m_matchedQuadsAllParts[rank];

    for(auto quad:matchedQuads)
    {
        auto findQuad = targetQuads.find(quad.first); 
        assert(findQuad!=targetQuads.end()); 
        auto targetMap = findQuad->second; 
        auto thisQuadMap = quad.second; 
        assert(thisQuadMap==targetMap); 

    }

}