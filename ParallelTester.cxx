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
ParallelTester::setInputTri (const vecTriHash& inTris)
{
    for(auto itri=0 ; itri<inTris.size(); itri++)
    {
        triHash Settri; 
        for(auto&iset:inTris[itri])
        {
            int nDivs    =  iset.getNumDivs(); 
            int partId   =  iset.getPartid(); 
            int remoteId =  iset.getRemotePartid(); 

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
ParallelTester::testInputTri (const int rank, const triHash& inTri) const{
    assert(inTri.size()==m_tris[rank].size()); 
   // for(const auto& itri:inTri){
        // How to check whether two unordered set are equal 

   // }
}
ParallelTester::~ParallelTester(){
    // delete m_vecPart; 
    // delete m_vecPartData; 
}