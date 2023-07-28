//#include "exa-defs"
#include "CubicMesh.h"
#include "UMesh.h"
#include "exa-defs.h"

#ifndef SRC_MPIIMPL_H_
#define SRC_MPIIMPL_H_

class mpiImpl 
{
private:

    const ExaMesh*   mCoarseMesh; // How to get pointer out of this ? 

public: 
    mpiImpl(const ExaMesh* const inCoarseMesh): mCoarseMesh(inCoarseMesh){};
    void refineMPI(const int numDivs);  
    ~mpiImpl(){}; 
  
};
#endif // Do I need this ? 