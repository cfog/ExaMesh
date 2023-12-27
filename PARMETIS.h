#include </home/kebriti/local/include/metis.h>

#if (HAVE_CGNS == 1)
#include "cgnslib.h"
#endif
//#include "ExaMesh.h"
#include "UMesh.h"
std::vector<std::vector<emInt>> partitionMetis(const std::unique_ptr<UMesh> &pEM, emInt iParts); 
void extractPartitions(); 
void mesh2MetisGraphs(const std::unique_ptr<UMesh> &pEM,
idx_t xadj[], idx_t adjncy[], idx_t adjwgt[]);

void printPart2Cell (const std::vector<std::vector<idx_t>>& part2cell); 