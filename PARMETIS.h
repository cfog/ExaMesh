//#include </home/kebriti/local/include/metis.h>
#include <metis.h>
#include "exa-defs.h"
#include "UMesh.h"
std::vector<std::vector<emInt>> partitionMetis(const ExaMesh *const pEM,
		emInt iParts, std::vector<emInt> &vaicelltopart);
void mesh2MetisGraphs(const ExaMesh *const pEM, idx_t xadj[], idx_t adjncy[],
		idx_t adjwgt[]);

void printPart2Cell(const std::vector<std::vector<emInt>> &part2cell);
