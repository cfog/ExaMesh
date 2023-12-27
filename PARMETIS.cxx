#include "PARMETIS.h"
#include <fstream>
void setMetisOptions(idx_t options[])
{
    int ierr; 
    ierr=METIS_SetDefaultOptions(options);
    assert(ierr==METIS_OK);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_NUMBERING] = 0;
    options[METIS_OPTION_CONTIG] = 1;
    // 0: Does not force contiguous partitions.
    // 1: Forces contiguous partitions.

    //options[METIS_OPTION_UFACTOR] = 25;

/*  options[METIS OPTION MINCONN]
    Specifies that the partitioning routines should try to minimize the maximum degree of the subdomain graph,
    i.e., the graph in which each partition is a node, and edges connect subdomains with a shared interface.
    0: Does not explicitly minimize the maximum connectivity.
    1:Explicitly minimize the maximum connectivity. */


/*     Specifies the maximum allowed load imbalance among the partitions. A value of x indicates that the
        allowed load imbalance is (1 + x)/1000. The load imbalance for the jth constraint is defined to be
        maxi (w[j, i])/t[j, i]), where w[j, i] is the fraction of the overall weight of the jth constraint that is as-
        signed to the ith partition and t[j, i] is the desired target weight of the jth constraint for the ith partition
        (i.e., that specified via -tpwgts). For -ptype=rb, the default value is 1 (i.e., load imbalance of 1.001) and for
        -ptype=kway, the default value is 30 (i.e., load imbalance of 1.03). */



    // Options valid for METIS_PartGraphway : 
 /*    METIS_OPTION_OBJTYPE, METIS_OPTION_CTYPE, METIS_OPTION_IPTYPE,
    METIS_OPTION_RTYPE, METIS_OPTION_NO2HOP, METIS_OPTION_NCUTS,
    METIS_OPTION_NITER, METIS_OPTION_UFACTOR, METIS_OPTION_MINCONN,
    METIS_OPTION_CONTIG, METIS_OPTION_SEED, METIS_OPTION_NUMBERING,
    METIS_OPTION_DBGLVL */
}
std::vector<std::vector<idx_t>> 
buildPart2Cell (const idx_t* aicelltopart, emInt numCells,emInt numParts)
{
    std::vector<std::vector<idx_t>> part2cell(numParts);

    for (auto cellID = 0; cellID < numCells; cellID++) 
    {
        idx_t partID = aicelltopart[cellID];
        part2cell[partID].push_back(cellID);
    }
    return part2cell; 
}
std::vector<std::vector<emInt>> partitionMetis(const std::unique_ptr<UMesh> &pEM, emInt iParts)
{
//   int METIS PartGraphKway(idx t *nvtxs, idx t *ncon,    idx t *xadj,   idx t *adjncy,
//                           idx t  *vwgt, idx t *vsize,   idx t *adjwgt, idx t *nparts, real t *tpwgts,
//                           real t ubvec, idx t *options, idx t *objval, idx t *part)

// nvtxs: The number of vertices in the graph -> pEM->nCells()
// ncon : The number of balacing constraints. It should be at least 1. ? Is it one? 
    idx_t options[METIS_NOPTIONS];
    // Pointer Declaration and Initialization 
    idx_t *xadj   = nullptr;
    idx_t *adjncy = nullptr;
    idx_t *adjwgt = nullptr;
    idx_t *aicelltopart = new idx_t[pEM->numCells()];

    idx_t ncon=1; 
    idx_t objval; 
    idx_t npart  = iParts; 

    
    idx_t nCells = pEM->numCells(); 

    xadj    = (idx_t*)calloc((size_t)(pEM->numCells()+1), sizeof(idx_t));
    adjncy  = (idx_t*)calloc((size_t)(pEM->numCells()*MAXADJ), sizeof(idx_t));
    adjwgt  = (idx_t*)calloc((size_t)(pEM->numCells()*MAXADJ), sizeof(idx_t));

    setMetisOptions(options); 
    mesh2MetisGraphs(pEM,xadj,adjncy,adjwgt); 
    METIS_PartGraphKway(&nCells, &ncon, xadj, adjncy,
                        NULL, NULL, adjwgt,&npart, NULL,
                        NULL, options, &objval,aicelltopart); 

    auto part2cell = buildPart2Cell(aicelltopart,pEM->numCells(),npart);

    printPart2Cell(part2cell); 




    free(xadj) ; 
    free(adjncy); 
    free(adjwgt); 
    delete [] aicelltopart;
    return part2cell; 
}; 
void extractPartitions(){}; 
void mesh2MetisGraphs(const std::unique_ptr<UMesh> &pEM, idx_t xadj[], idx_t adjncy[], idx_t adjwgt[])
{
    xadj [0] = 0 ;



    for (emInt icell=0 ; icell< pEM->numCells(); icell++)
    {
        xadj[icell+1] = xadj[icell]; 

        for (int iNeigh = 0 ; iNeigh < pEM->getCellConnSize(icell) ; iNeigh++)
        {
            adjncy[ xadj[icell+1] ] = pEM->getCellConn(icell,iNeigh);
            adjwgt[ xadj[icell+1] ] = 1;	
            xadj[icell+1]++;
        } 

    }




}; 


void calcMetisWeights()
{

}



void printPart2Cell (const std::vector<std::vector<idx_t>>& part2cell) 
{

    std::ofstream outFile("partition.txt");
  //  for (int partID = 0; partID < part2cell.size(); partID++) 
  //  {
    auto partID=0; 
        outFile << "Part " <<partID<< ": ";
        for (const auto& cell : part2cell[partID]) 
        {
            outFile << cell << " ";
        }
        outFile << std::endl;
 //   }

}

