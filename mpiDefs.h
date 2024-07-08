//#include <boost/mpi.hpp>
#include <mpi.h>
#include "exa-defs.h"
#include "Part.h"


inline MPI_Datatype register_mpi_type(Part const&){
	constexpr std::size_t num_members=9; 
	int lengths[num_members];
	for (std::size_t i=0 ; i<num_members; i++){
		lengths[i]= 1; 
	}
	Part dummy(1,1,1,1,1,1,1,1,1);
	MPI_Aint baseadress ; 
	MPI_Aint offsets [num_members];
	MPI_Get_address(&dummy,&baseadress); 
	MPI_Get_address(&dummy.m_xmin,&offsets[0]); 
	MPI_Get_address(&dummy.m_xmax,&offsets[1]); 
	MPI_Get_address(&dummy.m_ymin,&offsets[2]); 
	MPI_Get_address(&dummy.m_ymax,&offsets[3]); 
	MPI_Get_address(&dummy.m_zmin,&offsets[4]); 
	MPI_Get_address(&dummy.m_zmax,&offsets[5]); 
	MPI_Get_address(&dummy.m_first,&offsets[6]); 
	MPI_Get_address(&dummy.m_last,&offsets[7]); 
	MPI_Get_address(&dummy.m_nParts,&offsets[8]); 

	for (std::size_t i=0 ; i<num_members ;i++){
		offsets[i]=MPI_Aint_diff(offsets[i],baseadress);
	}
    MPI_Datatype types[num_members]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
	MPI_DOUBLE,MPI_DOUBLE,MPI_INT32_T,MPI_INT32_T,MPI_INT32_T};
	MPI_Datatype type;
	MPI_Type_create_struct(num_members,lengths,offsets,types,&type);
	MPI_Type_commit(&type);
	return type; 
};

inline MPI_Datatype register_mpi_type(CellPartData const&){
	constexpr std::size_t lengths=3; 
	int block_length[lengths]= {1,1,3};
	//int block_length[lengths]= {1,1};
	MPI_Aint displacements [lengths]; 
	//MPI_Aint base_address; 
	//CellPartData dummy (1,1,1,1,1) ; 
    
    displacements[0]=offsetof(CellPartData,m_index); 
	displacements[1]=offsetof(CellPartData,m_cellType);
	displacements[2]=offsetof(CellPartData,m_coords);

	MPI_Datatype types[3] = { MPI_INT32_T, MPI_INT32_T, MPI_DOUBLE };
	//MPI_Datatype types[2] = { MPI_INT32_T, MPI_INT32_T};
	MPI_Datatype type; 
    MPI_Type_create_struct(lengths,block_length, displacements, types, &type);
    MPI_Type_commit(&type);
	return type; 
}
