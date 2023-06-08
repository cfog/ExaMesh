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
template <typename T>
inline void SetToVector(const std::unordered_set<T>& sourceSet, 
std::vector<T>& destinationVector) {
    destinationVector.clear();
	// Change to the assignment 
   	for (const auto& element : sourceSet) {
        destinationVector.push_back(element);
    }
}

template <typename T>
inline void vectorToSet(const std::vector<T>& sourceVector, 
std::unordered_set<T>& destinationSet) {
    destinationSet.clear();
  
    for (const auto& element : sourceVector) {
        destinationSet.insert(element);
    }
	assert(destinationSet.size()==sourceVector.size()); 
}

// inline void registerTypes(MPI_Datatype &Tri, MPI_Datatype &Quad ,
// MPI_Datatype &TypevecCPD ,
// MPI_Datatype &TypePart ){
// 	emInt globalTri[3]={10,20,30}; 
// 	emInt globalQuad [4]= {10,20,30,40}; 
// 	QuadFaceVerts dummyQuad (1,globalQuad,1,1,true); 
// 	TriFaceVerts dummyTri(1,globalTri,1,1,false); 
// 	//Tri= register_mpi_type(dummyTri); 
// 	//Quad= register_mpi_type(dummyQuad);
// 	Part DummyPart(1,1,1,1,1,1,1,1,1); 

//  	CellPartData DummyCellPartData(1,1,1.0,1.0,1.0); 
// 	TypePart= register_mpi_type(DummyPart);
// 	TypevecCPD= register_mpi_type(DummyCellPartData); 

// }
// inline MPI_Datatype register_mpi_type(TriFaceVerts const&){
// 	constexpr std::size_t num_members= 18 ; 

// 	emInt number_m_intVerts = 
// 	(MAX_DIVS - 1)*(MAX_DIVS - 1);
// 	//std::cout<<number_m_intVerts<<std::endl; 

// 	emInt number_m_param_st = 
// 	(MAX_DIVS + 1)*(MAX_DIVS + 1)*(2); 

// 	//std::cout<<number_m_param_st<<std::endl; 

// 	emInt number_m_param_uvw = 
// 	(MAX_DIVS + 1)*(MAX_DIVS + 1)*(3); 

// 	//std::cout<<number_m_param_uvw<<std::endl; 



// 	int block_length[num_members]= {4,4,4,4,4,4,12,1,1,
// 	number_m_intVerts,number_m_param_st,
// 	number_m_param_uvw,1,1,1,1,1,1};


	
// 	MPI_Aint displacements [num_members]; 
// 	MPI_Aint base_address; 
	
//     displacements[0]=offsetof(TriFaceVerts,global_corners); 
// 	displacements[1]=offsetof(TriFaceVerts,global_sorted);
// 	displacements[2]=offsetof(TriFaceVerts,remoteIndices);
// 	displacements[3]=offsetof(TriFaceVerts,sortedRemoteIndices);
// 	displacements[4]=offsetof(TriFaceVerts,m_corners);
// 	displacements[5]=offsetof(TriFaceVerts,m_sorted);
// 	displacements[6]=offsetof(TriFaceVerts,m_cornerUVW);
// 	displacements[7]=offsetof(TriFaceVerts,m_nCorners);
// 	displacements[8]=offsetof(TriFaceVerts,m_nDivs);
// 	displacements[9]=offsetof(TriFaceVerts,m_intVerts);
// 	displacements[10]=offsetof(TriFaceVerts,m_param_st);
// 	displacements[11]=offsetof(TriFaceVerts,m_param_uvw);
// 	displacements[12]=offsetof(TriFaceVerts,m_volElem);
// 	displacements[13]=offsetof(TriFaceVerts,m_volElemType);
// 	displacements[14]=offsetof(TriFaceVerts,m_bothSidesDone);
// 	displacements[15]=offsetof(TriFaceVerts,partid);
// 	displacements[16]=offsetof(TriFaceVerts,remotePartid);
// 	displacements[17]=offsetof(TriFaceVerts,m_globalComparison);
	
// 	MPI_Datatype types[num_members] = { MPI_INT32_T, MPI_INT32_T,MPI_INT32_T,
// 	MPI_INT32_T,MPI_INT32_T,MPI_INT32_T,
// 	MPI_DOUBLE,
// 	MPI_INT32_T,MPI_INT32_T,
// 	MPI_INT32_T,
// 	MPI_DOUBLE,MPI_DOUBLE,
// 	MPI_INT32_T,MPI_INT32_T,
// 	MPI::BOOL,
// 	MPI_INT32_T,MPI_INT32_T,
// 	MPI::BOOL};

// 	MPI_Datatype type; 
//     MPI_Type_create_struct(num_members,block_length, displacements, types, &type);
//     MPI_Type_commit(&type);
// 	return type; 
// }
// inline MPI_Datatype register_mpi_type(QuadFaceVerts const&){
// 	constexpr std::size_t num_members= 18 ; 
// 	emInt number_m_intVerts = 
// 	(MAX_DIVS - 1)*(MAX_DIVS - 1);

// 	emInt number_m_param_st = 
// 	(MAX_DIVS + 1)*(MAX_DIVS + 1)*(2); 

// 	emInt number_m_param_uvw = 
// 	(MAX_DIVS + 1)*(MAX_DIVS + 1)*(3); 



// 	int block_length[num_members]= {4,4,4,4,4,4,12,1,1,
// 	number_m_intVerts,number_m_param_st,
// 	number_m_param_uvw,1,1,1,1,1,1};


	
// 	MPI_Aint displacements [num_members]; 
// 	MPI_Aint base_address; 
	
//     displacements[0]=offsetof(QuadFaceVerts,global_corners); 
// 	displacements[1]=offsetof(QuadFaceVerts,global_sorted);
// 	displacements[2]=offsetof(QuadFaceVerts,remoteIndices);
// 	displacements[3]=offsetof(QuadFaceVerts,sortedRemoteIndices);
// 	displacements[4]=offsetof(QuadFaceVerts,m_corners);
// 	displacements[5]=offsetof(QuadFaceVerts,m_sorted);
// 	displacements[6]=offsetof(QuadFaceVerts,m_cornerUVW);
// 	displacements[7]=offsetof(QuadFaceVerts,m_nCorners);
// 	displacements[8]=offsetof(QuadFaceVerts,m_nDivs);
// 	displacements[9]=offsetof(QuadFaceVerts,m_intVerts);
// 	displacements[10]=offsetof(QuadFaceVerts,m_param_st);
// 	displacements[11]=offsetof(QuadFaceVerts,m_param_uvw);
// 	displacements[12]=offsetof(QuadFaceVerts,m_volElem);
// 	displacements[13]=offsetof(QuadFaceVerts,m_volElemType);
// 	displacements[14]=offsetof(QuadFaceVerts,m_bothSidesDone);
// 	displacements[15]=offsetof(QuadFaceVerts,partid);
// 	displacements[16]=offsetof(QuadFaceVerts,remotePartid);
// 	displacements[17]=offsetof(QuadFaceVerts,m_globalComparison);
	
// 	MPI_Datatype types[num_members] = { MPI_INT32_T, MPI_INT32_T,MPI_INT32_T,
// 	MPI_INT32_T,MPI_INT32_T,MPI_INT32_T,
// 	MPI_DOUBLE,
// 	MPI_INT32_T,MPI_INT32_T,
// 	MPI_INT32_T,
// 	MPI_DOUBLE,MPI_DOUBLE,
// 	MPI_INT32_T,MPI_INT32_T,
// 	MPI::BOOL,
// 	MPI_INT32_T,MPI_INT32_T,
// 	MPI::BOOL};

// 	MPI_Datatype type; 
//     MPI_Type_create_struct(num_members,block_length, displacements, types, &type);
//     MPI_Type_commit(&type);
// 	return type; 
// }
// struct foo {
//     int x;
//     int y[2][3];

//     foo(int _x) : x(_x) {
//     }

//     foo(int _x, int _y[2][3]) : x(_x) {
//         for (int i = 0; i < 2; i++) {
//             for (int j = 0; j < 3; j++) {
//                 y[i][j] = _y[i][j];
//             }
//         }
//     }
// 	void print (){
// 		std::cout<<"x: "<<x<< std::endl; 
// 		for (int i = 0; i < 2; i++) {
//             for (int j = 0; j < 3; j++) {
//                std::cout<<"i: "<<i<<" j: "<<j<<" y: "<<y[i][j]<<
// 			   std::endl; 
//             }
//         }
		

// 	}
// };

// inline MPI_Datatype register_mpi_type(foo const&){

// 	constexpr std::size_t lengths=1; 

// 	int block_length[lengths]= {1};
	
// 	MPI_Aint displacements [lengths]; 
// 	MPI_Aint base_address; 

    
//     displacements[0]=offsetof(foo,x); 
// 	//displacements[1]=offsetof(foo,y);
	

// 	MPI_Datatype types[lengths] = { MPI_INT32_T};
	
// 	MPI_Datatype type; 
//     MPI_Type_create_struct(lengths,block_length, displacements, types, &type);
//     MPI_Type_commit(&type);
// 	return type; 
// }