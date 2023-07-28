//  Copyright 2019 by Carl Ollivier-Gooch.  The University of British
//  Columbia disclaims all copyright interest in the software ExaMesh.//
//
//  This file is part of ExaMesh.
//
//  ExaMesh is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as
//  published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  ExaMesh is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with ExaMesh.  If not, see <https://www.gnu.org/licenses/>.

/*
 * refine.cxx
 *
 *  Created on: Jul. 10, 2019
 *      Author: cfog
 */

#include <unistd.h>
#include <cstdio>

#include "ExaMesh.h"
#include "CubicMesh.h"
#include "UMesh.h"
#include "mpiImpl.h"

int main(int argc, char* const argv[]) {
	char opt = EOF;
	emInt nDivs = 1;
	emInt nTestParts=2; 
	emInt maxCellsPerPart = 1000000;
	char type[10];
	char infix[10];
	char inFileBaseName[1024];
	char cgnsFileName[1024];
	char outFileName[1024];
	bool isInputCGNS = false, isParallel = false, isMPI=false;

	sprintf(type, "vtk");
	sprintf(infix, "b8");
	sprintf(outFileName, "/dev/null");
	sprintf(inFileBaseName, "/need/a/file/name");
	sprintf(cgnsFileName, "/need/a/file/name");

	while ((opt = getopt(argc, argv, "s:c:i:m:n:o:pt:u:q")) != EOF) {
		switch (opt) {
			case 's':
				sscanf(optarg, "%d", &nTestParts);
				break;
			case 'c':
				sscanf(optarg, "%1023s", cgnsFileName);
				isInputCGNS = true;
				break;
			case 'i':
				sscanf(optarg, "%1023s", inFileBaseName);
				break;
			case 'n':
				sscanf(optarg, "%d", &nDivs);
				break;
			case 'm':
				sscanf(optarg, "%d", &maxCellsPerPart);
				break;
			case 'o':
				sscanf(optarg, "%1023s", outFileName);
				break;
			case 'p':
				isParallel = true;
				break;
			case 't':
				sscanf(optarg, "%9s", type);
				break;
			case 'u':
				sscanf(optarg, "%9s", infix);
				break;
			case 'q':
				isParallel=true; 
				isMPI=true; 
				break;
		}
	}

	if (isInputCGNS) {
#if (HAVE_CGNS == 1)
		CubicMesh CMorig(cgnsFileName);
		if (isParallel){
			if(isMPI){
#ifndef NDEBUG
				ParallelTester* tester= new ParallelTester(); 
				CMorig.TestMPI(nDivs,nTestParts,tester);
				//CMorig.refineForMPI(nDivs,tester);
#endif
				//CMorig.refineForMPI(); 
				//CMorig.TestMPI(nDivs,nTestParts); 
				//CMorig.refineForParallel(nDivs, maxCellsPerPart);
			}
			else{
				CMorig.refineForParallel(nDivs, maxCellsPerPart);
			}
		}else {
			double start = exaTime();
			UMesh UMrefined(CMorig, nDivs);
			double time = exaTime() - start;
			size_t cells = UMrefined.numCells();
			fprintf(stderr, "\nDone serial refinement.\n");
			fprintf(stderr, "CPU time for refinement = %5.2F seconds\n", time);
			fprintf(stderr,
							"                          %5.2F million cells / minute\n",
							(cells / 1000000.) / (time / 60));

//			UMrefined.writeUGridFile("/tmp/junk.b8.ugrid");
			UMrefined.writeVTKFile("/tmp/junk.vtk");
		}
#else
		fprintf(stderr, "Not compiled with CGNS; curved meshes not supported.\n");
		exit(1);
#endif
	}
	else {
		UMesh UMorig(inFileBaseName, type, infix);
		if (isParallel){
			if(isMPI)
			{
#ifndef NDEBUG
				ParallelTester* tester= new ParallelTester(); 
				
				UMorig.TestMPI(nDivs,nTestParts,tester); 
				UMorig.refineForMPI(nDivs,tester);
				
				//mpiImpl MPIRefine (&UMorig); // Is this right ? 
				//MPIRefine.refineMPI(nDivs); 	
				//UMorig.refineMPI(nDivs); 
				//UMorig.TestForMPI(nDivs,nTestParts); 


#endif				
				//UMorig.refineForMPI(nDivs,tester); 
				

				//
				
			}
			//else{
				//UMorig.refineForParallel(nDivs, maxCellsPerPart);
			//}

		}
		if (!isParallel) {
			double start = exaTime();
			UMesh UMrefined(UMorig, nDivs);
			double time = exaTime() - start;
			size_t cells = UMrefined.numCells();
			fprintf(stderr, "\nDone serial refinement.\n");
			fprintf(stderr, "CPU time for refinement = %5.2F seconds\n", time);
			fprintf(stderr,
							"                          %5.2F million cells / minute\n",
							(cells / 1000000.) / (time / 60));
			//UMrefined.writeUGridFile(outFileName);
			//UMrefined.writeVTKFile(outFileName);
		}
	}

	printf("Exiting\n");
	exit(0);
}

