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
#include "resultGenerator.cxx"
#include <chrono>

int main(int argc, char *const argv[]) {
	double startAppTime = exaTime();
	char opt = EOF;
	emInt nDivs = 1;
	int nTestParts = 2;
	emInt maxCellsPerPart = 1000000;
	std::string baseFileName("/need/a/file/name");
	std::string fileInfix("b8");
	std::string fileSuffix("unknown");
	std::string outFileName("/dev/null");
	bool isInputCGNS = false, isParallel = false, isMPI = false;
	bool meshScanning = false;

	double satrtScanningTime;
	double scanningTime;

	double wallTime;

	while ((opt = getopt(argc, argv, "g:s:c:i:m:n:o:pt:u:q")) != EOF) {
		switch (opt) {
		case 'g':
			meshScanning = true;
			break;
		case 's':
			sscanf(optarg, "%d", &nTestParts);
			break;
		case 'c':
			baseFileName = optarg;
			isInputCGNS = true;
			fileSuffix = "cgns";
			break;
		case 'i':
			baseFileName = optarg;
			break;
		case 'n':
			sscanf(optarg, "%d", &nDivs);
			break;
		case 'm':
			sscanf(optarg, "%d", &maxCellsPerPart);
			break;
		case 'o':
			outFileName = optarg;
			break;
		case 'p':
			isParallel = true;
			isMPI = false;
			break;
		case 't':
			fileSuffix = optarg;
			break;
		case 'u':
			fileInfix = optarg;
			break;
		case 'q':
			isParallel = true;
			isMPI = true;
			break;
		}
	}

	size_t lastSlashPos = baseFileName.find_last_of('/');
	std::string mshName = baseFileName.substr(lastSlashPos + 1);

	if (isMPI) {
		refineForMPI(baseFileName, fileSuffix, fileInfix, nDivs, mshName);
	} else {
		if (isInputCGNS) {
#if (HAVE_CGNS == 1)
			CubicMesh CMorig(baseFileName);
			double start = exaTime();
			auto refined = CMorig.subdivideMesh(nDivs);
			double time = exaTime() - start;
			size_t cells = refined->numCells();

			fprintf(stderr, "\nDone serial refinement.\n");
			fprintf(stderr, "CPU time for refinement = %5.2F seconds\n", time);
			fprintf(stderr,
					"                          %5.2F million cells / minute\n",
					(cells / 1000000.) / (time / 60));

			//			UMrefined.writeUGridFile("/tmp/junk.b8.ugrid");
			//			UMrefined.writeVTKFile("/tmp/junk.vtk");

#else
		fprintf(stderr, "Not compiled with CGNS; curved meshes not supported.\n");
		exit(1);
#endif
		}
		else {
			UMesh UMorig(baseFileName, fileSuffix, fileInfix);

			double start = exaTime();
			auto refined = UMorig.subdivideMesh(nDivs);
			double time = exaTime() - start;
			size_t cells = refined->numCells();

			fprintf(stderr, "CPU time for refinement = %5.2F seconds\n", time);
			fprintf(stderr,
					"                          %5.2F million cells / minute\n",
					(cells / 1000000.) / (time / 60));
			//UMrefined.writeUGridFile(outFileName);
			//UMrefined.writeVTKFile(outFileName);

		}
		printf("Exiting\n");
		exit(0);
	}

}

