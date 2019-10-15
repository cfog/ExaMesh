/*
 * examesh.cxx
 *
 *  Created on: Jul. 10, 2019
 *      Author: cfog
 */

#include <unistd.h>
#include <cstdio>

#include "examesh.h"
#include "CubicMesh.h"
#include "UMesh.h"

int main(int argc, char* const argv[]) {
	char opt = EOF;
	int nDivs = 1;
	char type[10];
	char infix[10];
	char inFileBaseName[1024];
	char cgnsFileName[1024];
	char outFileName[1024];
	bool isInputCGNS = false;

	sprintf(type, "vtk");
	sprintf(infix, "b8");
	sprintf(outFileName, "/dev/null");
	sprintf(inFileBaseName, "/need/a/file/name");
	sprintf(cgnsFileName, "/need/a/file/name");

	while ((opt = getopt(argc, argv, "c:i:n:o:t:u:")) != EOF) {
		switch (opt) {
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
			case 'o':
				sscanf(optarg, "%1023s", outFileName);
				break;
			case 't':
				sscanf(optarg, "%9s", type);
				break;
			case 'u':
				sscanf(optarg, "%9s", infix);
				break;
		}
	}

	if (isInputCGNS) {
		CubicMesh CMorig(cgnsFileName);
		UMesh UMrefined(CMorig, nDivs);
		UMrefined.writeUGridFile("/tmp/junk.b8.ugrid");
		UMrefined.writeVTKFile("/tmp/junk.vtk");
	}
	else {
		UMesh UMorig(inFileBaseName, type, infix);
		UMesh UMrefined(UMorig, nDivs);
		UMrefined.writeUGridFile("/tmp/junk.b8.ugrid");
		UMrefined.writeVTKFile("/tmp/junk.vtk");
	}

	exit(0);
}

