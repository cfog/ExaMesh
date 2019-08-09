/*
 * examesh.cxx
 *
 *  Created on: Jul. 10, 2019
 *      Author: cfog
 */

#include <unistd.h>

#include "examesh.h"
#include "UMesh.h"

int main(int argc, char* const argv[]) {
	char opt = EOF;
	int nDivs = 1;
	char type[10];
	char infix[10];
	char inFileBaseName[1024];
	char outFileName[1024];

	sprintf(type, "vtk");
	sprintf(infix, "b8");
	sprintf(outFileName, "/dev/null");
	sprintf(inFileBaseName, "/need/a/file/name");

	while ((opt = getopt(argc, argv, "i:n:o:t:u:")) != EOF) {
		switch (opt) {
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

	UMesh UMorig(inFileBaseName, type, infix);
	UMesh UMrefined(UMorig, nDivs);
	UMrefined.writeUGridFile("/tmp/junk.b8.ugrid");
	UMrefined.writeVTKFile("/tmp/junk.vtk");

	exit(0);
}

