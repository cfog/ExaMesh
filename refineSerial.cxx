/*
 * refineSerial.cxx
 *
 *  Created on: Jul. 4, 2019
 *      Author: cfog
 */

int main(int nArgs, char *args[]) {
	char strInputMeshFileName[FILE_NAME_LEN],
			strOutputMeshFileName[FILE_NAME_LEN], strQualFileName[FILE_NAME_LEN];
	bool gotInputFileName = false, gotOutputFileName = false;
	int charRead = EOF, numLevels = -1;

	while ((charRead = getopt(nArgs, args, "i:o:n:")) != EOF) {
		switch (charRead) {
			case 'i':
				strncpy(strInputMeshFileName, optarg, FILE_NAME_LEN - 1);
				gotInputFileName = true;
				break;
			case 'n':
				sscanf(optarg, "%d", &numLevels);
				break;
			case 'o':
				strncpy(strOutputMeshFileName, optarg, FILE_NAME_LEN - 1);
				gotOutputFileName = true;
				break;
			default:
				printUsage();
				break;
		}
	}
	if (!(gotInputFileName)) {
		printUsage();
	}
	if (!(gotOutputFileName)) {
		sprintf(strOutputMeshFileName, "output");
	}

	// Simplest thing that can possibly be useful:  read a mesh, split it,
	// write it.
	GRUMMPInit("examesh");
	openMessageFile("examesh");

	VolMesh VM(2);
	VM.readFromFile(strInputMeshFileName);
	assert(VM.isMeshValid());
	VM.evaluateQuality();
	writeVTKLegacy(VM, strOutputMeshFileName, "-init");

	double timeBefore = clock() / double(CLOCKS_PER_SEC);
	subdividePartMesh(&VM, numLevels);
	double timeAfter = clock() / double(CLOCKS_PER_SEC);
	double elapsed = timeAfter - timeBefore;

	assert(VM.isMeshValid());
	VM.evaluateQuality();
	// Confirm that all cells are positively oriented, in the sense of having
	// positive volumes at least.
#ifndef NDEBUG
	double totalVolume = 0;
	for (GR_index_t ii = 0; ii < pVM->getNumCells(); ii++) {
		Cell *pC = pVM->getCell(ii);
		if (pC->isDeleted()) continue;
		double vol = pC->calcSize();
//		if (vol < 0) {
//			logMessage(MSG_MANAGER, "Neg volume!  Cell %d, %G\n", ii, vol);
//		}
//		assert(vol > 0);
		totalVolume += vol;
	}
	logMessage(MSG_MANAGER, "Total mesh volume: %f\n", totalVolume);
#endif

	makeFileName(strQualFileName, "%s.qual", strOutputMeshFileName,
								"main() [refinePart.cxx]");
	VM.writeQualityToFile(strQualFileName);

	setlocale(LC_ALL, "");
	GR_index_t totalCells = VM.getNumCells();
	logMessage(
			MSG_GLOBAL,
			"Final mesh has:\n %'15u verts,\n %'15u tets,\n %'15u pyramids,\n %'15u prisms,\n %'15u hexes,\n %'15u cells total\n",
			VM.getNumVerts(), VM.getNumTetCells(), VM.getNumPyrCells(),
			VM.getNumPrismCells(), VM.getNumHexCells(), totalCells);
	logMessage(MSG_GLOBAL, "CPU time for refinement = %5.2F seconds\n", elapsed);
	logMessage(MSG_GLOBAL,
							"                        = %5.2F million cells / minute\n",
							(totalCells / 1000000.) / (elapsed / 60));

	logMessage(MSG_GLOBAL, "Writing output files...\n");

	writeVTKLegacy(VM, strOutputMeshFileName);
	//	writeNative(VM, strOutputMeshFileName);
	exit(0);
}

static void printUsage() {
	fprintf(stderr, "Usage: examesh -i mesh_filename -o output_mesh_filename -n "
					"split_levels\n");
	exit(1);
}




