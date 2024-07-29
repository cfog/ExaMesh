#include <unistd.h>
#include <fstream> 
#include <iostream> 
#include <vector>
#include <numeric>

std::string writeParallelJobTemplate(int nodes, int tasks, std::string jobName,
		int numDivs, std::string cgnsOrUmesh, std::string meshName) {
	// std::string fullName = "Jobs/"+ jobName+".sh";
	std::string fullName = jobName + ".sh";
	std::ofstream out(fullName);
	out << "#!/bin/bash" << std::endl;
	out << std::endl;
	out << "#SBATCH --account=def-cfog" << std::endl;
	out << "#SBATCH --time=00:30:00" << std::endl;
	out << "#SBATCH --nodes=" << nodes << std::endl;
	out << "#SBATCH --ntasks-per-node=" << tasks << std::endl;
	out << "#SBATCH --mem=0" << std::endl;
	out << "srun ./refine -" << cgnsOrUmesh << " TestCases/" << meshName
			<< " -n " << numDivs << " -q" << std::endl;
	return fullName;
}

void writeSerialJobTemplate(std::string jobName, int numDivs,
		std::string cgnsOrUmesh, std::string meshName) {

	std::ofstream out(jobName);
	out << "#!/bin/bash" << std::endl;
	out << std::endl;
	out << "#SBATCH --account=def-cfog" << std::endl;
	out << "#SBATCH --time=02:00:00" << std::endl;
	out << "#SBATCH --nodes=" << 1 << std::endl;
	out << "#SBATCH --ntasks-per-node=" << 1 << std::endl;
	out << "#SBATCH --mem=0" << std::endl;
	out << "./refine -" << cgnsOrUmesh << " TestCases/" << meshName << " -n "
			<< numDivs << std::endl;

}

int main(int argc, char *const argv[]) {

	std::string fileName = argv[1];
	std::string cgnsorUmesh = argv[2];
	std::cout << fileName << std::endl;
	std::cout << cgnsorUmesh << std::endl;

	std::ofstream MasterJobSubmiiter("submitAlljobs.sh");

	MasterJobSubmiiter << "#!/bin/bash" << std::endl;

	std::vector<int> array = { 50 };

	std::vector<int> nProcessors = { 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024 };
	// for (int num = 2; num <= 64; num += 2)
	// {
	//     nProcessors .push_back(num);
	// }

	std::vector<std::string> files;

	for (auto i = 0; i < array.size(); i++) {

		for (auto k = 0; k < nProcessors.size(); k++) {
			int task;
			int nodes;
			if (nProcessors[k] <= 32) {
				nodes = 1;
				task = nProcessors[k];
			} else {
				nodes = 8;
				task = nProcessors[k] / nodes;
			}
			// else
			// {
			//     nodes=2 ;
			//     task = nProcessors[k]/2;
			// }

			std::string jobName = fileName + "_N" + std::to_string(array[i])
					+ "_nProc" + std::to_string(nProcessors[k]);
			auto jbName = writeParallelJobTemplate(nodes, task, jobName,
					array[i], cgnsorUmesh, fileName);
			files.emplace_back(jbName);
		}
		std::string jobName = fileName + "_N" + std::to_string(array[i])
				+ "_nProc" + std::to_string(1) + ".sh";
		writeSerialJobTemplate(jobName, array[i], cgnsorUmesh, fileName);
		files.emplace_back(jobName);

	}

	for (auto i = 0; i < files.size(); i++) {
		MasterJobSubmiiter << "sbatch " << files[i] << std::endl;
	}

	return 0;

}

