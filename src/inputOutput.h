#ifndef _eqINPUTOUTPUT_H
#define _eqINPUTOUTPUT_H
#include <iostream>
#include <array>
#include <queue>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <algorithm>    // std::replace
#include <sstream>      // std::stringstream

struct initParams{
	unsigned int timeSinceEpoch, randSeed;
	double timeStep;
	double trapDiffusionRate;
	std::array<double,2> gridSizeMicrons;
	std::size_t gridNodesPerMicron;
};

class inputOutput
{
	public:	
	inputOutput() = default;

	int parseInputLine(int argc, char* argv[]);
//    std::string initOutputFiles();
    std::string initOutputFiles(std::string imageFilesRoot);
    long getTimeStamp(){return timeStamp;}

	bool isLocalComputer = true;
    bool isAugsburgCluster = false;
    bool isOpuntiaCluster = false;
    bool isArrayLocal = false;
    size_t  slurmArrayIndex=0;
    size_t  localArrayIndex=0;
	std::string fbase;
    long timeStamp;

    std::string launchData;
    std::string timeString;
    std::string dateString;

private:
	std::stringstream sstream;
	std::ofstream logFile;
	std::string fpath;

};
#endif
