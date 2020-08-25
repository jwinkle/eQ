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
#include <iomanip>      // std::setw()

#include <boost/filesystem.hpp>

#include "eQ.h"

class inputOutput
{
	public:	
    inputOutput(long time);

	int parseInputLine(int argc, char* argv[]);
    std::string initOutputFiles(std::string &imageFilesRoot);
    long getTimeStamp(){return timeStamp;}
    void writeParametersToFile(std::string paramsRoot, size_t simNumber, eQ::data::parametersType &params);
    void setSimulationNumber(size_t simNumber);

	bool isLocalComputer = true;
    bool isAugsburgCluster = false;
    bool isOpuntiaCluster = false;
    bool isArrayCluster = false;
    bool isArrayLocal = false;
    size_t  slurmArrayIndex=0;
    size_t  localArrayIndex=0;
	std::string fbase;
    long timeStamp;

    std::string launchData;
    std::string timeString;
    std::string dateString;
    std::string uniqueString;

private:
	std::stringstream sstream;
	std::ofstream logFile;
	std::string fpath;
    std::string froot;

};
#endif
