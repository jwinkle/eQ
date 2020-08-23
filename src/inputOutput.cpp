
#include "inputOutput.h"

void inputOutput::writeParametersToFile(std::string paramsRoot, size_t simNumber, eQ::data::parametersType &params)
{
    std::ofstream logFile;

    std::string fname =
            paramsRoot
            + uniqueString;

    if(isOpuntiaCluster || isAugsburgCluster)
        fname += ("-" + std::to_string(slurmArrayIndex));
    if(isArrayLocal)
        fname += ("-" + std::to_string(localArrayIndex));


    fname += "_" + std::to_string(simNumber) + "-"
            + "parameters.json";

    logFile.open(fname, std::ios::trunc);
        logFile<<std::setw(4)<<params<<std::endl;
//    logFile<<jfile<<std::endl;
    logFile.flush();
    logFile.close();
    std::cout<<"\t parameters data written to file: "<<fname<<std::endl;
    std::cout<<std::endl;
}

int inputOutput::parseInputLine(int argc, char* argv[])
{
    launchData.clear();
	if(argc > 1)
	{
        if(std::string(argv[1]) == std::string("cluster.math.uni-augsburg.de"))
        {
          isAugsburgCluster = true;
          isLocalComputer = false;
        }
        else if(std::string(argv[1]) == std::string("opuntia.cacds.uh.edu"))
        {
          isOpuntiaCluster = true;
          isLocalComputer = false;
        }
        else if(std::string(argv[1]) == std::string("test"))
        {
            std::cout<<"TESTING ACKNOWLEDGED"<<std::endl;
            return 0;//test return value
        }
	}
    if( (argc > 2) && (isOpuntiaCluster || isAugsburgCluster) )
	{
        isArrayCluster = true;
        slurmArrayIndex = size_t(atoi(argv[2]));
		std::cout<<"SLURM_ARRAY_TASK_ID = "<<slurmArrayIndex<<std::endl;
        if(argc > 3)
        {
           launchData = std::string(argv[3]);
        }
	}	
    else if(argc > 2)
    {
        if(std::string(argv[1]) == std::string("eQarray"))
        {
            localArrayIndex = size_t(atoi(argv[2]));
            isArrayLocal = true;
            return 1;
        }
    }
    return 1;//normal return value
}

void inputOutput::setSimulationNumber(size_t simNumber)
{
    if(isAugsburgCluster){
      fbase.assign("/homes/gast/winkle/images/");
      fbase += uniqueString;
        fbase += ("-" + std::to_string(slurmArrayIndex));
      fbase += ("_" + std::to_string(simNumber));
      fbase += "/";
    }
    else if(isOpuntiaCluster){
//        fbase.assign("./images/");
        fbase.assign("/project/josic/winkle/job/images/");
//                fbase.assign("/home/jjwinkle/images/");
      fbase += uniqueString;
        fbase += ("-" + std::to_string(slurmArrayIndex));
      fbase += ("_" + std::to_string(simNumber));
      fbase += "/";
    }
    else
    {
        fbase.assign(imageFilesRoot);
        fbase += uniqueString;
        if(isArrayLocal)
          fbase += ("-" + std::to_string(localArrayIndex));
        fbase += ("_" + std::to_string(simNumber));
      fbase += "/";
    }

}
std::string inputOutput::initOutputFiles(std::string root)
{

    imageFilesRoot = root;
    timeString = __TIME__;//only populates at compile time
    std::replace( timeString.begin(), timeString.end(), ':', '_'); // replace all 'x' to 'y'
    dateString = __DATE__;//only populates at compile time
	std::replace( dateString.begin(), dateString.end(), ' ', '_'); // replace all 'x' to 'y'

    uniqueString =
            std::to_string(timeStamp)
            + "-"
            + timeString
            + "-"
            + dateString;

    setSimulationNumber(0);
	return fbase;
}
