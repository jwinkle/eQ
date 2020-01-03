
#include "inputOutput.h"

int inputOutput::parseInputLine(int argc, char* argv[])
{
    launchData.clear();
	std::cout<<"argc= "<<argc<<std::endl;
//	if(argc == 1) isLocalComputer = true;
	if(argc > 1)
	{
		std::cout<<"argv[] = ";
		for (int i=0; i<argc; i++)
			std::cout<<argv[i]<<" ";
		std::cout<<std::endl;

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

std::string inputOutput::initOutputFiles(std::string imageFilesRoot)
{
    auto timeSinceEpoch =  time(nullptr);
    timeStamp =  timeSinceEpoch;


    timeString = __TIME__;//only populates at compile time
    std::replace( timeString.begin(), timeString.end(), ':', '_'); // replace all 'x' to 'y'
    dateString = __DATE__;//only populates at compile time
	std::replace( dateString.begin(), dateString.end(), ' ', '_'); // replace all 'x' to 'y'

    std::string uniqueString =
            std::to_string(timeSinceEpoch)
            + "-"
            + timeString
            + "-"
            + dateString;

            if(isAugsburgCluster){
              fbase.assign("/homes/gast/winkle/images/");
              fbase += uniqueString;
              if(slurmArrayIndex > 0)
			    fbase += ("-" + std::to_string(slurmArrayIndex));        
			  fbase += "/";
			}   
			else if(isOpuntiaCluster){    
                fbase.assign("/project/josic/images/");
//                fbase.assign("/home/jjwinkle/images/");
              fbase += uniqueString;
              if(slurmArrayIndex > 0)
			    fbase += ("-" + std::to_string(slurmArrayIndex));      
			  fbase += "/";
			}   
			else
			{
                fbase.assign(imageFilesRoot);
                fbase += uniqueString;
                if(isArrayLocal)
                  fbase += ("-" + std::to_string(localArrayIndex));
              fbase += "/";
			}

	return fbase;
}

	/*
    std::string initOutputFiles(struct initParams &p, std::shared_ptr<jNonlinearVariationalSolver> solver)
    {
		sstream << dateString<<std::endl; 
		sstream <<"time(NULL) = "<<p.timeSinceEpoch<<std::endl; 
				sstream<<"TimeStepSize: "<<TimeStepSize<<std::endl;
				sstream<<"MaximumTimeSteps: "<<MaximumTimeSteps<<std::endl; 
				sstream<<"slurmArrayIndex: "<<slurmArrayIndex<<std::endl;

				sstream<<"wL/S: "<<p.wL<<", "<<p.wS<<std::endl;
				sstream<<"BetaL/S: "<<p.betaL<<", "<<p.betaS<<std::endl;
				sstream<<"alphaL/R: "<<p.alphaL<<", "<<p.alphaR<<std::endl;
				sstream<<"betL/R: "<<p.betL<<", "<<p.betR<<std::endl;
				sstream<<"gammaL/R: "<<p.gammaL<<", "<<p.gammaR<<std::endl;
				sstream<<"deltaL/R: "<<p.deltaL<<", "<<p.deltaR<<std::endl;

				sstream<<"randSeed: "<<p.randSeed<<std::endl;

				sstream<<"MPhi: "<<p.MPhi<<std::endl;
				sstream<<"epsPhi: "<<p.epsPhi<<std::endl;
				sstream<<"msym: "<<p.mPhiSym<<std::endl;
				sstream<<"s0: "<<p.s0<<std::endl;

				sstream<<"MTheta: "<<p.MTheta<<std::endl;
				sstream<<"H: "<<p.H<<std::endl;
				sstream<<"T: "<<p.T<<std::endl;
				sstream<<"KappaTheta: "<<p.KappaTheta<<std::endl;
				sstream<<"MConc: "<<p.MConc<<std::endl;


	    sstream<<"SOLVER PARAMETERS:"<<std::endl;
	         // Parameters p = SolverPhi->parameters("newton_solver");
	         Parameters psolver = solver->parameters("newton_solver");
	          sstream<<
	            "linear_solver: "<<std::string(psolver["linear_solver"]) <<std::endl<<
	            "preconditioner: "<<std::string(psolver["preconditioner"]) <<std::endl<<
	            "maximum_iterations: "<<int(psolver["maximum_iterations"]) <<std::endl<<
	            "relative_tolerance: "<<double(psolver["relative_tolerance"]) <<std::endl<<
	            "absolute_tolerance: "<<double(psolver["absolute_tolerance"]) <<std::endl<<
	            "convergence_criterion: "<<std::string(psolver["convergence_criterion"]) <<std::endl<<
	            // std::string(p["report"]) <<std::endl<<
	            "error_on_nonconvergence: "<<(true == bool(psolver["error_on_nonconvergence"]) ? "true":"false")
	          << std::endl;
	    std::cout << sstream.str();
	  
		fpath.assign(fbase + "logFile.txt");
	    logFile.open(fpath, std::ios::trunc);
	    	logFile << sstream.str();
	    logFile.flush();
	    logFile.close();

	    return fbase;
	}
	*/

	/*
	void finalizeLogFile(double timingData[][5])
	{	
	    logFile.open(fpath, std::ios::app);
	    for (int i=0; i< timingCounter;i++)
	    {
	      for (int j=0; j<5; j++)
	        logFile << timingData[i][j] << ",";
	      logFile<<std::endl;
	    }
	    logFile.flush();
	    logFile.close();
	}
	*/
//}	  
