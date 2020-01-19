    #include <dolfin.h>

#include "eQ.h"
#include "simulation.h"
#include "inputOutput.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp> // used here for printing


#include "../fenics/AD1Dss.h"

//#include "../histogram-develop/include/boost/histogram.hpp"


//#include <mshr.h>
//using namespace dolfin;
//using namespace mshr;

//static std::string abortPath = "../abort.txt";
static std::string abortPath = "./abort.txt";
static auto abortFlagBoost = boost::filesystem::path(abortPath);   // p reads clearer than argv[1] in the following code


double       g_lambdaFitness;
double       g_compressionLimit;
unsigned int g_indexValue;
unsigned int g_overComp;


////////////////////////////////////////////////////////////////////////////////
static volatile std::sig_atomic_t gSignalStatus;
void signal_handler(int signal)
{
  gSignalStatus = signal;
}
bool signalReceived()
{  
    if( (SIGINT == gSignalStatus) || (SIGTERM == gSignalStatus) )
    {
        std::cout << "SignalValue: " << gSignalStatus << '\n';
        return true;
    }
    else if(boost::filesystem::exists(abortFlagBoost))
    {
        std::cout<<"file: "<<abortPath<<" exists...terminating..."<<std::endl;
        return true;
    }
    else return false;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class jsonRecording
{
public:
    jsonRecording(int whichNode, int whichSim)
        : nodeID(whichNode), simNumber(whichSim) {}
//    void init(std::shared_ptr<Simulation> pSim, long timeStamp)
    void init(std::shared_ptr<Simulation> pSim, long timeStamp, std::string cTime)
    {
        compileTime = cTime;
        sim=pSim;
        timeSinceEpoch=timeStamp;
        clear();
        jfile["jsonSimVersion"] =  "v00_01";//needs to match  a spec. document
        jfile["timeSinceEpoch"] =  timeStamp;
        jfile["nodeID"]         =  nodeID;
        jfile["simNumber"]      =  simNumber;
    }
    void recordFrame()
    {
        thisFrame.clear();
        getCellData();
        thisFrame["simTime"] = sim->simTime;
        thisFrame["divisions"] = sim->ABM->divisionList;
        thisFrame["cells"] = jcells;
        jframes.push_back(thisFrame);
        sim->ABM->divisionList.clear();
    }
    void getCellData()
    {
        jcells.clear();
        auto icells   = sim->ABM->cellList.begin();
        while(icells != sim->ABM->cellList.end())
        {
            thisCell.clear();
            std::shared_ptr<eColi> cell = (*icells);
            std::vector<double> cellData ={
                cell->getCenter_x()
                , cell->getCenter_y()
                , cell->getAngle()
                , cell->getLengthMicrons()
                , cell->getSpringCompression()
                , eQ::proteinNumberToNanoMolar(cell->strain->getProteinNumber(H),  cell->getLengthMicrons())
                , eQ::proteinNumberToNanoMolar(cell->strain->getProteinNumber(I),  cell->getLengthMicrons())
                , eQ::proteinNumberToNanoMolar(cell->strain->getProteinNumber(L),  cell->getLengthMicrons())
            };
            thisCell["i"] =  cell->getCellID();
            thisCell["d"] =  cellData;
            thisCell["p"] =  (eQ::strainType::ACTIVATOR == cell->Params.strainType) ? 0 : 1;
            jcells.push_back(thisCell);
            ++icells;
        }
    }
    void finalize()
    {
        jfile["frames"] = jframes;
        jfile["parameters"] = eQ::parameters;
//        std::cout<<std::setw(4)<<jfile<<std::endl;
        std::ofstream logFile;
        std::string fname =
                compileTime
                + std::to_string(timeSinceEpoch)
                +"_" + std::to_string(nodeID)
                +"-" + std::to_string(simNumber)
                + "_jsonData" + ".json";
        logFile.open(fname, std::ios::trunc);
//        logFile<<std::setw(4)<<jfile<<std::endl;
        logFile<<jfile<<std::endl;
        logFile.flush();
        logFile.close();
        std::cout<<"_jsonData data written to file: "<<fname<<std::endl;
    }
    void clear()
    {
        jfile.clear();
        jframes.clear();
    }
private:
    int nodeID, simNumber;
    json jfile, jframes, jcells, thisFrame, thisCell;
    std::shared_ptr<Simulation> sim;
    long timeSinceEpoch;
    std::string compileTime;
};


//static class members that must be defined outside of main:
eQ::parametersType      eQ::parameters;//main parameters structure
eQ::dataStringsType     eQ::dataStrings[eQ::dataParameter::NUM_DATAPARAMETERS];//file names for recording data
bool                    eQ::initializedParameters = false;
std::default_random_engine eQ::generator;
std::vector<std::shared_ptr<std::lognormal_distribution<double>>>
    eQ::lndistributions;
eQ::uniformRandomNumber eQ::zeroOne;

std::vector<double> eQ::solutionVector;

void eQ::populateFileStrings()
{
    eQ::dataStrings[eQ::dataParameter::SPRING_COMPRESSION]  =  std::string("comp.pvd");
    eQ::dataStrings[eQ::dataParameter::CELL_ANGLE]          =  std::string("angle.pvd");
    eQ::dataStrings[eQ::dataParameter::CELL_VELOCITY]       =  std::string("vel.pvd");
    eQ::dataStrings[eQ::dataParameter::C4]                  =  std::string("c4int.pvd");
    eQ::dataStrings[eQ::dataParameter::C14]                 =  std::string("c14int.pvd");
    eQ::dataStrings[eQ::dataParameter::FP]                  =  std::string("fp.pvd");
    eQ::dataStrings[eQ::dataParameter::CFP]                 =  std::string("cfp.pvd");
    eQ::dataStrings[eQ::dataParameter::YFP]                 =  std::string("yfp.pvd");
    eQ::dataStrings[eQ::dataParameter::SYN]                 =  std::string("synth.pvd");
    eQ::dataStrings[eQ::dataParameter::LACI]                =  std::string("laci.pvd");
    eQ::dataStrings[eQ::dataParameter::AIIA]                =  std::string("aiia.pvd");
    eQ::dataStrings[eQ::dataParameter::MFP]                 =  std::string("mfp.pvd");
    eQ::dataStrings[eQ::dataParameter::QTENSOR]            =  std::string("Qtensor.pvd");
    eQ::dataStrings[eQ::dataParameter::GRADVEL_ISOTROPIC]   =  std::string("gradV_iso.pvd");
    eQ::dataStrings[eQ::dataParameter::MODULUS_S]   =  std::string("modulusS.pvd");
    eQ::dataStrings[eQ::dataParameter::MODULUS_H]   =  std::string("modulusH.pvd");
    eQ::dataStrings[eQ::dataParameter::MODULUS_R]   =  std::string("modulusR.pvd");
    eQ::dataStrings[eQ::dataParameter::DTENSOR_11]   =  std::string("dtensor11.pvd");
    eQ::dataStrings[eQ::dataParameter::DTENSOR_22]   =  std::string("dtensor22.pvd");
    eQ::dataStrings[eQ::dataParameter::DTENSOR_12]   =  std::string("dtensor12.pvd");
}

eQ::databaseClass eQ::database;

// to run with two mpi nodes:
// /usr/bin/mpirun -n 2 -display-allocation  ./eQ
////////////////////////////////////////////////////////////////////////////////
//                          MAIN()
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    eQ::init();

    //need to pass sequential seed values to each mpi node to ensure non-parallel rn generation from same seed.
    eQ::database.init_default_random_engine();

    inputOutput fileIO;
    auto checkTest = fileIO.parseInputLine(argc, argv);
    if(0 == checkTest)
    {
            std::cout<<"Testing, localArrayIndex = "
                    << fileIO.localArrayIndex
                    <<std::endl;
            return 0;
    }

    fileIO.initOutputFiles("./images/");//pass path to write relative to root path

    //global MPI initialization
    MPI_Comm world = MPI_COMM_WORLD;
    int        npes;                // number of PEs
    int        my_PE_num;           // my PE number
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(world, &my_PE_num);
    MPI_Comm_size(world, &npes);
    if(1 == npes)
    {
        std::cout<<"Only one MPI process launched...no parallelization!."
                <<std::endl;
        //to abort instead:
//        std::cout<<"Only one MPI process launched...no parallelization...exiting.  my_PE_num = "
//                << my_PE_num
//                <<std::endl;
//        MPI_Finalize();
//        return 0;
    }

    bool isController = (0 == my_PE_num);
    if(isController)
    {
        //todo:  need to check if slurm array job and switch on array number to only check once
        // p reads clearer than argv[1] in the following code
        auto p = boost::filesystem::path(abortPath);//path defined at top of this file
        if (boost::filesystem::exists(p))
        {
          std::cout<<"file: "<<abortPath<<" exists. Deleting this file..."<<std::endl;
          boost::filesystem::remove(p);
        }
    }
    MPI_Barrier(world);

    // Install a signal handler
    //todo:  synchronize processes to handle Ctrl-C or use abort path exclusively
    std::signal(SIGINT, signal_handler);
    std::signal(SIGTERM, signal_handler);


//**************************************************************************//
//              SET SIMULATION PARAMETERS:
//**************************************************************************//

    std::shared_ptr<Simulation> simulation;
    struct Simulation::params params;
    params.argc = argc;
    params.argv = argv;
    params.fileIO = &fileIO;

    double simTimeTare;
    int simulationStepsMax, stepsPerMin, stepsPerHour;

    auto setSimulationTimeStep = [&] (double dt)
    {
        eQ::parameters["dt"] = dt;
        stepsPerMin = int(round(1.0/double(eQ::parameters["dt"])));
        stepsPerHour = 60*stepsPerMin;
    };
    setSimulationTimeStep(0.01);//default value;  changed when setting sim. parameters based on dt, below

//========================================================================================================================================//
//  lambda functions for initialization per simulation:
//========================================================================================================================================//
    auto getIndex = [&](size_t index, size_t size1D, size_t size2D)
    {
        size_t imutant, iwt;

        if(index < (size1D*size2D))
        {
            iwt = index/size1D;
            imutant = index%size1D;
            return std::make_pair(imutant, iwt);
        }
        else
        {
            return std::make_pair(size_t(0), size_t(0));
        }
    };
    auto computeSimulationScalings = [&]()
    {
        //convert physical height to simulation height: divide by scaling factor
        eQ::parameters["simulationTrapHeightMicrons"] = size_t(round(
                double(eQ::parameters["physicalTrapHeight_Y_Microns"])
                / double(eQ::parameters["lengthScaling"])));
        eQ::parameters["simulationTrapWidthMicrons"] = size_t(round(
                double(eQ::parameters["physicalTrapWidth_X_Microns"])
                / double(eQ::parameters["lengthScaling"])));

        eQ::parameters["diffusionScaling"] =
                1.0/(double(eQ::parameters["lengthScaling"])*double(eQ::parameters["lengthScaling"]));

    };
    auto computeEffectiveDegradationRates = [&](std::map<std::string, double> &physicalDiffusionRates)
    {

        //use an average cell size of 3um to compute volume fraction ratio
        const double rhoe_by_rhoi = eQ::computeVolumeRatio_ExtraToIntra(3.0);

        //the effective degradation rate due to open-boundary flux is: pe/pi * 8D/h^2 (uses peak value for flux, not mean-field 12D/h^2)
        double gammaE_C4 = rhoe_by_rhoi*(8.0 * physicalDiffusionRates["C4"])
                / (double(eQ::parameters["physicalTrapHeight_Y_Microns"]) * double(eQ::parameters["physicalTrapHeight_Y_Microns"]));
        double gammaE_C14 = rhoe_by_rhoi*(8.0 * physicalDiffusionRates["C14"])
                / (double(eQ::parameters["physicalTrapHeight_Y_Microns"]) * double(eQ::parameters["physicalTrapHeight_Y_Microns"]));
        //membrane diffusion rate
        double gammadC4 = 3.0;
        double gammadC14 = 2.1;
        //total degradation as seen intra-cellularly is "parallel" combination of gammas
        double gammaT_C4 = 1./(1./gammaE_C4 + 1./gammadC4);
        double gammaT_C14 = 1./(1./gammaE_C14 + 1./gammadC14);

        eQ::parameters["gammaT_C4"]    = gammaT_C4;
        eQ::parameters["gammaT_C14"]    = gammaT_C14;
        eQ::parameters["rhoe_by_rhoi"] = rhoe_by_rhoi;

    };

    auto assignSimulationParameters = [&](size_t simNum)
    {
        eQ::parameters["diffusionScaling"] = 1.0;

//        eQ::parameters["dt"]            = 0.005;
        eQ::parameters["dt"]            = 0.01;
//        eQ::parameters["dt"]            = 0.02;//simple sender/receiver`
//        eQ::parameters["dt"]            = 0.05;

        eQ::parameters["physicalTrapHeight_Y_Microns"]    = 100;
//        eQ::parameters["physicalTrapHeight_Y_Microns"]    = 80;
//        eQ::parameters["physicalTrapWidth_X_Microns"]     = 1000;
        eQ::parameters["physicalTrapWidth_X_Microns"]     = 500;
//        eQ::parameters["physicalTrapWidth_X_Microns"]     = 400;
//        eQ::parameters["physicalTrapHeight_Y_Microns"]    = 80;
//        eQ::parameters["physicalTrapWidth_X_Microns"]     = 640;

//        eQ::parameters["lengthScaling"] = 10.0;
//        eQ::parameters["lengthScaling"] = 5.0;
        eQ::parameters["lengthScaling"] = 4.0;

        computeSimulationScalings();

        eQ::parameters["trapChannelLinearFlowRate"] = 100.0 * 60.0;//microns/sec * 60sec/min;
        //resolution of the HSL diffusion grid: #lattice points per micron
        eQ::parameters["nodesPerMicronSignaling"] = 2;
        eQ::parameters["trapType"]      = "NOWALLED";
        eQ::parameters["boundaryType"]  = "DIRICHLET_0";
        //scale factors are relative to "WT" division length of ecoli, defined in ecoli.h:
        eQ::parameters["mutantAspectRatioScale"]       = 1.0;
        eQ::parameters["defaultAspectRatioFactor"]     = 1.0;

        eQ::parameters["numberSeedCells"] = 8;
        eQ::parameters["cellInitType"] = "AB_HALF";

        eQ::parameters["modelType"] = "OFF_LATTICE_ABM";
//        eQ::parameters["modelType"] = "ON_LATTICE_STATIC";



        simulationStepsMax = 1 * stepsPerHour;//default value;  changed when setting sim. parameters, below


//****************************************************************************************
                                //ASPECTRATIO_INVASION
//****************************************************************************************

        eQ::parameters["simType"]       = "ASPECTRATIO_INVASION";
        eQ::parameters["modelType"]       = "OFF_LATTICE_ABM";
        eQ::parameters["nodesPerMicronSignaling"] = 2;

//        eQ::parameters["trapChannelLinearFlowRate"] =    100.0 * 60.0;//microns/sec * 60sec/min;
//        eQ::parameters["trapChannelLinearFlowRate"] = 10.0 * 60.0;//microns/sec * 60sec/min;
        eQ::parameters["trapChannelLinearFlowRate"] = 5.0 * 60.0;//microns/sec * 60sec/min;
        eQ::parameters["boundaryType"] = "DIRICHLET_UPDATE";

//        eQ::parameters["trapType"]      = "TWOWALLED";

//        eQ::parameters["physicalTrapHeight_Y_Microns"]    = 100;
//        eQ::parameters["physicalTrapWidth_X_Microns"]     = 400;

        //        eQ::parameters["physicalTrapHeight_Y_Microns"]    = 200;
//        eQ::parameters["physicalTrapWidth_X_Microns"]     = 500;
//        eQ::parameters["physicalTrapWidth_X_Microns"]     = 1000;
//        eQ::parameters["physicalTrapWidth_X_Microns"]     = 300;


        eQ::parameters["physicalTrapHeight_Y_Microns"]    = 100;
        eQ::parameters["physicalTrapWidth_X_Microns"]     = 500;
        eQ::parameters["lengthScaling"] = 5.0;//150mins
//        eQ::parameters["lengthScaling"] = 4.0;//150mins



        eQ::parameters["divisionNoiseScale"] = 0.05;// = +/- 0.025
//        eQ::parameters["divisionNoiseScale"] = 0.0;


        computeSimulationScalings();//sets simulation trap h,w and diffusion scaling from above


        //SET SIMULATION TIME HERE:
//        setSimulationTimeStep(0.01);//updates parameters and stepsPerHour multiplier
        setSimulationTimeStep(0.05);//updates parameters and stepsPerHour multiplier
//        setSimulationTimeStep(0.1);//updates parameters and stepsPerHour multiplier
//        simulationStepsMax = 5*stepsPerHour/2;//150 mins
//        simulationStepsMax = 10*stepsPerHour/2;//300 mins
//        simulationStepsMax = 10*stepsPerHour;//600 mins
//        simulationStepsMax = 30*stepsPerHour;//1800 mins

//        simulationStepsMax = 60*stepsPerHour;//1800 mins
        simulationStepsMax = 180*stepsPerHour;//1800 mins



        eQ::parameters["hslSignaling"]  = true;
        std::map<std::string, double> physicalDiffusionRates = {
            {"C4", 3.0e4}, {"C14", 1.6e4}
        };
        eQ::parameters["physicalDiffusionRates"] = physicalDiffusionRates;

        eQ::parameters["D_HSL"]			= {
                physicalDiffusionRates["C4"] * double(eQ::parameters["diffusionScaling"])//C4
                , physicalDiffusionRates["C14"] * double(eQ::parameters["diffusionScaling"])//C14
//                , physicalDiffusionRates["C4"] * double(eQ::parameters["diffusionScaling"])//C4
//                , physicalDiffusionRates["C14"] * double(eQ::parameters["diffusionScaling"])//C14
        };
        eQ::parameters["membraneDiffusionRates"] = {//for C4, C14 HSL and E.coli membrane diff. rates, see Pai and You (2009)
                3.0
                , 2.1
        };


        eQ::parameters["AnisotropicDiffusion_Axial"] = 1.0;
        //    eQ::parameters["AnisotropicDiffusion_Axial"] = 0.1;
        //        eQ::parameters["AnisotropicDiffusion_Transverse"] = 0.01;
//        eQ::parameters["AnisotropicDiffusion_Transverse"] = 0.2;
                eQ::parameters["AnisotropicDiffusion_Transverse"] = 1.0;




        computeEffectiveDegradationRates(physicalDiffusionRates);

        double hslPeakValue = 1.0e3;

        eQ::parameters["hslProductionRate_C4"]
                = (hslPeakValue * double(eQ::parameters["gammaT_C4"]));//
        eQ::parameters["hslProductionRate_C14"]
                = (hslPeakValue * double(eQ::parameters["gammaT_C14"]));//



        eQ::parameters["numberSeedCells"] = 1000;
//        eQ::parameters["numberSeedCells"] = 20;
        eQ::parameters["cellInitType"] = "AB_HALF";

        eQ::parameters["mutantAspectRatioScale"]        = 0.9;
//        eQ::parameters["mutantAspectRatioScale"]        = 0.8;
        eQ::parameters["defaultAspectRatioFactor"]      = 1.0;

//        eQ::parameters["aspectRatioThresholdHSL"]       = 800.0;
//        eQ::parameters["aspectRatioThresholdHSL"]       = 500.0;
        eQ::parameters["aspectRatioThresholdHSL"]       = 220.0;
//        eQ::parameters["aspectRatioThresholdHSL"]       = 150.0;

        //convert to simulation units here:
        eQ::parameters["channelLengthMicronsLeft"] = 1000/double(eQ::parameters["lengthScaling"]);
        eQ::parameters["channelLengthMicronsRight"] = 1000/double(eQ::parameters["lengthScaling"]);



//****************************************************************************************
    //DATA RECORDING SETUP:
//****************************************************************************************
        eQ::parameters["nodesPerMicronData"]       = 1;
        if(isController)
        {
            //populate here the data to be recorded:
            std::vector<eQ::dataParameter> dataToRecord = {
                eQ::dataParameter::LACI,
//                eQ::dataParameter::DTENSOR_11,
//                eQ::dataParameter::DTENSOR_22,
//                eQ::dataParameter::DTENSOR_12,

//                eQ::dataParameter::MODULUS_S,
//                eQ::dataParameter::MODULUS_H,
//                eQ::dataParameter::MODULUS_R
////                eQ::dataParameter::SPRING_COMPRESSION,
                eQ::dataParameter::C4,
//                eQ::dataParameter::C4RHL,//transcription factor
//                eQ::dataParameter::RHL_T,//total rhlR
//                eQ::dataParameter::SYN,
                eQ::dataParameter::C14,
////                eQ::dataParameter::MFP,
//                eQ::dataParameter::LACI,
//                eQ::dataParameter::FP,
//                eQ::dataParameter::C4,
//                eQ::dataParameter::C14,
            };
            params.dataFiles = eQ::initDataRecording(dataToRecord);
        }
        params.argc = argc;
        params.argv = argv;

        eQ::parameters["recordingInterval"] = 10;//minutes between snapshots

////////////////////////////////////////////////////////////////////////////////
//              SLURM ARRAY INDEXING
////////////////////////////////////////////////////////////////////////////////
    // double noiseVars[NUM_ARRAY_VALS] = {32.,64.,128.,256.};
    //    double noiseVars[NUM_ARRAY_VALS] = {0.0, 16., 32.,64.};
        std::vector<double> wtRatios = {1.0, 1.25, 1.42, 1.67};
        std::vector<double> mRatios = {1.0, 0.8, 0.7, 0.6};
        // double noiseVars[NUM_ARRAY_VALS] = {0.0, 1024., 512., 256.};
//        if(fileIO.isOpuntiaCluster || fileIO.isAugsburgCluster || (simNum >= 1))
        if(false)
        {
            auto indices = (simNum >= 1) ?  //if simNum is sequenced, assume it has control over indexing (ignore slurm arrayID)
                    getIndex(simNum, mRatios.size(), wtRatios.size())
                  : getIndex(fileIO.slurmArrayIndex, mRatios.size(), wtRatios.size());
            auto im = indices.first;
            auto iwt = indices.second;
    //        double g_noiseScale = 0.;
    //        g_noiseScale = noiseVars[fileIO.slurmArrayIndex];
    //        g_aspectRatio = ratios[fileIO.slurmArrayIndex];
    //          eQ::parameters["mutantAspectRatioScale"] = ratios[fileIO.slurmArrayIndex];
            eQ::parameters["mutantAspectRatioScale"] = mRatios[im];
            eQ::parameters["defaultAspectRatioFactor"] = wtRatios[iwt];
        }
        //else defaults

    };//end assignSimulationParameters() lambda function

//========================================================================================================================================//
//========================================================================================================================================//


    auto timeStart = MPI_Wtime();//get mpi wall time


//    int numSimulations = 250;
    size_t numSimulations = 1;
//    size_t numSimulations = 16;


    if(isController)
    {
        std::cout<<"TimeSinceEpoch: "<<fileIO.getTimeStamp()<<std::endl;
        std::cout<<"Beign simulations...total number: "<<numSimulations<<std::endl;
    }


//========================================================================================================================================//
//          SIMULATION SEQUENCER:
//========================================================================================================================================//
    std::shared_ptr<jsonRecording> jsonRecord;
    for(size_t simulationNumber(0); simulationNumber < numSimulations; simulationNumber++)
    {
        if(signalReceived()) break;

        if(isController)
        {
            std::cout<<"Initializing simulation number: "<<simulationNumber+1<<" of "<<numSimulations<<std::endl;
        }

        simTimeTare = MPI_Wtime();//get mpi wall time

        assignSimulationParameters(simulationNumber);


        if(isController)
        {//output the json parameter data structure:
            std::cout<<"Simulation Parameters:"<<std::endl
                    <<std::setw(4)<<eQ::parameters
                   <<std::endl;
        }

//**************************************************************************//
//                          SIMULATION CREATE:
//**************************************************************************//

    //=============================================================================
    //  1. CREATE SIMULATION CLASS:
    //=============================================================================
    simulation = std::make_shared<Simulation>(MPI_COMM_WORLD, params);

    //=============================================================================
    //  2. CREATE HSL GRID, SET UP DATA STRUCTURES AND TRANSFER MAPPINGS:
    //=============================================================================
    if(true == eQ::parameters["hslSignaling"])
    {
        simulation->create_HSLgrid(argc, argv);

        if(isController)
            std::cout<<"Simulation:  create_HSLgrid completed..."<<std::endl;
    }
    if(isController)
    {
        int nodeStamp = 0;
        if(fileIO.isLocalComputer)
        {
            nodeStamp = int(fileIO.localArrayIndex);
        }
        else if( (fileIO.isOpuntiaCluster) || (fileIO.isAugsburgCluster) )
        {
            nodeStamp = int(fileIO.slurmArrayIndex);
        }
        jsonRecord = std::make_shared<jsonRecording>(nodeStamp, simulationNumber);
        jsonRecord->init(simulation, fileIO.getTimeStamp(), fileIO.timeString);

        if(eQ::parameters["modelType"] == "ON_LATTICE_STATIC")
        {
            //DSO ON-LATTICE MODEL FOR TESING:
            double tare = MPI_Wtime();
                simulation->create_DSOgrid();//lattice for dso tests
            tare = MPI_Wtime() - tare;
            std::cout<<"create_DSOgrid() MPI_Wtime() = "<<tare<<std::endl;
        }
        else
        {
    //=============================================================================
    //  3. INITIALIZE THE ABM:
    //=============================================================================
//            std::cout<<"\tController:  simulation->init_ABM begin..."<<std::endl;
            simulation->init_ABM(int(eQ::parameters["numberSeedCells"]));
//            std::cout<<"simulation->init_ABM completed..."<<std::endl;
        }
    }
    //DATA RECORDING GRID CREATION:
    if(!params.dataFiles.empty())
    {
        simulation->create_DataGrid();
        std::cout<<"create_DataGrid completed..."<<std::endl;
//        simulation->writeDataFiles();//initial t=0 data
//        std::cout<<"Initial data write complete... "<<std::endl;
    }

    double computeTimer = 0.0;
    simulation->diffusionTimer = 0.0;
    MPI_Barrier(world);
//    simulation->writeHSLFiles();//initial t=0 data (will not write if no HSL created); switches on "isDiffusionNode"
//    MPI_Barrier(world);

    bool timerFlagThrown = false;
    std::vector<double> channelHSLTimeSeries;

    if(isController)
    {
        std::cout<<"Starting simulation loop... "<<std::endl;
        std::cout<<"stepsPerMin = "<<stepsPerMin<<std::endl;
        std::cout<<"MAX SIMULATION TIME (mins) = "<<simulationStepsMax/stepsPerMin<<std::endl;
    }

//**************************************************************************//
//                          SIMULATION MAIN LOOP:
//**************************************************************************//

    for(int timeSteps=0; timeSteps<=simulationStepsMax; timeSteps++)
    {
        if(signalReceived()) break;


//        if("SINGLEMUTANT_ASPECTRATIO" == eQ::parameters["simType"])
//            if(checkMutantFixation(simulation->simTime))
//                break;

//        if("STATIC_ASPECTRATIO" == eQ::parameters["simType"])
//        {
//            double sr = simulation->ABM->strainRatio;
//            if((0.0 == sr) || (1.0 == sr) )
//            {
//                std::cout<<"Strain fixation of strain: "<<((1.0 == sr) ? "A" : "B")<<std::endl;
//                break;
//            }
//        }


        if(isController)
        {
            if("ASPECTRATIO_INVASION" == eQ::parameters["simType"])
            {
                if( (250.0 < simulation->simTime) && (!timerFlagThrown) )
                {
                    timerFlagThrown = true;
                    simulation->ABM->aspectRatioInduction = true;
                    std::cout<<"Trigger of aspect ratio change to: "
                            <<eQ::parameters["mutantAspectRatioScale"]
                            <<" at simTime="<<simulation->simTime<<std::endl;
                }
            }

            if("SENDER_RECEIVER" == eQ::parameters["simType"])
            {
//                if(250.0 < simulation->simTime)
//                    eQ::parameters["rhlRValue"] = 1.0e6;
//                else
//                    eQ::parameters["rhlRValue"] = 1.0e5;
//                if(200.0 < simulation->simTime)
//                    eQ::parameters["hslProductionRate_C4"] = 0.0;

            }

            if( ("INDUCED_SENDER_RECEIVER" == eQ::parameters["simType"])
                || ("INDUCED_DYNAMIC_ASPECTRATIO" == eQ::parameters["simType"]) )
            {
//                if( (250.0 < simulation->simTime) && (!timerFlagThrown) )
                if( (600.0 < simulation->simTime) && (!timerFlagThrown) )
//                if(false)
                {
                    timerFlagThrown = true;
                    simulation->ABM->aspectRatioInduction = true;
                    std::cout<<"Trigger of aspect ratio change to: "
                            <<eQ::parameters["mutantAspectRatioScale"]
//                            << simulation->ABM->Params.simData.aspectRatioFactor_A
                            <<" at simTime="<<simulation->simTime<<std::endl;
                }
            }
//            if( (250.0 < simulation->simTime) && (!timerFlagThrown) )
//            {
//                timerFlagThrown = true;
//                simulation->ABM->pressureInductionFlag = true;
//                std::cout<<"Trigger of pressureInductionFlag at simTime="<<simulation->simTime<<std::endl;
//            }
        }
//        if(false == isController)  //only HSL signaling nodes
//        {
//            if( (250.0 < simulation->simTime) && (!timerFlagThrown) )
//            {
//                timerFlagThrown = true;
//                simulation->fenicsDiffusion->boundaryDecayRate *= 0.1;
//                std::cout<<"Trigger of boundaryDecayRate at simTime="<<simulation->simTime<<std::endl;
//            }
//        }


//        std::cout<<"simulation->stepSimulation()"<<std::endl;
        simulation->stepSimulation();//parallel chipmunk + HSL
//        std::cout<<"simulation->updateCells()"<<std::endl;
        simulation->updateCells();

        //NOTE:  DATA XFER TO HSL WORKER NODES IS STILL OPEN HERE;  ALL NODES CONTINUE WITHOUT A BARRIER...


        if(timeSteps%(10*stepsPerMin) == 0)
//        if(timeSteps%(1*stepsPerMin) == 0)
        {
            if(isController)
            {
                if(simulation->simulateABM)
                {
                    simulation->ABM->timeSeriesDataTrigger = true;
                    jsonRecord->recordFrame();
                }
            }
            else
            {
                channelHSLTimeSeries.push_back(simulation->boundaryWellConcentration);
            }
        }

//        if(timeSteps%(1*stepsPerMin) == 0)
//        if(timeSteps%(25*stepsPerMin) == 0)
        if(timeSteps%(10*stepsPerMin) == 0)
        {
            if(isController)
            {
                std::cout<<"step: "<<timeSteps<<" = "<<timeSteps*double(eQ::parameters["dt"]);
                if(simulation->simulateABM)
                {
                    std::cout
                            <<" minutes. Cell count: "
                           <<simulation->ABM->cellCount
                          <<", erased: "<<simulation->ABM->eraseCounter
                         <<";  STRAIN RATIO: "<<simulation->ABM->strainRatio
                        <<";  MaxSpringCompression: "<<simulation->ABM->maxSpringCompression
                    <<";  averagePointsPerCell: "<<simulation->ABM->averagePointsPerCell;
                }
                std::cout<<std::endl;
            }
            for(int i(0); i<npes; i++)
            {
                if(i == my_PE_num)
                {
                    std::cout<<"NODE #: "<<my_PE_num<<": ";
                    if(isController)
                    {
                        std::cout<<simulation->physicsTimer<<", "<<simulation->waitTimer<<std::endl;
                    }
                    else
                    {
                        std::cout<<simulation->diffusionTimer<<", "<<simulation->waitTimer<<std::endl;
                    }
                }
                MPI_Barrier(world);
            }
            simulation->diffusionTimer=0.0;
            simulation->physicsTimer=0.0;
            simulation->waitTimer=0.0;

            simulation->writeHSLFiles();
            if(!params.dataFiles.empty())
                simulation->writeDataFiles();

            if(false)
//                if(isController)
            {
                std::cout<<"angle bin data: ";
                for(auto bin : simulation->ABM->binBuffer)
                    std::cout << bin << ",";
                std::cout<<std::endl;
                    std::cout<<"angle bin data2: ";
                    for(auto bin : simulation->ABM->binBuffer2)
                        std::cout << bin << ",";
                    std::cout<<std::endl;
            }
        }

        simulation->stepFinalize();
        MPI_Barrier(world);
    }//end main for loop

    simulation->finalizeDataRecording();

    if(false)
//        if(!isController)
    {
        std::stringstream sstream;
        std::ofstream logFile;
        for(auto conc : channelHSLTimeSeries)
            sstream << conc << std::endl;

//        logFile.open(fileIO.fbase + "_channelHSL_" + std::to_string(my_PE_num) + ".txt", std::ios::trunc);
        logFile.open("_channelHSL_" + std::to_string(fileIO.slurmArrayIndex) + "-" + std::to_string(my_PE_num)
                     + ".txt", std::ios::trunc);
        logFile << sstream.str();
        logFile.flush();
        logFile.close();
        std::cout<<"_channelHSL_ data written to file..."<<std::endl;
    }

    MPI_Barrier(world);
    if(isController)
    {
        std::cout<<"computeTimer() = "<<computeTimer<<std::endl;
        jsonRecord->finalize();
    }
    else
        std::cout<<"simulation->diffusionTimer() = "<<simulation->diffusionTimer<<std::endl;

    MPI_Barrier(world);
    if(!fileIO.isLocalComputer)
    {
        unsigned int sleepTime = 60;
        std::cout<<"sleeping for: "<<sleepTime<<std::endl;
        sleep(sleepTime);
    }
    jsonRecord.reset();
    simulation.reset();
    }//end simulation sequencer loop

//========================================================================================================================================//
    if("SINGLEMUTANT_ASPECTRATIO" == eQ::parameters["simType"])
    {
        std::string dateString = __TIME__;
        std::replace( dateString.begin(), dateString.end(), ':', '_'); // replace all 'x' to 'y'
        dateString += "-";
        dateString += __DATE__;
        std::replace( dateString.begin(), dateString.end(), ' ', '_'); // replace all 'x' to 'y'
        auto timeSinceEpoch =  time(nullptr);

//        std::stringstream sstream;
//        std::ofstream logFile;
//        for(auto data : mutantData)
//            sstream<<std::get<0>(data)<<","<<std::get<1>(data)<<","<<std::get<2>(data)<<","<<std::get<3>(data)<<std::endl;
//        std::string fname;
//        if(fileIO.isAugsburgCluster)
//            fname = std::to_string(timeSinceEpoch) + "_mutantData_" + std::to_string(fileIO.slurmArrayIndex) + ".txt";
//        else
//            fname = std::to_string(timeSinceEpoch) + "_mutantData_" + std::to_string(my_PE_num) + ".txt";

//        logFile.open(fname, std::ios::trunc);
//        logFile << sstream.str();
//        logFile.flush();
//        logFile.close();
//        std::cout<<"_mutantData_ data written to file: "<<fname<<std::endl;

        auto totalTime = (MPI_Wtime()-timeStart)/60.0;//minutes
        std::cout<<"\nTotal simulation time: "<<totalTime<<" : average time per sim. = "<<double(totalTime/numSimulations)<<std::endl;
    }
//========================================================================================================================================//


    //Cleanup MPI:
    int finalized;
    MPI_Finalized(&finalized);
    std::cout<<"BEGIN MPI_Finalized() = "<<finalized<<std::endl;
    if(0 == finalized)
    {
        MPI_Barrier(world);
        MPI_Finalize();
    }
    MPI_Finalized(&finalized);
    std::cout<<"END MPI_Finalized() = "<<finalized<<std::endl;

 return 0;

}//end main()
////////////////////////////////////////////////////////////////////////////////
//                          END MAIN()
////////////////////////////////////////////////////////////////////////////////
//void stepCells(size_t start, size_t end, Simulation *sim)
//{
//    std::shared_ptr<eQ::gridFunction<std::shared_ptr<Strain>>>
//            AGrid = sim->dsoGrid[0];
//    std::shared_ptr<eQ::gridFunction<std::shared_ptr<Strain>>>
//            RGrid = sim->dsoGrid[1];
//    size_t dofC4, dofC14;

//    for(size_t i(start); i<end; i++)
////        for(size_t i(0); i<sim->globalNodesH; i++)
//        for(size_t j(0); j<sim->globalNodesW; j++)
//        {
//            dofC4 = sim->dof_from_grid[0]->grid[i][j];
//            dofC14 = sim->dof_from_grid[1]->grid[i][j];
//            double &c4 = sim->HSLGrids[0][dofC4];
//            double &c14 = sim->HSLGrids[1][dofC14];
//            AGrid->grid[i][j]->computeProteins(c4, c14, 1.0);
//            RGrid->grid[i][j]->computeProteins(c4, c14, 1.0);
//        }
//}
//void threadTest(double *dptr, size_t N, double *mySum)
//{
//    *mySum = 0.0;
//    double *d = dptr;
////    for(size_t k(0);k<N; k++)//square the loop to kill time
//    // for(size_t j(0);j<N; j++)//square the loop to kill time
//    {
//        d = dptr;//reset pointer to top
//        for(size_t i(0); i<N; i++)
//        {
//            *mySum += log(1.0+exp(*d));
//            d++;
//        }
//    }
//}
