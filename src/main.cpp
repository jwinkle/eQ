#include <dolfin.h>

#include "eQinit.h"

#include "simulation.h"
#include "inputOutput.h"
#include "Strain.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp> // used here for printing


static std::string abortPath = "./abort.txt";
static auto abortFlagBoost = boost::filesystem::path(abortPath);   // p reads clearer than argv[1] in the following code

namespace eQ {
    eQ::data::parametersType &operator<<(eQ::data::parametersType &os, const eQ::boundaryCondition &bc)
    {
        eQ::data::parameters["boundaries"] = bc._bcs;
        return os;
    }
    bool							eQ::data::isControllerNode;
    bool							eQ::data::initializedParameters;
    eQ::data::parametersType        eQ::data::parameters;

    //HSL type=key, {D,d} (diffusion constants in media D and membrane rate d)
    //note: initialized in main
    std::map<std::string, std::vector<double>> eQ::data::physicalDiffusionRates;

}



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
        std::shared_ptr<eColi> cell;


        //TRAVERSE LIST OF CELLS: (uses C++ STL forward list data structure)
        sim->ABM->cellList.beginIteration();
        while(sim->ABM->cellList >> cell)
        {
            thisCell.clear();
            //SENDER_RECEIVER
            std::vector<double> cellData ={
                cell->getCenter_x()
                , cell->getCenter_y()
                , cell->getAngle()
                , cell->getLengthMicrons()
                , cell->getSpringCompression(),
//                , eQ::proteinNumberToNanoMolar(cell->strain->getProteinNumber(Strain::concentrations::FP),  cell->getLengthMicrons())
            };
            std::vector<double> ihsl;
            if(sim->HSL_signaling)
            {
                for(auto &hsl : cell->strain->getHSLvec())
                {
                    cellData.push_back(
                                eQ::Cell::proteinNumberToNanoMolar(hsl,  cell->getLengthMicrons()));
                }
            }

            thisCell["i"] =  cell->getCellID();
            thisCell["d"] =  cellData;
            thisCell["p"] =  (eQ::Cell::strainType::ACTIVATOR == cell->params.strain->getStrainType()) ? 0 : 1;
            jcells.push_back(thisCell);
            ++sim->ABM->cellList;
        }
    }
    void finalize()
    {
        jfile["frames"] = jframes;
        jfile["parameters"] = eQ::data::parameters;
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


// to run with two mpi nodes:
// /usr/bin/mpirun -n 2 -display-allocation  ./eQ
////////////////////////////////////////////////////////////////////////////////
//                          MAIN()
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    //global MPI initialization
    int        npes;                // number of PEs
    int        my_PE_num;           // my PE number
    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(world, &my_PE_num);
    MPI_Comm_size(world, &npes);

    eQ::data::isControllerNode = (0 == my_PE_num);



    sleep(my_PE_num);
    long timeSinceEpoch = time(nullptr);
    std::cout<<"MPI node #"<<my_PE_num<<" initially reads: time(nullptr) = "<<timeSinceEpoch<<std::endl;

    eQ::mpi dataTransferFromController;
    dataTransferFromController.init(world, 0);

        dataTransferFromController >> eQ::mpi::method::BROADCAST >> timeSinceEpoch;

    std::cout<<"MPI node #"<<my_PE_num<<" after MPI transfer = "<<timeSinceEpoch<<std::endl;

    //may return 0 when not able to detect
    const auto processor_count = std::thread::hardware_concurrency();
    std::cout<<"MPI node #"<<my_PE_num<<" reads: std::thread::hardware_concurrency() = "<<processor_count<<std::endl;
    std::cout<<std::endl;
    MPI_Barrier(world);


    inputOutput fileIO(timeSinceEpoch);//pass synchronized timestamp for data files
    auto checkTest = fileIO.parseInputLine(argc, argv);
    if(0 == checkTest)
    {
            std::cout<<"Testing, localArrayIndex = "
                    << fileIO.localArrayIndex
                    <<std::endl;
            return 0;
    }
    fileIO.initOutputFiles("./images/");//pass path to write relative to root path


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

    if(eQ::data::isControllerNode)
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
//              SET INITIAL SIMULATION PARAMETERS:
//**************************************************************************//

    eQ::initDefaultParameters();
    eQ::initDiffusionRates();


    std::shared_ptr<Simulation>     simulation;
    Simulation::Params              params;
    params.argc         =       argc;
    params.argv         =       argv;
    params.fileIO       =       &fileIO;
    //pass seed here for repeatability:
    params.zeroOne      =       std::make_shared<eQ::uniformRandomNumber>();
    params.dataFiles    =       std::make_shared<eQ::data::files_t>();


    eQ::simulationTiming    simulationTimer;

    std::vector<std::shared_ptr<Strain>> strainTypes;

//========================================================================================================================================//
//  lambda functions for initialization per simulation:
//========================================================================================================================================//
    auto setSimulationTimeStep = [&] (double dt)
    {
        eQ::data::parameters["dt"] = dt;
        simulationTimer.dt(dt);
    };
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
        eQ::data::parameters["simulationTrapHeightMicrons"] = size_t(round(
                double(eQ::data::parameters["physicalTrapHeight_Y_Microns"])
                / double(eQ::data::parameters["lengthScaling"])));
        eQ::data::parameters["simulationTrapWidthMicrons"] = size_t(round(
                double(eQ::data::parameters["physicalTrapWidth_X_Microns"])
                / double(eQ::data::parameters["lengthScaling"])));

        eQ::data::parameters["diffusionScaling"] =
                1.0/(double(eQ::data::parameters["lengthScaling"])*double(eQ::data::parameters["lengthScaling"]));


        double channelLengthLeft    = double(eQ::data::parameters["mediaChannelMicronsLeft"])/double(eQ::data::parameters["lengthScaling"]);
        double channelLengthRight    = double(eQ::data::parameters["mediaChannelMicronsRight"])/double(eQ::data::parameters["lengthScaling"]);
        eQ::data::parameters["simulationChannelLengthLeft"] = channelLengthLeft;//channel length is in simulation units
        eQ::data::parameters["simulationChannelLengthRight"] = channelLengthRight;//channel length is in sROBINimulation units
    };
    auto computeEffectiveDegradationRates = [&](std::map<std::string, std::vector<double>> &physicalDiffusionRates)
    {
        //use an average cell size of 3um to compute volume fraction ratio
        const double rhoe_by_rhoi
                = eQ::Cell::computeVolumeRatio_ExtraToIntra(0.75 * eQ::Cell::DEFAULT_DIVISION_LENGTH_MICRONS);

        //the effective degradation rate due to open-boundary flux is: pe/pi * 8D/h^2 (uses peak value for flux, not mean-field 12D/h^2)
        double gammaE_C4 = rhoe_by_rhoi*(8.0 * physicalDiffusionRates["C4"][0])
                / (double(eQ::data::parameters["physicalTrapHeight_Y_Microns"]) * double(eQ::data::parameters["physicalTrapHeight_Y_Microns"]));
        double gammaE_C14 = rhoe_by_rhoi*(8.0 * physicalDiffusionRates["C14"][0])
                / (double(eQ::data::parameters["physicalTrapHeight_Y_Microns"]) * double(eQ::data::parameters["physicalTrapHeight_Y_Microns"]));
        //membrane diffusion rate
        double gammadC4 = physicalDiffusionRates["C4"][1];
        double gammadC14 = physicalDiffusionRates["C14"][1];
        //total degradation as seen intra-cellularly is "parallel" combination of gammas
        double gammaT_C4 = 1./(1./gammaE_C4 + 1./gammadC4);
        double gammaT_C14 = 1./(1./gammaE_C14 + 1./gammadC14);

        eQ::data::parameters["gammaT_C4"]    = gammaT_C4;
        eQ::data::parameters["gammaT_C14"]    = gammaT_C14;
        eQ::data::parameters["rhoe_by_rhoi"] = rhoe_by_rhoi;

    };
    auto setBoundaryConditions = [&](double flowRate)
    {
        eQ::data::parameters["trapChannelLinearFlowRate"] = flowRate;//microns/sec;
        eQ::data::parameters["simulationFlowRate"] = flowRate * 60.0/double(eQ::data::parameters["lengthScaling"]);//microns/sec * 60sec/min;

        using type = eQ::boundaryCondition::type;
        using wall = eQ::boundaryCondition::wall;

        eQ::boundaryCondition bcs;
        if("DIRICHLET_0" == eQ::data::parameters["boundaryType"])
        {
            bcs << type::DIRICHLET_0
                    << wall::ALL;
        }
        else if("MICROFLUIDIC_TRAP" == eQ::data::parameters["boundaryType"])
        {
            // aH' + bH = cS
            if("H_TRAP" == eQ::data::parameters["trapType"])
            {
                //Robin BC for left/right walls: u'/u = 1/length = b/a
                bcs << type::ROBIN
                        << wall::LEFT
                        << wall::RIGHT;

                bcs << type::NEUMANN_0
                        << wall::TOP
                        << wall::BOTTOM;
            }
            else if("TWOWALLED" == eQ::data::parameters["trapType"])
            {
                bcs << type::NEUMANN_0
                        << wall::LEFT
                        << wall::RIGHT;

                bcs << type::CHANNELUPDATE
                        << wall::TOP
                        << wall::BOTTOM;
            }
            else
            {
                //Robin BC for left/right walls: u'/u = 1/length = b/a
                bcs << type::ROBIN
                        << wall::LEFT
                        << wall::RIGHT;

                bcs << type::CHANNELUPDATE
                        << wall::TOP
                        << wall::BOTTOM;
            }
        }

        eQ::data::parameters << bcs;

    };
    auto checkAdvectionDiffusionStability = [&](double D)
    {
        //nodes per micron must be adjust by length scaling to get "physical" node size
        double h = double(eQ::data::parameters["lengthScaling"])/double(eQ::data::parameters["nodesPerMicronSignaling"]);
        double v = double(eQ::data::parameters["trapChannelLinearFlowRate"]) * 60.0;//use physical rate (um/min)
        double dt = double(eQ::data::parameters["dt"]);

        double h_stability = 2.0*D/(v * h);
        double t_stability = h*h/(D*dt);//Crank-N scheme

        size_t timeLoops = ceil(dt/t_stability);
        eQ::data::parameters["channelSolverNumberIterations"] = timeLoops;

        std::cout<<"\t  ADVECTION-DIFFUSION STABILITY (h,t) > 1 => loops for channel solver: "
                <<h_stability<<", "<<t_stability
               <<"  ==>  "<<timeLoops<<std::endl;

        return( (h_stability > 1.0) && (t_stability > 1.0) );
    };
//****************************************************************************************
            //assignSimulationParameters: SENDER_RECEIVER
//****************************************************************************************
    auto assignSimulationParameters = [&](size_t simNum)
    {

        eQ::data::parameters["simType"]       = "SENDER_RECEIVER";
        int numberOfDiffusionNodes      = 1;

//        eQ::data::parameters["simType"]       = "MODULUS_2";
//        int numberOfDiffusionNodes      = 1;

//        eQ::data::parameters["simType"]       = "INDUCED_DYNAMIC_ASPECTRATIO";
//        int numberOfDiffusionNodes      = 0;

//        eQ::data::parameters["simType"]       = "ASPECTRATIO_INVASION";
//        int numberOfDiffusionNodes      = 2;


//****************************************************************************************
                        //GEOMETRY
//****************************************************************************************
        eQ::data::parameters["modelType"]         = "OFF_LATTICE_ABM";

        eQ::data::parameters["trapType"]          = "NOWALLED";
//        eQ::data::parameters["trapType"]      = "TWOWALLED";
//        eQ::data::parameters["trapType"]      = "H_TRAP";

//        eQ::data::parameters["boundaryType"]  = "DIRICHLET_0";
//        eQ::data::parameters["boundaryType"]  = "DIRICHLET_UPDATE";
        eQ::data::parameters["boundaryType"]      = "MICROFLUIDIC_TRAP";

//        eQ::data::parameters["PETSC_SIMULATION"] =  true;
        eQ::data::parameters["PETSC_SIMULATION"]  =  false;

//        double trapFlowRate = 0.0;//um/sec
//        double trapFlowRate = 1.0;//um/sec
//        double trapFlowRate = 5.0;//um/sec
//        double trapFlowRate = 10.0;//um/sec
//        double trapFlowRate = 25.0;//um/sec
//        double trapFlowRate = 50.0;//um/sec
        double trapFlowRate = 100.0;//um/sec
//        double trapFlowRate = 150.0;
//        double trapFlowRate = 250.0;

//        eQ::data::parameters["lengthScaling"] = 2.0;//150mins
        eQ::data::parameters["lengthScaling"] = 5.0;//150mins
//        eQ::data::parameters["lengthScaling"] = 4.0;//150mins


        eQ::data::parameters["physicalTrapHeight_Y_Microns"]    = 100;
        eQ::data::parameters["physicalTrapWidth_X_Microns"]     = 500;



        eQ::data::parameters["mediaChannelMicronsLeft"] = 100;
        eQ::data::parameters["mediaChannelMicronsRight"] = 100;

        setBoundaryConditions(trapFlowRate);
        computeSimulationScalings();//sets simulation trap h,w and diffusion scaling from above

//****************************************************************************************
                    //SIMULATION TIME:
//****************************************************************************************
//        simulationStepsMax = 30*stepsPerHour;//1800 mins
    setSimulationTimeStep(0.1);//resets the timer object
//    setSimulationTimeStep(0.01);

    simulationTimer.setSimulationTimeMinutes(10);
//    simulationTimer.setSimulationTimeMinutes(100);
//    simulationTimer.setSimulationTimeMinutes(300);
//    simulationTimer.setSimulationTimeMinutes(600);

//****************************************************************************************
                        //DIFFUSION
//****************************************************************************************

        eQ::data::parameters["hslSignaling"]  = (numberOfDiffusionNodes > 0);

        //resolution of the HSL diffusion grid: #lattice points per micron
        eQ::data::parameters["nodesPerMicronSignaling"]   = 2;
        eQ::data::parameters["physicalDiffusionRates"]    = eQ::data::physicalDiffusionRates;
        eQ::data::parameters["membraneDiffusionRates"] = {//for C4, C14 HSL and E.coli membrane diff. rates, see Pai and You (2009)
                eQ::data::physicalDiffusionRates["C4"][1],  //3.0 min^-1
                eQ::data::physicalDiffusionRates["C14"][1]//1.7 min^-1
                };
//                eQ::data::physicalDiffusionRates["C4"][1],  //3.0 min^-1
//                eQ::data::physicalDiffusionRates["C14"][1]};//1.7 min^-1

        eQ::data::parameters["D_HSL"].clear();
        //will need to add to this for other HSL:
        for(int i(0); i<numberOfDiffusionNodes; ++i)
        {
            if(0 == i)
                eQ::data::parameters["D_HSL"].push_back(eQ::data::physicalDiffusionRates["C4"][0] * double(eQ::data::parameters["diffusionScaling"]));
            else if(1 == i)
                eQ::data::parameters["D_HSL"].push_back(eQ::data::physicalDiffusionRates["C14"][0] * double(eQ::data::parameters["diffusionScaling"]));
//            else if(2 == i)
//                eQ::data::parameters["D_HSL"].push_back(eQ::data::physicalDiffusionRates["C4"][0] * double(eQ::data::parameters["diffusionScaling"]));
//            else if(3 == i)
//                eQ::data::parameters["D_HSL"].push_back(eQ::data::physicalDiffusionRates["C14"][0] * double(eQ::data::parameters["diffusionScaling"]));
        }

        if(bool(eQ::data::parameters["hslSignaling"]))
        {
            computeEffectiveDegradationRates(eQ::data::physicalDiffusionRates);
            checkAdvectionDiffusionStability(eQ::data::physicalDiffusionRates["C4"][0]);//send max. D, sets channel loop iterations

            //target max HSL in bulk:
            double hslPeakValue = 1.0e4;

            eQ::data::parameters["hslProductionRate_C4"]
                    = (hslPeakValue * double(eQ::data::parameters["gammaT_C4"]));//
            eQ::data::parameters["hslProductionRate_C14"]
                    = (hslPeakValue * double(eQ::data::parameters["gammaT_C14"]));//
        }


        eQ::data::parameters["AnisotropicDiffusion_Axial"]        = 1.0;
        eQ::data::parameters["AnisotropicDiffusion_Transverse"]   = 1.0;
        //    eQ::data::parameters["AnisotropicDiffusion_Axial"] = 0.1;
        //        eQ::data::parameters["AnisotropicDiffusion_Transverse"] = 0.01;
//        eQ::data::parameters["AnisotropicDiffusion_Transverse"] = 0.2;



//****************************************************************************************
                        //INITIAL CELLS
//****************************************************************************************
        //scale factors are relative to "WT" division length of ecoli, defined in ecoli.h:
        //DEFAULTS:  (may be changed below)
        eQ::data::parameters["mutantAspectRatioScale"]       = 1.0;
        eQ::data::parameters["defaultAspectRatioFactor"]     = 1.0;

        auto dt         = double(eQ::data::parameters["dt"]);
        auto npm        = double(eQ::data::parameters["nodesPerMicronSignaling"]);
        auto numHSL     = size_t(eQ::data::parameters["D_HSL"].size());


        eQ::data::parameters["numberSeedCells"] = 16;
//        eQ::data::parameters["numberSeedCells"] = 32;

//        eQ::data::parameters["cellInitType"] = "RANDOM";
        eQ::data::parameters["cellInitType"] = "BANDED";


        Strain::Params strainA {eQ::Cell::strainType::ACTIVATOR,
                    dt, npm, numHSL, eQ::Cell::DEFAULT_PROMOTER_DELAY_TIME_MINUTES, nullptr};//nullptr is for cell data, not yet created the cell
        Strain::Params strainB {eQ::Cell::strainType::REPRESSOR,
                    dt, npm, numHSL, eQ::Cell::DEFAULT_PROMOTER_DELAY_TIME_MINUTES, nullptr};//nullptr is for cell data, not yet created the cell

        //ASPECT RATIO:
//        strainTypes.push_back(std::make_shared<aspectRatioInvasionStrain>(strainA));
//        strainTypes.push_back(std::make_shared<aspectRatioInvasionStrain>(strainB));
//        eQ::data::parameters["mutantAspectRatioScale"]       = 0.6;

        //SENDER-RECEIVER:
        strainTypes.push_back(std::make_shared<sendRecvStrain>(strainA));
        strainTypes.push_back(std::make_shared<sendRecvStrain>(strainB));

        //MODULUS:
//        strainTypes.push_back(std::make_shared<MODULUSmodule>(dt, npm, numHSL));



        eQ::data::parameters["divisionNoiseScale"] = 0.1;// = +/- 0.05
//        eQ::data::parameters["divisionNoiseScale"] = 0.05;// = +/- 0.025
    //        eQ::data::parameters["divisionNoiseScale"] = 0.0;


        eQ::data::parameters["rhlRValue"] = 1.0e5;
        eQ::data::parameters["senderScale"] = 64.0;


//****************************************************************************************
                        //DATA RECORDING SETUP:
//****************************************************************************************
        eQ::data::parameters["nodesPerMicronData"]       = 1;
        eQ::data::parameters["recordingInterval"] = 10;//minutes between snapshots
//        if(false)
        if(eQ::data::isControllerNode)
        {
            //populate here the data to be recorded:
                auto n = size_t(eQ::data::parameters["nodesPerMicronData"]);
                auto y = size_t(eQ::data::parameters["simulationTrapHeightMicrons"]);
                auto x = size_t(eQ::data::parameters["simulationTrapWidthMicrons"]);

                params.dataFiles->push_back(eQ::data::record
                    {eQ::dataParameterType::HSL, Strain::hsl::C4, std::string("c4int.pvd"),
                     std::make_shared<eQ::data::tensor>(n,y,x,eQ::data::tensor::rank::SCALAR)});

//                params.dataFiles->push_back(eQ::gridData
//                    {eQ::dataParameterType::HSL, Strain::hsl::C14, std::string("c14int.pvd"),
//                     std::make_shared<eQ::tensorData>(n,y,x,eQ::tensorData::rank::SCALAR)});

//                params.dataFiles->push_back(eQ::gridData
//                    {eQ::dataParameterType::PROTEIN, Strain::concentrations::CFP, std::string("cfp.pvd"),
//                     std::make_shared<eQ::tensorData>(n,y,x,eQ::tensorData::rank::SCALAR)});
//                params.dataFiles->push_back(eQ::gridData
//                    {eQ::dataParameterType::PROTEIN, Strain::concentrations::H, std::string("h.pvd"),
//                     std::make_shared<eQ::tensorData>(n,y,x,eQ::tensorData::rank::SCALAR)});
//                params.dataFiles->push_back(eQ::gridData
//                    {eQ::dataParameterType::PROTEIN, Strain::concentrations::L, std::string("laci.pvd"),
//                     std::make_shared<eQ::tensorData>(n,y,x,eQ::tensorData::rank::SCALAR)});

        }
    };//end assignSimulationParameters() lambda function

//========================================================================================================================================//
//========================================================================================================================================//

//==========================================================================================//
//                              SIMULATION SEQUENCER:
//==========================================================================================//

    simulationTimer.tare();

//    size_t numSimulations = 1;
//    size_t numSimulations = 2;
    size_t numSimulations = 4;
    //    size_t numSimulations = 16;

    if(eQ::data::isControllerNode)
    {
        std::cout<<std::endl<<"TimeSinceEpoch: "<<time(nullptr)<<std::endl;
        std::cout<<"Beign simulations...total number: "<<numSimulations<<std::endl;
    }


//**************************************************************************//
//                          SEQUENCER LOOP:
//**************************************************************************//
    std::shared_ptr<jsonRecording> jsonRecord;
    for(size_t simulationNumber(0); simulationNumber < numSimulations; simulationNumber++)
    {
        if(signalReceived()) break;

        if(eQ::data::isControllerNode)
        {
            std::cout<<std::endl;
            std::cout<<std::endl;
            std::cout<<"Initializing simulation number: "<<simulationNumber+1<<" of "<<numSimulations<<std::endl;
        }

        std::cout<<std::endl;
        assignSimulationParameters(simulationNumber);

        if(eQ::data::isControllerNode)
        {//output the json parameter data structure:
            std::cout<<std::endl
                    <<"Simulation Parameters:"<<std::endl
                    <<std::setw(4)<<eQ::data::parameters
                      <<std::endl<<std::endl;

            fileIO.writeParametersToFile("./", simulationNumber, eQ::data::parameters);
        }

        fileIO.setSimulationNumber(simulationNumber);

        MPI_Barrier(world);


    //=============================================================================
    //  1. CREATE SIMULATION CLASS: (creates ABM)
    //=============================================================================

    simulation = std::make_shared<Simulation>(MPI_COMM_WORLD, params);

    //=============================================================================
    //  2. CREATE HSL GRID, SET UP DATA STRUCTURES AND TRANSFER MAPPINGS:
    //=============================================================================
    if(true == eQ::data::parameters["hslSignaling"])
    {
        simulation->create_HSLgrid();

        if(eQ::data::isControllerNode)
            std::cout<<"Simulation:  create_HSLgrid completed..."<<std::endl<<std::endl;
    }
    if(eQ::data::isControllerNode)
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

        if(eQ::data::parameters["modelType"] == "ON_LATTICE_STATIC")
        {
//            //DSO ON-LATTICE MODEL FOR TESING:
//            double tare = simulationTimer.getTime();
//            simulation->create_DSOgrid();//lattice for dso tests
//            std::cout<<"create_DSOgrid() MPI_Wtime() = "
//                    <<simulationTimer.getTime() - tare
//                    <<std::endl;
            std::cout<<"create_DSOgrid() MPI_Wtime() DISABLED in main()"
                    <<std::endl;
        }
    //=============================================================================
    //  3. INITIALIZE THE ABM:
    //=============================================================================
        else
        {
            simulation->init_ABM(int(eQ::data::parameters["numberSeedCells"]), strainTypes);
        }
        //DATA RECORDING GRID CREATION:
        if(!params.dataFiles->empty())
        {
            simulation->create_DataGrid();
            std::cout<<"create_DataGrid completed..."<<std::endl<<std::endl;
        }
    }


    //=============================================================================
    //  3. INITIALIZE TIMER FLAGS:
    //=============================================================================

    simulation->resetTimers();

    if( ("INDUCED_SENDER_RECEIVER" == eQ::data::parameters["simType"])
        || ("INDUCED_DYNAMIC_ASPECTRATIO" == eQ::data::parameters["simType"]) )
    {
//        simulationTimer.setTimerFlag("Aspect ratio induction", 600);
        simulationTimer.setTimerFlag("Aspect ratio induction", 300);
    }

    auto checkTimerFlags = [&]()
    {
        if("SENDER_RECEIVER" == eQ::data::parameters["simType"])
        {
//                if(250.0 < simulation->simTime)
//                    eQ::data::parameters["rhlRValue"] = 1.0e6;
//                else
//                    eQ::data::parameters["rhlRValue"] = 1.0e5;
//                if(200.0 < simulation->simTime)
//                    eQ::data::parameters["hslProductionRate_C4"] = 0.0;

        }
        if( ("INDUCED_SENDER_RECEIVER" == eQ::data::parameters["simType"])
            || ("INDUCED_DYNAMIC_ASPECTRATIO" == eQ::data::parameters["simType"]) )
        {
            if(simulationTimer.flagThrown("Aspect ratio induction"))
            {
                simulationTimer.flagIgnore("Aspect ratio induction");
                simulation->ABM->aspectRatioInduction = true;
                std::cout<<std::endl<<std::endl;
                std::cout<<"Trigger of aspect ratio change to: "
                        <<eQ::data::parameters["mutantAspectRatioScale"]
//                            << simulation->ABM->Params.simData.aspectRatioFactor_A
                        <<" at simTime="<<simulation->simTime<<std::endl;
                std::cout<<std::endl<<std::endl;
            }
        }
    };



    if(eQ::data::isControllerNode)
    {
        std::cout<<std::endl<<"\t Starting simulation loop... "<<std::endl;
        std::cout<<"\t stepsPerMin = "<<simulationTimer.stepsPerMin<<std::endl;
        std::cout<<"\t MAX SIMULATION TIME (mins) = "
                <<simulationTimer.getSimulationTimeMinutes()
                <<std::endl<<std::endl;
    }
    auto displayDataStep = [&]()
    {
        if(eQ::data::isControllerNode)
        {
            std::cout<<"step: "
                    <<simulationTimer.steps()<<" = "
                   <<simulationTimer.simTime();
            if(simulation->simulateABM)
            {
                std::cout
                        <<" minutes. Cell count: "
                       <<simulation->ABM->cellList.cellCount()
                      <<", erased: "<<simulation->ABM->eraseCounter
                     <<";  STRAIN RATIO: "<<simulation->ABM->strainRatio
                    <<";  MaxSpringCompression: "<<simulation->ABM->maxSpringCompression
                <<";  averagePointsPerCell: "<<simulation->ABM->averagePointsPerCell;
            }
            std::cout<<std::endl;
        }
        MPI_Barrier(world);
        std::cout<<std::flush;
        for(int i(0); i<npes; i++)
        {
            MPI_Barrier(world);
            if(i == my_PE_num)
            {
                std::cout<<"\t NODE #: "<<my_PE_num<<": \t";
                if(eQ::data::isControllerNode)
                {
                    std::cout<<simulation->physicsTimer<<", \t"<<simulation->waitTimer<<std::endl;
                }
                else
                {
                    std::cout<<simulation->diffusionTimer<<", \t"<<simulation->waitTimer<<std::endl;
//                    std::cout<<"\t PETSC TIMER: "<<simulation->petscTimer<<", \t"<<simulation->waitTimer<<std::endl;
                }
            }
        }
        std::cout<<std::flush;
        MPI_Barrier(world);

    };
    auto recordFrame = [&]()
    {
        if(eQ::data::isControllerNode)
        {
            if(simulation->simulateABM)
            {
                jsonRecord->recordFrame();
            }
        }
    };
//**************************************************************************//
//                          SIMULATION MAIN LOOP:
//**************************************************************************//
    displayDataStep();
    recordFrame();
    simulation->writeHSLFiles();
    simulation->writeDataFiles();
    MPI_Barrier(world);

    while(simulationTimer.stepTimer())
    {
        if(signalReceived()) break;

        if(eQ::data::isControllerNode)
        {
            checkTimerFlags();
        }

        simulation->stepSimulation();//parallel chipmunk + HSL

        //NOTE:  DATA XFER BACK TO HSL WORKER NODES IS STILL OPEN HERE;
        //ALL NODES CONTINUE WITHOUT A BARRIER...

//        if(timeSteps%(10*stepsPerMin) == 0)
        if(simulationTimer.periodicTimeMinutes(10))
        {
            recordFrame();
        }
        MPI_Barrier(world);

        if(simulationTimer.periodicTimeMinutes(10))
        {
            displayDataStep();

            simulation->resetTimers();
            simulation->writeHSLFiles();
            simulation->writeDataFiles();

        }

        simulation->stepFinalize();
        MPI_Barrier(world);
    }//end main for loop


    MPI_Barrier(world);
    if(eQ::data::isControllerNode)
    {
        jsonRecord->finalize();
    }

    MPI_Barrier(world);
    if(!fileIO.isLocalComputer)
    {
//        unsigned int sleepTime = 60;
        unsigned int sleepTime = 30;
        std::cout<<"sleeping for: "<<sleepTime<<std::endl;
        sleep(sleepTime);
    }
    jsonRecord.reset();
    simulation.reset();

//    sleep(1);
    MPI_Barrier(world);

    }//end simulation sequencer loop

//========================================================================================================================================//
//========================================================================================================================================//


    sleep(1);
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
