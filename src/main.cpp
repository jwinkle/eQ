#include <dolfin.h>

#include <version.h>

#include "eQinit.h"

#include "simulation.h"
#include "inputOutput.h"
#include "Strain.h"



//auto-generated branch/hash in "version.h" file (see CMakeLists.txt)
static std::string gitBranch    = std::string(GIT_BRANCH);
static std::string gitHash      = std::string(GIT_COMMIT_HASH);


using event_t       = eQ::simulationTiming::triggerEvent;
using params_t      = eQ::data::parametersType;

namespace eQ {
    eQ::data::parametersType &
    operator<<(eQ::data::parametersType &os, const eQ::boundaryCondition &bc)
    {
        eQ::data::parameters["boundaries"] = bc._bcs;
        return os;
    }
    params_t &
    operator<<(params_t &params, const eQ::simulationTiming &timer)
    {
        for(auto flag : timer.flags)
        {
            params["timers"][flag.first] = flag.second.when;
        }
        return params;
    }

    bool							eQ::data::isControllerNode;
    bool							eQ::data::initializedParameters;
    eQ::data::parametersType        eQ::data::parameters;


    std::vector<std::shared_ptr<event_t>>
                                    eQ::simulationTiming::triggerEvent::list;

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
    else return false;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

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

    eQ::mpi(world, 0) >> eQ::mpi::method::BROADCAST >> timeSinceEpoch;

    std::cout<<"MPI node #"<<my_PE_num<<" after MPI transfer = "<<timeSinceEpoch<<std::endl;

    //may return 0 when not able to detect
    const auto processor_count = std::thread::hardware_concurrency();
    std::cout<<"MPI node #"<<my_PE_num<<" reads: std::thread::hardware_concurrency() = "<<processor_count<<std::endl;
    std::cout<<std::endl;
    MPI_Barrier(world);


    // Install a signal handler
    //todo:  synchronize processes to handle Ctrl-C or use abort path exclusively
    std::signal(SIGINT, signal_handler);
    std::signal(SIGTERM, signal_handler);




    inputOutput fileIO(timeSinceEpoch);//pass synchronized timestamp for data files
    auto checkTest = fileIO.parseInputLine(argc, argv);
    if(0 == checkTest)
    {
            std::cout<<"Testing, localArrayIndex = "
                    << fileIO.localArrayIndex
                    <<std::endl;
            return 0;
    }

    fileIO.initOutputFiles(std::string("./"), gitBranch);//pass path to write relative to root path

//**************************************************************************//
//              SET INITIAL SIMULATION PARAMETERS:
//**************************************************************************//

    eQ::initDefaultParameters();
    eQ::initDiffusionRates();


    //NOTE:  THIS CAN BE OVER-RIDDEN IN assignSimulationParameters()
    size_t numSimulations = 1;

//========================================================================================================================================//
//========================================================================================================================================//

//==========================================================================================//
//                              SIMULATION SEQUENCER:
//==========================================================================================//



    if(eQ::data::isControllerNode)
    {
        std::cout<<std::endl<<"TimeSinceEpoch: "<<time(nullptr)<<std::endl;
        std::cout<<"Beign simulations...total number: "<<numSimulations<<std::endl;
    }


//**************************************************************************//
//                          SEQUENCER LOOP:
//**************************************************************************//
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

        MPI_Barrier(world);



    std::vector<size_t> nodes(npes, 1);
    double sum=0;
    int n = 1e6;
    int nLoops = 1e3;

    eQ::mpi                 mpiController;
    std::vector<eQ::mpi>    mpiHSL;
    size_t numHSLGrids=npes-1;
    std::vector<double> initData{0};
    std::vector<std::vector<double>> data(npes, initData);

    if(eQ::data::isControllerNode)
    {
        mpiHSL.assign(numHSLGrids, mpiController);//init vector by default copy constructing with (uninitialized) mpiController member
        for(size_t i(0); i<numHSLGrids; ++i)
        {
            int thisNode = int(i)+1;//iterate through worker nodes (1,...) but use world comm (0,...)
            mpiHSL[i].init(world, thisNode);
        }
    }
    else
    {
        mpiController.init(world, 0);
    }



    auto timeMinutes = 5;

    if(fileIO.isArrayCluster)    timeMinutes = fileIO.timeMinutes;

    if(eQ::data::isControllerNode) std::cout<<"Time minutes = "<<timeMinutes<<std::endl;
    if(eQ::data::isControllerNode) std::cout<<"MPI_Wtime() = "<<MPI_Wtime()<<std::endl;

    auto simulationTimer = eQ::simulationTiming(timeMinutes * 60);

    simulationTimer.wallStart();
    while(simulationTimer.wallTimer())
    {
        if(signalReceived()) break;

        simulationTimer.tare();
        sum=0;
        for(auto loops(0); loops<nLoops; ++loops)
        {
            //processing
             for (int i = 0; i <= n; i++)
              {
                sum = sum + pow(-1, i)/(2*i+1);
              }
        }
        double thisTime = simulationTimer.timer();

         auto pi = 4 * sum;
//         std::cout << "Using " << n << " terms to calculate pi, we have exstimated pi at " << pi << std::endl;

         if(eQ::data::isControllerNode)
         {
             data[0][0]=pi;
             for(auto &node : mpiHSL)
             {
                 node >> eQ::mpi::method::RECV;
                 node >> data[node.index()];
             }
         }
         else
         {
             auto thisData = std::vector<double>{pi};
             mpiController << eQ::mpi::method::SEND;
             mpiController << thisData;
         }


         MPI_Barrier(world);
         if(eQ::data::isControllerNode)  std::cout<<"=========================="<<std::endl;
         sleep(1);
         MPI_Barrier(world);

         for(size_t i(0); i<npes; ++i)
         {
             if(int(i) == my_PE_num)
             {
                 std::cout<<"Process #"<<i<<"  compute time = "<<thisTime<<std::endl;
             }
             MPI_Barrier(world);
         }


        MPI_Barrier(world);
    }//end main for loop


    if(eQ::data::isControllerNode) std::cout<<"MPI_Wtime() = "<<MPI_Wtime()<<std::endl;

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
