#ifndef SIMULATION_H
#define SIMULATION_H

#include "eQ.h"
#include "eQmpi.h"

#include "inputOutput.h"
#include "Strain.h"
#include "fHSL.h"
#include "./abm/eQabm.h"
#include "../diffuclass.h"

class Simulation
{    
protected:


public:
    eQ::mpi                 mpiController;
    std::vector<eQ::mpi>    mpiHSL;

    struct Params
	{
        int                                 argc;
        char                                **argv;
        inputOutput                         *fileIO;
        std::shared_ptr<eQ::data::files_t>    dataFiles;
        std::shared_ptr<eQ::uniformRandomNumber>    zeroOne;
    };
    Simulation::Params params;


    Simulation(MPI_Comm commWorld, const Simulation::Params &);
    ~Simulation();

    void writeHSLFiles();
    void writeDataFiles();
    void computeGridParameters();
    void create_DataGrid();
    void create_HSLgrid();
//    void create_DSOgrid();
    void init_kinetics_DSO();

//    void init_ABM(int);
    void init_ABM(int numSeedCells, std::vector<std::shared_ptr<Strain>> &strains);

	void stepSimulation();
    void stepFinalize();

    void waitMPI();
    void printData();



    std::vector<std::vector<int>>       gridDofs;
    std::vector<std::vector<double>>    gridCoords;

    std::vector<double>                 D11Grid;
    std::vector<double>                 D22Grid;
    std::vector<double>                 D12Grid;
    std::vector<
                std::shared_ptr<
                    eQ::gridFunction<std::shared_ptr<Strain>>>
                > dsoGrid;

	//array of unknown rates;  swept over range in Chen
    double  basalRate[basalRates::numBasalRates];
    //strong, med, weak rates table;
	struct rates  rates;  
	struct dso_parameters pA, pR;

    double          simTime, sim_dt;
    double          computeTimer, diffusionTimer, physicsTimer, waitTimer, petscTimer;

    void            resetTimers()
    {
        computeTimer    = 0.0;
        diffusionTimer  = 0.0;
        physicsTimer    = 0.0;
        waitTimer       = 0.0;
        petscTimer      = 0.0;
    }

    bool            simulateABM     = false;
    bool            createdHSLgrid  = false;
    bool            HSL_signaling   = false;

    eQabm::Params                       paramsABM;
    std::shared_ptr<eQabm>              ABM;
    std::shared_ptr<fenicsInterface>    diffusionSolver;
    std::shared_ptr<diffusionPETSc>     diffusionSolver2;


    double  boundaryWellConcentration;
    int     boundaryUnderFlow;
    double  boundaryDecayRate, wellScaling;

private:
    std::vector<MPI_Group>          mpiGroups;
    std::vector<MPI_Comm>           mpiComms;
    std::vector<std::vector<int>>   mpiRanks;
    size_t whichHSL;

    MPI_Comm                    world, workers, controllerComm;
    MPI_Group                   world_group, worker_group, controller_group;
    int         npes;                // number of PEs
    int         my_PE_num;           // my PE number
    size_t      numHSLGrids;

    size_t      mpiNodesPerDiffusionLayer;


    eQ::diffusionSolver::params fenicsParams;

    bool                isControllerNode, isDiffusionNode;
    eQ::nodeType        globalNodes, globalNodesH, globalNodesW;
    eQ::nodeType        nodesHighData, nodesWideData;
    eQ::nodeType        nodesForChannels;



    void    stepChipmunk();
    void    mpiAssignCommunicators();
    void    initBoundaryWell();
    void    computeBoundaryWell();

};

#endif // SIMULATION_H
