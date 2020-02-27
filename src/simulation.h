#ifndef SIMULATION_H
#define SIMULATION_H

#include "eQ.h"
#include "./abm/eQabm.h"

#include "Strain.h"
#include "fHSL.h"
#include "inputOutput.h"

#include "../diffuclass.h"

class Simulation
{
public:


	struct params
	{
		int argc;
		char **argv;
		inputOutput *fileIO;
        eQ::dataFiles_t dataFiles;
    };


    Simulation(MPI_Comm commWorld, const struct Simulation::params &);
    ~Simulation();

    Simulation::params Params;

    void writeHSLFiles();
    void writeDataFiles();
    void updateDataFiles();

    void computeGridParameters();
    void create_DataGrid();
    void create_HSLgrid(int argc, char* argv[]);
    void create_DSOgrid();

    void init_ABM(int);
    void init_kinetics_DSO();

	void stepSimulation();
    void updateCells();
    void stepFinalize();

    void waitMPI();
    void writebackGridData();
    void printData();
    void finalizeDataRecording(void);

    eQ::diffusionSolver::params fenicsParams;
//    fenicsInterface::params fenicsParams;

    bool isControllerNode, isDiffusionNode;
    size_t globalNodes, globalNodesH, globalNodesW;
    size_t nodesHighData, nodesWideData;

    size_t nodesForChannels;

    double diffusionTimer, physicsTimer, waitTimer, petscTimer;
    bool HSL_signalingTrue;


    std::vector<
                std::vector<int>
            > gridDofs;
    std::vector<
                std::vector<double>
            > gridCoords;
//    std::vector<
//                std::vector<double>
//            > HSLGrids;

    std::vector<double>D11Grid;
    std::vector<double>D22Grid;
    std::vector<double>D12Grid;
    std::vector<
                std::shared_ptr<eQ::gridFunction<std::shared_ptr<Strain>>>
            > dsoGrid;
//    std::vector<
//                std::shared_ptr<eQ::gridFunction<size_t>>
//            > dof_from_grid;
    std::vector<
                std::vector<std::pair<double,double>>
            > coords;

	//array of unknown rates;  swept over range in Chen
	double  basalRate[static_cast<unsigned int>(basalRates::numBasalRates)];
	//strong, med, weak rates table;
	struct rates  rates;  
	struct dso_parameters pA, pR;

    double simTime, sim_dt;
    bool            simulateABM = false;
    bool            createdHSLgrid = false;

    eQabm::params                   paramsABM;
    std::shared_ptr<eQabm>          ABM;
    std::shared_ptr<fenicsInterface>    diffusionSolver;
    std::shared_ptr<diffusionPETSc> diffusionSolver2;


    double  boundaryWellConcentration;
    int     boundaryUnderFlow;
    double  boundaryDecayRate, wellScaling;

private:

    std::vector<MPI_Group> mpiGroups;
    std::vector<MPI_Comm> mpiComms;
    std::vector<std::vector<int>> mpiRanks;
    size_t whichHSL;

    MPI_Comm world, workers, controllerComm;
    MPI_Group world_group, worker_group, controller_group;
    std::vector<MPI_Request> mpiRequest;
    std::vector<MPI_Request> petscRequest;
    std::vector<MPI_Request> mpiRequestA, mpiRequestB, mpiRequestC;
    int         npes;                // number of PEs
    int         my_PE_num;           // my PE number
    size_t      numHSLGrids;

    size_t      mpiNodesPerDiffusionLayer;

    double  computeTimer;

//    std::vector<int> nodeIds;//vector of fenics worker node IDs

    void stepChipmunk();

};

#endif // SIMULATION_H
