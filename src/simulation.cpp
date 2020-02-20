#include "simulation.h"

//static const double linearFlowRate = 100.0 * 60.0;//microns/sec * 60sec/min
static const double linearFlowRate = 10.0 * 60.0;//microns/sec * 60sec/min

Simulation::Simulation(MPI_Comm commWorld, const struct Simulation::params &initParams)
    : Params(initParams), world(commWorld)
{
    // 1.  CREATE MPI COMMUNICATORS FOR CONTROLLER AND "WORKDERS" (DIFFUSION NODES)
    // 2.  INITIALIZE TIMERS, HSL GRID#, #NODES FOR DATA RECORDING

//    int pargc = 2;
//    char pargv1 [] = "./x";
//    char pargv2 [] = "-log_view";
//    char *pargv[] = {pargv1, pargv2};
//    dolfin::SubSystemsManager::singleton().init_petsc(pargc, pargv);

    //Simulation constructor is called by all nodes;  setup MPI partitions first:
    //**************************************************************************//
    //                  MPI INITIALIZATION:
    //**************************************************************************//
    MPI_Comm_rank(world, &my_PE_num);
    MPI_Comm_size(world, &npes);
    MPI_Barrier(world);

    //Generate workers MPI group and communicator
    numHSLGrids = size_t(eQ::parameters["D_HSL"].size());//size of diffusion vector is # of grids (known globally)

    if( (numHSLGrids > 0) && (size_t(npes-1) % numHSLGrids == 0) )
    {//each layer gets same number of mpi processing elements:
        mpiNodesPerDiffusionLayer = size_t(npes-1) / numHSLGrids;
        HSL_signalingTrue = true;
        if(0 == my_PE_num)
        {
            std::cout<<"numHSLGrids = "<<numHSLGrids<<" and mpiNodesPerDiffusionLayer = "<<mpiNodesPerDiffusionLayer<<std::endl;
        }
    }
    else
    {
        std::cout<<"Error: Diffusion list size does not divide number of MPI worker nodes!"<<std::endl;
        HSL_signalingTrue = false;
        mpiNodesPerDiffusionLayer = 0;
    }

    if(1==npes)
    {
        std::cout<<"Simulation class initializing with 1 pe."<<std::endl;
        workers = world;
        controllerComm = world;
        isControllerNode = true;
        //verify this is correct to do:
//        isDiffusionNode = true;
        isDiffusionNode = false;
        whichHSLNode=0;
    }
    else
    {//SAME # PEs (>=1) FOR EACH DIFFUSION LAYER (ALWAYS ONE FOR CONTROLLER LAYER):
        //ITERATE THROUGH EACH DIFFUSION LAYER RANK (1...N) FOR N+1 TOTAL RANKS (INCLUDING CONTROLLER RANK=0)
        int thisRank=0;
        for(size_t i(0); i<numHSLGrids; ++i)
        {
            mpiRanks.push_back(std::vector<int>());//creates mpiRanks[i] vector
            for(size_t j(0); j<mpiNodesPerDiffusionLayer; ++j)
            {
                mpiRanks[i].push_back(++thisRank);//starts at rank=1 for grid=0
                if(my_PE_num == thisRank)
                {
                    whichHSLNode = i;//index from 0; defines which HSL layer for each mpi process
                    std::cout<<"Node "<<my_PE_num<<" belongs to HSL layer "<<whichHSLNode<<std::endl;
                }
            }
        }

        int ranks[1] = {0};

        //create two MPI comms: one with server only (the controller), the other with all but server (the workers):
        MPI_Comm_group(world, &world_group);//extract the group from the comm
        MPI_Group_excl(world_group, 1, ranks, &worker_group);//create a new group, which excludes the server
        MPI_Group_incl(world_group, 1, ranks, &controller_group);//create a new group, only the server

            MPI_Comm_create(world, worker_group, &workers);//create a new comm from new group
            MPI_Comm_create(world, controller_group, &controllerComm);//create a new comm from new group

            //now form sub-comms for each HSL grid to pass to the diffusion solver
            for(size_t i(0); i<numHSLGrids; ++i)
            {
                MPI_Group thisGroup;
                MPI_Comm thisComm;
                MPI_Group_incl(world_group, int(mpiRanks[i].size()), mpiRanks[i].data(), &thisGroup);//create a new group with HSL layer[i] world ranks
                MPI_Comm_create(world, thisGroup, &thisComm);//create a new comm from new group (uses world ranks)
                mpiComms.push_back(MPI_Comm(thisComm));
                MPI_Group_free(&thisGroup);
            }

                MPI_Group_free(&worker_group);//free the groups
                MPI_Group_free(&controller_group);
                MPI_Group_free(&world_group);

        //deterimine who is controller, who is worker and verify # workers = number of grids requested
        isControllerNode = (0 == my_PE_num);
        isDiffusionNode = !isControllerNode;
    }
    //result: MPI nodes partitioned into controllerComm (root node), workers (diffusion solver nodes)

    diffusionTimer = 0.0;
    physicsTimer = 0.0;
    waitTimer = 0.0;
    simTime = 0.0;


        MPI_Barrier(world);

    //compute number of nodes for data recording (uses different resolution of nodes/micron than for HSL signaling):
    nodesHighData = size_t(eQ::parameters["simulationTrapHeightMicrons"])*size_t(eQ::parameters["nodesPerMicronData"])+ 1;
    nodesWideData = size_t(eQ::parameters["simulationTrapWidthMicrons"])*size_t(eQ::parameters["nodesPerMicronData"])+ 1;

    if(isControllerNode)
    {

        std::cout << "HSL_signalingTrue: "<<numHSLGrids<<" grids."<<std::endl;
        std::cout << "\tmpiNodesPerDiffusionLayer: "<<mpiNodesPerDiffusionLayer<<std::endl;

        int thisGrid=0;
        for(auto grid: mpiRanks)
        {
            std::cout<<"grid: "<<thisGrid<<"  ";
            for(auto rank:grid)
                std::cout<<rank<<" ";
            std::cout<<std::endl;
        }
        for(auto comm: mpiComms)
        {
            std::cout<<"comm: "<<comm<<"  ";
        }

        std::cout << "Data Grid Recording: "<<Params.dataFiles.size()<<" grids."<<std::endl;
        std::cout << "MPI_COMM_WORLD has "<<npes<<" nodes"<< std::endl;
        std::cout << "FENICS VERSION:" << std::endl;
            std::cout << "\t dolfin_version(): " << dolfin::dolfin_version() << std::endl;
            std::cout << "\t ufc_signature(): " << dolfin::ufc_signature() << std::endl;
            std::cout << "\t git_commit_hash(): " << dolfin::git_commit_hash() << std::endl;
        std::cout << "argc= "<<Params.argc<<std::endl;
        std::cout<<"argv[] = ";
        for (int i=0; i<Params.argc; i++)
            std::cout<<Params.argv[i]<<" ";
        std::cout<<std::endl;

        char procName[256];
        int resultLength;
        MPI_Get_processor_name(procName, &resultLength);
        std::cout<<"MPI_Get_processor_name: "<<procName<<std::endl;

        double tick = MPI_Wtick();
        std::cout<<"MPI_Wtime has MPI_Wtick: "<<tick<<" seconds"<<std::endl;
        std::cout<<std::endl;
        std::cout<<std::endl;
    }
    MPI_Barrier(world);
}
void Simulation::computeGridParameters()
{
    globalNodesH = size_t(eQ::parameters["simulationTrapHeightMicrons"])*size_t(eQ::parameters["nodesPerMicronSignaling"]) + 1;
    globalNodesW = size_t(eQ::parameters["simulationTrapWidthMicrons"])*size_t(eQ::parameters["nodesPerMicronSignaling"]) + 1;
    globalNodes = globalNodesH * globalNodesW;

    fenicsParams.filePath           = Params.fileIO->fbase;
    fenicsParams.dt                = double(eQ::parameters["dt"]);
    fenicsParams.trapWidthMicrons  = double(eQ::parameters["simulationTrapWidthMicrons"]);
    fenicsParams.trapHeightMicrons = double(eQ::parameters["simulationTrapHeightMicrons"]);

    auto uleft = double(eQ::parameters["channelLengthMicronsLeft"]);
    auto uright = double(eQ::parameters["channelLengthMicronsRight"]);
    nodesForChannels = unsigned(
                ceil( (double(eQ::parameters["simulationTrapWidthMicrons"]) + uleft + uright) * double(eQ::parameters["nodesPerMicronSignaling"]) ));
    nodesForChannels++;//always one more vertex for each dimension vs. # elements
}
void Simulation::create_HSLgrid(int argc, char* argv[])
{//all nodes (controller and workers) call this function

    MPI_Barrier(world);

    //DETERMINE NUMBER OF GRID POINTS:
    computeGridParameters();

    //SET INTERNAL NODE ID FOR CONTROLLER:
//    int workerNodeID = -1;//this is over-written below for all diffusion (worker) nodes; thus -1 is ID of controller

    if(isControllerNode)
    {
        std::cout<<"Creating trap (hxw): "
                <<eQ::parameters["simulationTrapHeightMicrons"]<<" x "
               <<eQ::parameters["simulationTrapWidthMicrons"]
              <<" um^2 ==> "
             <<globalNodes<<" global grid nodes"
            <<std::endl;
    }



//=============================================================================
    //1.  DIFFUSION NODES CREATE DIFFUSION CLASS AND TRANSFER MAPPINGS
//=============================================================================
    if(isDiffusionNode)
    {//only workers compute fenics:
        //TODO:  need to switch on the trap type for BC
        //create the fenics interface class:
    //**************************************************************
    //              CREATE FENICS CLASS:
    //**************************************************************
        //CREATE MAIN FENICS CLASS:
        //create diffusion interface class instance (subclass of template defined in fHSL.h)
//        diffusionSolver = std::make_shared<fenics_HSL<hslD::FunctionSpace, hslD::LinearForm, hslD::BilinearForm>>();//create instance
        diffusionSolver = std::make_shared<fenicsInterface>();//create instance

//        diffusionSolver = std::make_shared<fenicsClass>();//create instance
//        diffusionSolver2 = std::make_shared<diffusionPETSc>();//create instance
        std::cout<<"\nDiffusion solver created...\n";

        MPI_Barrier(workers);
//        MPI_Barrier(world);

        //TODO: should populate this perhaps in main, but left here for now
        std::vector<std::string> hslFilePaths;
        for(size_t i(0); i<numHSLGrids; i++)
        {
            auto thisFileString = Params.fileIO->fbase + "H" + std::to_string(i+1) + ".pvd";
            hslFilePaths.push_back(thisFileString);
        }

        auto vecD = std::vector<double>(eQ::parameters["D_HSL"].get<std::vector<double>>());

        fenicsParams.uniqueID  = whichHSLNode;
        fenicsParams.comm      = mpiComms[whichHSLNode];
        fenicsParams.D_HSL     = vecD[whichHSLNode];
        fenicsParams.filePath  = hslFilePaths[whichHSLNode];
        fenicsParams.nodesPerMicron    = double(eQ::parameters["nodesPerMicronSignaling"]);

        //the following are unique to each HSL grid, but common among nodes per HSL grid

        int layerNodes;
        MPI_Comm_size(mpiComms[whichHSLNode], &layerNodes);
        std::cout<<"\ndiffusionSolver->initDiffusion("<<whichHSLNode<<") comm = "<<mpiComms[whichHSLNode]
                <<"  with MPI_Comm_size = "<<layerNodes
                <<std::endl;

        MPI_Barrier(workers);
//        MPI_Barrier(world);

        //diffusion template initialization call
        diffusionSolver->initDiffusion(fenicsParams);
//        diffusionSolver->initDiffusion(whichHSLNode, thisComm, thisPath, thisD,  dt, argc, argv);
//        diffusionSolver2->initDiffusion(workers, hslFilePaths, argc, argv);

        std::cout<<"Diffusion solver "<<whichHSLNode<<" initialized..."<<std::endl;

        //NOTE:  Fenics specific verification and data transfer here
        //TODO:  switch on implementation and/or define common interface for translation of data sent in main loop
        //VERIFY:
        if(diffusionSolver->shell->mesh->num_vertices() != globalNodes)
        {
            std::cout<<"ERROR VERTICES VS. numNodes!"<<std::endl;
        }
//        workerNodeID = diffusionSolver->myRankMPI;


            //  INITIAL MPI DATA TRANSFER
            MPI_Status mpiStatus;

            //SEND the mesh coordinate-to-vertex, and vertex-to-dof maps:
            for(size_t i(0); i<numHSLGrids; i++)
            {
                int thisNode = int(i)+1;//iterate through worker nodes (1,...) but use world comm (0,...)
                if(thisNode == my_PE_num)//each worker node (1,...) sends node data to the controller (node 0)
                {//note:  controller node launches corresponding MPI_Recv to catch each MPI_Send below
                    int mpiDest = 0;  int mpiTag = 0;//destination is node zero, tag is dummy=0

                    //SEND MESH COORDINATE DATA AND VERTEX-TO-DOF MAPPINGS
                    MPI_Send(diffusionSolver->shell->mesh_coords.data(), 2*globalNodes, MPI_DOUBLE,
                             mpiDest, mpiTag, world);
                    MPI_Send(diffusionSolver->shell->dof_from_vertex.data(), globalNodes, MPI_INT,
                             mpiDest, mpiTag, world);
                    std::cout<<"Sent dof/coords data via MPI_Send..."<<std::endl;

                    //RECEIVE THE INITIAL SETTINGS FOR THE DIFFUSION TENSOR
                    MPI_Recv(diffusionSolver->D11->data(), int(globalNodes), MPI_DOUBLE,
                             mpiDest, mpiTag, world, &mpiStatus);
                    MPI_Recv(diffusionSolver->D22->data(), int(globalNodes), MPI_DOUBLE,
                             mpiDest, mpiTag, world, &mpiStatus);
                    MPI_Recv(diffusionSolver->D12->data(), int(globalNodes), MPI_DOUBLE,
                             mpiDest, mpiTag, world, &mpiStatus);
                    std::cout<<"Received diffusion Tensor data via MPI_Recv..."<<std::endl;
                }
                //wait other nodes to sync with receive
                MPI_Barrier(world);
            }

        //BOUNDARY WELL MODEL INITIALIZATION
        //initialization w.r.t. the trap boundary (if used)
        boundaryWellConcentration = 0.0;
        boundaryUnderFlow = 0;
        boundaryDecayRate =  //units: min^-1; = linearFlowRate/trapLength = 100/2000=0.05sec^-1
                (linearFlowRate/double(eQ::parameters["lengthScaling"]))
                * (1.0/double(eQ::parameters["simulationTrapWidthMicrons"]));

        //compute volume of flow channels (+ left/right sides in first approximation)
        wellScaling = 2.0 * 10.0 * (15.0/double(eQ::parameters["lengthScaling"])) //10um z-height, 15um y-height
                * (double(eQ::parameters["simulationTrapWidthMicrons"]) + double(eQ::parameters["simulationTrapHeightMicrons"]));//x channelThickness x channelHeight x #

        if((double(eQ::parameters["dt"]) * boundaryDecayRate) >= 1.0)
            std::cout<<"ERROR!  timestep too large for decay of boundary well concentration: "
                    <<(double(eQ::parameters["dt"]) * boundaryDecayRate)<<std::endl;

    }

//=============================================================================
    //2.  MPI RECEIVE HSL DATA, SEND ANISOTROPIC DIFFUSION TENSOR DATA:
//=============================================================================
    if(isControllerNode)
    {
        //set diffusion tensor to identity matrix (will be updated every time step, so just default values here):
        D11Grid.assign(globalNodes, 1.0);
        D22Grid.assign(globalNodes, 1.0);
        D12Grid.assign(globalNodes, 0.0);

            std::cout<<"\n\tController:  Waiting for MPI RECV...\n\n";


            MPI_Status mpiStatus;
            for(size_t i(0); i<numHSLGrids; i++)
            {
                int mpiSource = int(i)+1;  int mpiTag = 0;

                //RECEIVE MESH COORDINATE DATA AND VERTEX-TO-DOF MAPPINGS
                gridCoords.push_back(std::vector<double>(2*globalNodes, 0.0));//initialized vector to hold result
                MPI_Recv(gridCoords.back().data(), 2*globalNodes, MPI_DOUBLE,
                         mpiSource, mpiTag, world, &mpiStatus);
                gridDofs.push_back(std::vector<int>(globalNodes, 0.0));
                MPI_Recv(gridDofs.back().data(), globalNodes, MPI_INT,
                         mpiSource, mpiTag, world, &mpiStatus);

                //SEND THE INITIAL SETTINGS FOR THE DIFFUSION TENSOR
                MPI_Send(D11Grid.data(), globalNodes, MPI_DOUBLE,
                         mpiSource, mpiTag, world);
                MPI_Send(D22Grid.data(), globalNodes, MPI_DOUBLE,
                         mpiSource, mpiTag, world);
                MPI_Send(D12Grid.data(), globalNodes, MPI_DOUBLE,
                         mpiSource, mpiTag, world);

                MPI_Barrier(world);
            }
            std::cout<<"\n\tController:  Received MPI RECV...\n\n";
    }

    //re-sync after data transfer
    MPI_Barrier(world);

    //store the local nodeIDs of the sub-communicator vs. world:
    //NB: these seem to line-up with the world comm, but just to be sure
//    nodeIds.assign(size_t(npes),0);
//    MPI_Gather(&workerNodeID, 1, MPI_INT,  //send 1 int
//               nodeIds.data(), 1, MPI_INT, //only receiver needs this buffer (1 int from each gathered node)
//               0, world); //gathered in node 0, world comm.

//================================================================
    //3.  BUILD THE HSL LAYER DATA STRUCTURES TO PASS TO THE ABM:
//================================================================
    if(isControllerNode)
    {
        //print the sub-communicator node IDs returned vs. how MPI_COMM_WORLD sees them
        //note: the controller node is set to -1, above
//        std::cout<<"nodeID order: ";
//        for(size_t i(0); i<size_t(npes); i++)
//        {
//            std::cout<<nodeIds[i]<<" ";
//        }
//        std::cout<<std::endl;

        //create vector of HSL grids (numHSLGrids is defined in the Simulation constructor)
        for(size_t grid(0); grid<numHSLGrids; grid++)
        {
            //this will be the solution vector buffer, in dof order (for fenics)
//            HSLGrids.push_back(
//                        std::vector<double>(globalNodes));
            //CREATE THE SOLUTION VECTOR FOR THIS DIFFUSION LAYER:
            //write directly to the ABM parameters data structure:
//            paramsABM.hslSolutionVector.push_back(std::make_shared<std::vector<double>>(globalNodes));
            paramsABM.hslSolutionVector.push_back(std::vector<double>(globalNodes));



            //this will map x,y position directly to dof
//            dof_from_grid.push_back(
//                        std::make_shared<eQ::gridFunction<size_t>>(globalNodesH, globalNodesW));
            //CREATE THE LOOKUP TABLE FOR THIS DIFFUSION LAYER:
            paramsABM.dofLookupTable.push_back(std::make_shared<eQ::gridFunction<size_t>>(globalNodesH, globalNodesW));

            //create a map of vertices to absolute grid position
            //TODO:  see if this is actually used anywhere:
            coords.push_back(
                        std::vector<std::pair<double,double>>(globalNodes));

            //gridCoords is a vector (one entry per diffusion grid) of a vector of x,y pairs that map to the data vector for fenics
            //thus we step through every other entry to get x,y values iteratively (total size 2*globalNodes size)
            for (size_t i(0); i < globalNodes; i++)
            {//mesh_coords is a 1D array of (x,y) pairs in vertex order (not in dof order):
                //convert to a 1D array of std::pair
                auto x = gridCoords[grid][2*i];
                auto y = gridCoords[grid][2*i+1];
                coords.back()[i] = std::pair<double,double>{x,y};//x,y pair indexed by vertex#

                unsigned jx = unsigned(round(x * double(eQ::parameters["nodesPerMicronSignaling"])));
                unsigned iy = unsigned(round(y * double(eQ::parameters["nodesPerMicronSignaling"])));
                //the point of this all is to have this mapping from physical node # to the fenics dof #,
                //allows to look that up the dof directly for read/write in the main acquisition loop
//                dof_from_grid[grid]->grid[iy][jx] = gridDofs[grid][i];//note row,column = y,x

                //WRITE THE DOF LOOKUP TABLE:
                //write directly to the ABM parameters data structure:
                paramsABM.dofLookupTable[grid]->grid[iy][jx] = size_t(gridDofs[grid][i]);
            }
        }
        //to verify the dof are as expected (for small grids to test):
//        printData();
    }

    createdHSLgrid = true;
    MPI_Barrier(world);
}
void Simulation::init_ABM(int numSeedCells)
{
    //ONLY THE CONTROLLER OWNS THE ABM MODEL:
    if(isControllerNode)
    {
        //set the seed value for the ABM random number generator (set to same value for repeatability)
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();

        //set seed value to fixed value if it is set non-zero, else use timestamp seed
        paramsABM.seedValue = (0 == eQ::parameters["fixedSeedValue"])
                ? size_t(seed) : size_t(eQ::parameters["fixedSeedValue"]);

        //COMMON CODE FOR HSL SIGNALING SETUP:
        if(createdHSLgrid)
        {//init rates and strain parameters:
            //note:  these were first the dual-strain osc. parameters.
            //TODO:  define a class to pass rates as needed, or leave this to the Strain class
            init_kinetics_DSO();
            paramsABM.pA = &pA;
            paramsABM.pR = &pR;
        }

        //SWITCH ON SIMULATION TYPE:
/*
        if("DUALSTRAIN_OSCILLATOR" == eQ::parameters["simType"])
        {//need to check grids are init;  check size automatically, etc...
            if(createdHSLgrid)
            {//init rates and strain parameters:
                init_kinetics_DSO();
                params.pA = &pA;
                params.pR = &pR;
                params.c4grid = &HSLGrids[0];
                params.c14grid = &HSLGrids[1];
                //NOTE:  fenics specific;  need to switch on implementation or provide in template
                params.c4lookup = dof_from_grid[0];
                params.c14lookup = dof_from_grid[1];
            }
            else
            {
                eQ::parameters["simType"] = "NO_SIGNALING";
                std::cout<<"HSL grid not init!  Over-write of simType to eQ::simType::NO_SIGNALING"<<std::endl;
            }
        }
        else if(
                ("SENDER_RECEIVER" == eQ::parameters["simType"])
                || ("INDUCED_SENDER_RECEIVER" == eQ::parameters["simType"])
                || ("DUAL_SENDER_RECEIVER" == eQ::parameters["simType"])
                )
        {//need to check grids are init;  check size automatically, etc...
            if(createdHSLgrid)
            {//init rates and strain parameters:
                init_kinetics_DSO();
                params.pA = &pA;
                params.pR = &pR;
                params.c4grid = &HSLGrids[0];
                //NOTE:  fenics specific;  need to switch on implementation or provide in template
                params.c4lookup = dof_from_grid[0];
//                if ("DUAL_SENDER_RECEIVER" == eQ::parameters["simType"])
                {
                    params.c14grid = &HSLGrids[1];
                    params.c14lookup = dof_from_grid[1];
                }
            }
            else
            {
                eQ::parameters["simType"] = "NO_SIGNALING";
                std::cout<<"HSL grid not init!  Over-write of simType to eQ::simType::NO_SIGNALING"<<std::endl;
            }
        }
        else if("MODULUS_1" == eQ::parameters["simType"])
        {//need to check grids are init;  check size automatically, etc...
            if(createdHSLgrid)
            {//init rates and strain parameters:
                init_kinetics_DSO();
                params.pA = &pA;
//                params.pR = &pR;
                params.c4grid = &HSLGrids[0];
//                params.c14grid = &HSLGrids[1];
                //NOTE:  fenics specific;  need to switch on implementation or provide in template
                params.c4lookup = dof_from_grid[0];
//                params.c14lookup = dof_from_grid[1];
            }
        }
        else if("MODULUS_2" == eQ::parameters["simType"])
        {//need to check grids are init;  check size automatically, etc...
            if(createdHSLgrid)
            {//init rates and strain parameters:
                init_kinetics_DSO();
                params.pA = &pA;
                params.pR = &pR;
                params.c4grid = &HSLGrids[0];
                params.c14grid = &HSLGrids[1];
                //NOTE:  fenics specific;  need to switch on implementation or provide in template
                params.c4lookup = dof_from_grid[0];
                params.c14lookup = dof_from_grid[1];
            }
        }

*/

//=============================================================================
    //CREATE ABM CLASS:
//=============================================================================
        paramsABM.dataFiles = &Params.dataFiles;//a list of cell parameters to record using fenics
        ABM = std::make_shared<eQabm>(paramsABM);
        ABM->initCells(numSeedCells);
        simulateABM = true;
    }
    else
        simulateABM = false;
}
void Simulation::stepSimulation()
{
    MPI_Barrier(world);
    if(isControllerNode)
    {
        if(simulateABM)
        {
            double tare = MPI_Wtime();
        //==========================
        // STEP THE PHYSICS ENGINE
        //==========================
            ABM->stepChipmunk();
            physicsTimer += (MPI_Wtime() - tare);
        }

        //get the HSL grid data from the fenics worker nodes:
        if(HSL_signalingTrue)
        {
            mpiRequest.clear();
            mpiRequest.assign(size_t(numHSLGrids), nullptr);
            int mpiTag=0; int sourceNode;
            for(size_t i(0); i<numHSLGrids; i++)
            {
                sourceNode = int(i)+1;//controller is node 0 (destination of Irecv), workers are nodes 1,2,...(matching Isend is below)
//                MPI_Irecv(HSLGrids[i].data(), int(globalNodes), MPI_DOUBLE,
//                          sourceNode, mpiTag, world, &mpiRequest[i]);

//                MPI_Irecv( (*ABM->Params.hslSolutionVector[i].get()).data(), int(globalNodes), MPI_DOUBLE,
                MPI_Irecv(ABM->Params.hslSolutionVector[i].data(), int(globalNodes), MPI_DOUBLE,
                          sourceNode, mpiTag, world, &mpiRequest[i]);
                //matching wait is in waitMPI below
            }
        }
    }
    if(isDiffusionNode)
    {
        if(HSL_signalingTrue)
        {
            double tare = MPI_Wtime();
            //==========================
            // SOLVE HSL DIFFUSION:
            //==========================
            //if Robin system:
//            diffusionSolver->myParams.fenForms->setRobinBoundaryValue(diffusionSolver->boundaryWellConcentration);

            //diffusion template methods:
            diffusionSolver->stepDiffusion();
//            diffusionSolver2->stepDiffusion();

            eQ::parametersType flux    = diffusionSolver->getBoundaryFlux();
//            auto diffusionConstant     = diffusionSolver->getDiffusionConstant();

            //COMPUTE TRAP CHANNEL CONCENTRATION:
            boundaryWellConcentration
//                    += 0.1*double(flux["totalFlux"])/wellScaling;//scale #c from flux integral
                    += double(flux["totalFlux"])/wellScaling;//scale #c from flux integral
            boundaryWellConcentration
                    -= double(eQ::parameters["dt"])*boundaryDecayRate*boundaryWellConcentration;

            if(boundaryWellConcentration < 0.0)
            {
                boundaryWellConcentration = 0.0;
                boundaryUnderFlow++;
            }
            if("DIRICHLET_UPDATE" == eQ::parameters["boundaryType"])
            {
                //note:  set "allBoundaries" to false if not using all boundaries with same value!
                eQ::parametersType bvals;
                bvals["allBoundaries"] = bool(true);
                bvals["boundaryValue"] = boundaryWellConcentration;
                diffusionSolver->setBoundaryValues(bvals);
            }

            diffusionTimer += (MPI_Wtime() - tare);
            //ship the new HSL grids to the controller node:
            mpiRequest.clear();
            mpiRequest.push_back(nullptr);
            int destNode=0;  int mpiTag=0;

            MPI_Isend(diffusionSolver->solution_vector.data(), int(globalNodes), MPI_DOUBLE,
                      destNode, mpiTag, world, &mpiRequest[0]);
            //matching wait is in waitMPI below
        }
    }
    simTime += double(eQ::parameters["dt"]);
//    MPI_Barrier(world);
}

//TODO:  GENERALIZE THIS FOR EACH TYPE OF GENE CIRCUIT
void Simulation::create_DataGrid()
{//for using Fenics for data recording only (populating grid functions)
//    MPI_Barrier(world);
    if(isControllerNode)//only controller makes data grid:
    {
        computeGridParameters();
        //TODO:  need to switch on the trap type for BC
        //create the fenics interface class:
        //**************************************************************
        //              CREATE FENICS CLASS:
        //**************************************************************
        fenicsParams.comm       = controllerComm;//single mpi node for the controller
        fenicsParams.dataFiles = Params.dataFiles;
        fenicsParams.nodesPerMicron = double(eQ::parameters["nodesPerMicronData"]);
        //uses non-template version of diffusion class for fenics
        diffusionSolver = std::make_shared<fenicsInterface>(fenicsParams);
        std::cout<<"Fenics for data recording Class created..."<<std::endl;
    }
//    MPI_Barrier(world);
}

void Simulation::updateCells()
{//UPDATES CELL POSITIONS, DIVIDES CELLS, AND REMOVES CELLS OUTSIDE THE TRAP
    if(isControllerNode && simulateABM)
    {
//        std::cout<<"ABM->updateCellPositions()"<<std::endl;
        ABM->updateCellPositions(simTime);//single threaded; divides and/or removes cells
//        std::cout<<"ABM->updateCellPositions() end"<<std::endl;
    }

    //WAIT FOR MPI DATA TRANSFER TO COMPLETE (HSL DATA TO CONTROLLER NODE)
    if(HSL_signalingTrue)
    {
        waitMPI();//wait for HSL data transfer to complete (all nodes)
    }

    if(isControllerNode)
    {
        //SETUP THE CELL LIST AND DIVIDE TO THREADS:
        double tare = MPI_Wtime();

//        std::cout<<"updateCellModels"<<std::endl;
        if(simulateABM)
            ABM->updateCellModels(size_t(npes));//multi-threaded, passes world pe count (total)
//        std::cout<<"updateCellModels end"<<std::endl;

        computeTimer += (MPI_Wtime() - tare);
    }

    //non-controller nodes will wait here (but their cpu cores will be used in shared memory above)
    MPI_Barrier(world);

    //SEND THE MODIFIED HSL DATA BACK TO THE HSL WORKER NODES:
    if(HSL_signalingTrue)
    {
        writebackGridData();
    }
}
void Simulation::writebackGridData()
{//this writes-back data from the controller to the hsl grids:
    MPI_Barrier(world);
    if(isControllerNode)
    {
        //ship the modified HSL grid data to the fenics worker nodes:
        if(HSL_signalingTrue)
        {
            mpiRequest.clear();
            mpiRequest.assign(size_t(numHSLGrids), nullptr);
            int mpiTag=0; int destNode;

            mpiRequestA.assign(size_t(numHSLGrids), nullptr);
            mpiRequestB.assign(size_t(numHSLGrids), nullptr);
            mpiRequestC.assign(size_t(numHSLGrids), nullptr);

            eQ::linearize2Dgrid<double>(ABM->D11grid, D11Grid);
            eQ::linearize2Dgrid<double>(ABM->D22grid, D22Grid);
            eQ::linearize2Dgrid<double>(ABM->D12grid, D12Grid);

            for(size_t i(0); i<numHSLGrids; i++)
            {
                destNode = int(i)+1;
//                MPI_Isend(HSLGrids[i].data(), int(globalNodes), MPI_DOUBLE,
//                          destNode, mpiTag, world, &mpiRequest[i]);
//                MPI_Isend((*ABM->Params.hslSolutionVector[i].get()).data(), int(globalNodes), MPI_DOUBLE,
                MPI_Isend(ABM->Params.hslSolutionVector[i].data(), int(globalNodes), MPI_DOUBLE,
                          destNode, mpiTag, world, &mpiRequest[i]);

                //JW:  petsc Isend here (use unique mpiRequest vector)...
                //...
                //matching wait is in waitMPI below

                if(true)//todo: switch on a flag to send
                {
                    MPI_Isend(D11Grid.data() , int(globalNodes), MPI_DOUBLE,
                              destNode, mpiTag, world, &mpiRequestA[i]);
                    MPI_Isend(D22Grid.data(), int(globalNodes), MPI_DOUBLE,
                              destNode, mpiTag, world, &mpiRequestB[i]);
                    MPI_Isend(D12Grid.data(), int(globalNodes), MPI_DOUBLE,
                              destNode, mpiTag, world, &mpiRequestC[i]);
                }
            }
        }
    }
    if(isDiffusionNode)
    {
        if(HSL_signalingTrue)
        {
            //get the new HSL grids from the controller node:
            mpiRequest.clear();
            mpiRequest.push_back(nullptr);
            int sourceNode=0;  int mpiTag=0;
            MPI_Irecv(diffusionSolver->solution_vector.data(), int(globalNodes), MPI_DOUBLE,
                      sourceNode, mpiTag, world, &mpiRequest[0]);
            //JW:  petsc Irecv here (use unique mpiRequest vector)...
            //...
            //matching wait is in waitMPI below

            if(true)//todo: switch on a flag to send
            {
                mpiRequestA.assign(size_t(numHSLGrids), nullptr);
                mpiRequestB.assign(size_t(numHSLGrids), nullptr);
                mpiRequestC.assign(size_t(numHSLGrids), nullptr);

                MPI_Irecv(diffusionSolver->D11->data(), int(globalNodes), MPI_DOUBLE,
                          sourceNode, mpiTag, world, &mpiRequestA[0]);
                MPI_Irecv(diffusionSolver->D22->data(), int(globalNodes), MPI_DOUBLE,
                          sourceNode, mpiTag, world, &mpiRequestB[0]);
                MPI_Irecv(diffusionSolver->D12->data(), int(globalNodes), MPI_DOUBLE,
                          sourceNode, mpiTag, world, &mpiRequestC[0]);
            }
        }
    }
}
void Simulation::stepFinalize()
{
    if(HSL_signalingTrue)
    {
        waitMPI();//wait for HSL data transfer to complete (all nodes)\n"
    }
}
void Simulation::waitMPI()
{
    if(HSL_signalingTrue)
    {
       if(isControllerNode)
       {
           double tare = MPI_Wtime();
           for(size_t i(0); i<numHSLGrids; i++)
           {
               MPI_Status thisStaus;
               MPI_Wait(&mpiRequest[i], &thisStaus);
               //JW: petsc wait here
//               MPI_Wait(&mpiRequestPETSC[i], &thisStaus);

               if(!mpiRequestA.empty())
               {
                   MPI_Wait(&mpiRequestA[i], &thisStaus);
                   MPI_Wait(&mpiRequestB[i], &thisStaus);
                   MPI_Wait(&mpiRequestC[i], &thisStaus);
               }
           }
           waitTimer += (MPI_Wtime() - tare);           
       }
       if(isDiffusionNode)
       {
           double tare = MPI_Wtime();
           MPI_Status thisStaus;
           MPI_Wait(&mpiRequest[0], &thisStaus);

           //JW: petsc wait here
//           MPI_Wait(&mpiRequestPETSC[0], &thisStaus);

           if(!mpiRequestA.empty())
           {
               MPI_Wait(&mpiRequestA[0], &thisStaus);
               MPI_Wait(&mpiRequestB[0], &thisStaus);
               MPI_Wait(&mpiRequestC[0], &thisStaus);
           }
           waitTimer += (MPI_Wtime() - tare);
       }
       mpiRequestA.clear();
       mpiRequestB.clear();
    }
}
void Simulation::writeHSLFiles()
{
    if((isDiffusionNode) && createdHSLgrid)
    {
//        std::cout<<"writing hsl files..."<<std::endl;

        diffusionSolver->writeDiffusionFiles(simTime);
//        diffusionSolver2->writeDiffusionFiles(simTime);

        //print the l2 norm and boundary flux of the HSL:
//        auto thisNorm = diffusionSolver->u->vector()->norm("l2");
        auto flux = diffusionSolver->getBoundaryFlux();
        std::cout<<"boundaryFlux assemble result: ("<<my_PE_num<<") "<<flux["totalFlux"]
//                    <<"  u norm: "<<thisNorm
                   <<" well conc.: "<<boundaryWellConcentration
//                  <<" boundary underflow counter: "<<boundaryUnderFlow
//                 <<" boundary finite counter: "<<boundaryFinite
                  <<std::endl;

    }
}
void Simulation::writeDataFiles()
{//when using fenics to record cell data over the grid (not for HSL data here)
    if(isControllerNode)
    {
        //compute average value:
//        for(size_t i(0); i<nodesHighData; i++)
//        {
//            for(size_t j(0); j<nodesWideData; j++)
//            {
//                if(ABM->gridDataCounter->grid[i][j] != 0)
//                {
//                    Params.dataFiles[0].first->dataGrid->grid[i][j]
//                            /= double(ABM->gridDataCounter->grid[i][j]);
////                    //vector data:
////                    Params.dataFiles[1].first->xdataGrid->grid[i][j]
////                            /= double(ABM->gridDataCounter->grid[i][j]);
////                    Params.dataFiles[1].first->ydataGrid->grid[i][j]
////                            /= double(ABM->gridDataCounter->grid[i][j]);
//                }
//            }
//        }

        //temp...comment out here:
        //don't divide FP snapshots...
        diffusionSolver->writeDataFiles(simTime);


        //revert:
//        for(size_t i(0); i<nodesHighData; i++)
//        {
//            for(size_t j(0); j<nodesWideData; j++)
//            {
//                if(ABM->gridDataCounter->grid[i][j] != 0)
//                {
//                    Params.dataFiles[0].first->dataGrid->grid[i][j]
//                            *= double(ABM->gridDataCounter->grid[i][j]);
////                    //vector data:
////                    Params.dataFiles[1].first->xdataGrid->grid[i][j]
////                            *= double(ABM->gridDataCounter->grid[i][j]);
////                    Params.dataFiles[1].first->ydataGrid->grid[i][j]
////                            *= double(ABM->gridDataCounter->grid[i][j]);
//                }
//            }
//        }
    }
}
void Simulation::finalizeDataRecording(void)
{//end-of-simulation final data dump
    if(isControllerNode && simulateABM)
    {
        ABM->finalizeDataRecording(Params.fileIO->fbase);
    }
}
Simulation::~Simulation()
{
    std::cout<<"BEGIN Simulation::~Simulation()"<<std::endl;

//**************************************************************************//
//                  SHUTDOWN MPI SEQUENCE:
//**************************************************************************//
    MPI_Barrier(world);
    if( (bool(eQ::parameters["hslSignaling"])) && !isControllerNode)
    {
        diffusionSolver->finalize();
//        diffusionSolver2->finalize();
    }
//    diffusionSolver.reset();
//    diffusionSolver2.reset();

//    auto p = mpiComms.data();
//    for(size_t i(0); i<mpiComms.size(); ++i)
//    {
//        if(whichHSLNode == i)
//            if(MPI_COMM_NULL != p[i])
//                MPI_Comm_free(&p[i]);
//    }


    if(world != workers)
        if(MPI_COMM_NULL != workers)
            MPI_Comm_free(&workers);
    if(world != controllerComm)
        if(MPI_COMM_NULL != controllerComm)
            MPI_Comm_free(&controllerComm);

    MPI_Barrier(world);
    dolfin::SubSystemsManager::singleton().finalize();
    std::cout<<"END Simulation::~Simulation()"<<std::endl;
}
void Simulation::printData()
{
    size_t count = 0;
    for(size_t rank=0; rank<numHSLGrids; rank++)
    {
        std::cout<<"GRID #"<<rank<<std::endl;
        for(int ii=int(globalNodesH)-1; ii>=0; ii--)//descend for pretty-printing
        {
            size_t i = size_t(ii);
            for(size_t j=0; j<globalNodesW; j++)
            {
//                std::cout<< std::setw(3)<<std::left
//                    <<dof_from_grid[rank]->grid[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
    }
}
void Simulation::init_kinetics_DSO()
{
    //writes the master table (defined in *strain.cpp file):
    initRates(basalRate, rates);
    //loads rates chosen for this simulation; pA,pR are the global rate tables:
    loadParams(eQ::strainType::ACTIVATOR,pA,rates);
    loadParams(eQ::strainType::REPRESSOR,pR,rates);

    if("DUALSTRAIN_OSCILLATOR" == eQ::parameters["simType"])
    {
        //make rhlI production in sender more sensitive to autoinduce:
        pA.K.H /= 20.0;//result ~=300 nM
        pR.K.H /= 20.0;//for promoter in receiver
        pA.K.L /= 10.0;
//        pR.K.L /= 10.0;//for promoter in receiver
        pA.K.I /= 40.0;
//        pR.K.I /= 10.0;//for promoter in receiver
    }
    else
    {
//        //make rhlI production in sender more sensitive to autoinduce:
//        pA.K.H /= 20.0;//result ~=300 nM
//        pR.K.H /= 20.0;//for promoter in receiver
    }

    //scale HSL production rate (due to diffusion boundaries)
//    double rateScale = 1.0e4;

    //zero rate scale:
//    rateScale = 1.0;
//    pA.phi    *= rateScale;
//    pR.phi    *= rateScale;
//    pA.K.H    *= rateScale;   //set to decrease response to  C4 in Activator (n=4)
//    pR.K.I    *= rateScale;    //sets higher thresh for repressor C14 repression (n=4)

    //when using ABM model for DSO:
//    pA.K.I    /= 100.0;   //set to increase response to low C14 in Activator (n=4)
//    pR.K.H    /= 100.0;   //set to increase response to low C4 in repressor (n=4)

    //test cutting C14 membrane diffusion rate:
//    pA.DI /= 10.0;
//    pR.DI /= 10.0;
}
//TODO:  need to make this pass the strain type, or switch on it to use for more than just DSO
void Simulation::create_DSOgrid()
{
    //only the controller node should call this:
    //note:  without ABM, the controller node will also be an hsl node.
    computeGridParameters();
    init_kinetics_DSO();

    std::cout<<"Creating dsoGrid of size H, W (dt):"<<globalNodesH<<", "<<globalNodesW
            <<" ("<<eQ::parameters["dt"]<<")"<<std::endl;

    //CREATE 2D GRID OF STRAIN POINTERS AND INSTANTIATE WITH STRAIN TYPE:
    //for DSO, use two grids (one activator, one repressor at each lattice point)
    dsoGrid.push_back(std::make_shared<eQ::gridFunction<std::shared_ptr<Strain>>>(globalNodesH, globalNodesW));
    dsoGrid.push_back(std::make_shared<eQ::gridFunction<std::shared_ptr<Strain>>>(globalNodesH, globalNodesW));
    for(size_t i(0); i<globalNodesH; i++)
    {
        for(size_t j(0); j<globalNodesW; j++)
        {
            dsoGrid[0]->grid[i][j] =
                        std::make_shared<Strain>(eQ::strainType::ACTIVATOR, &pA,
                                                 eQ::parameters["dt"], eQ::parameters["nodesPerMicronSignaling"]);
            dsoGrid[1]->grid[i][j] =
                        std::make_shared<Strain>(eQ::strainType::REPRESSOR, &pR,
                                                 eQ::parameters["dt"], eQ::parameters["nodesPerMicronSignaling"]);
        }
    }

    simulateABM = false;
}
