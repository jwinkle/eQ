#include "simulation.h"


Simulation::Simulation(MPI_Comm commWorld, const Simulation::Params &initParams)
    : params(initParams), world(commWorld)
{
    //Simulation constructor is called by all nodes;  setup MPI partitions first:
    //==========================================================================//
    //                  MPI INITIALIZATION:
    //==========================================================================//
    MPI_Comm_rank(world, &my_PE_num);
    MPI_Comm_size(world, &npes);
    MPI_Barrier(world);

    // 1.  CREATE MPI COMMUNICATORS FOR CONTROLLER AND "WORKERS" (DIFFUSION NODES)
    // 2.  INITIALIZE TIMERS, HSL GRID#, #NODES FOR DATA RECORDING
    // 3.  CREATE THE ABM CLASS

    //NOTE:  their are two input parameters: #HSL layers and #MPI nodes
    //  --determine their relationships for well-formed signaling via MPI:
    //  Possible states:
    //  1 MPI rank -- no MPI transfer, HSL signaling not allowed
    //  >1 MPI ranks -- HSL signaling allowed only if HSL vector size divides "worker" #MPI nodes

    numHSLGrids = size_t(eQ::data::parameters["D_HSL"].size());//size of diffusion vector is # of grids (known globally)

    if(1==npes)
    {
        std::cout<<"Simulation class initializing with 1 pe."<<std::endl;
        HSL_signaling   = false;
        workers         = world;
        controllerComm  = world;
        isControllerNode = true;
        //verify this is correct to do:
        isDiffusionNode = false;
        whichHSL=0;
    }
    else if( (numHSLGrids > 0) && (size_t(npes-1) % numHSLGrids == 0) )
    {//each layer gets same number of mpi processing elements:
        mpiNodesPerDiffusionLayer = size_t(npes-1) / numHSLGrids;
        HSL_signaling = true;
    }
    else
    {
        std::cout<<"Error: Diffusion list size does not divide number of MPI worker nodes!"<<std::endl;
        HSL_signaling = false;
        mpiNodesPerDiffusionLayer = 0;
    }

    if(isControllerNode)
    {
    //OUTPUT INFORMATIONAL DATA (TODO: use a log class to pipe information)
            std::cout<<std::endl;
            std::cout << "Data Grid Recording: "        <<params.dataFiles->size()<<" grids."<<std::endl;
            std::cout << "MPI_COMM_WORLD has "          <<npes<<" nodes"<< std::endl;
            std::cout << "FENICS VERSION:"              << std::endl;
                std::cout << "\t dolfin_version(): "    << dolfin::dolfin_version() << std::endl;
                std::cout << "\t ufc_signature(): "     << dolfin::ufc_signature() << std::endl;
                std::cout << "\t git_commit_hash(): "   << dolfin::git_commit_hash() << std::endl;
            std::cout << "argc= "                       <<params.argc<<std::endl;
            std::cout<<"argv[] = ";
            for (int i=0; i<params.argc; i++)
                std::cout<<params.argv[i]<<" ";
            std::cout<<std::endl;

            char procName[1024];
            int resultLength;
            MPI_Get_processor_name(procName, &resultLength);
            std::cout<<"MPI_Get_processor_name: "       <<procName<<std::endl;

            double tick = MPI_Wtick();
            std::cout<<"MPI_Wtime has MPI_Wtick: "      <<tick<<" seconds"<<std::endl;
            std::cout<<std::endl;
            std::cout<<std::endl;
    }


    if(HSL_signaling)
    {//SAME # PEs (>=1) FOR EACH DIFFUSION LAYER (ALWAYS ONE FOR CONTROLLER LAYER):
        mpiAssignCommunicators();
        //deterimine who is controller, who is worker and verify # workers = number of grids requested
        isControllerNode = (0 == my_PE_num);
        isDiffusionNode = !isControllerNode;
        if(isControllerNode)
        {
            int n=0;
            std::cout<<"CONTROLLER:"<<std::endl;
            std::cout<<"\t World: \t"<<world<<std::endl;
            std::cout<<"\t Controller: \t"<<controllerComm<<std::endl;
            std::cout<<"\t Workers: \t"<<workers<<std::endl;
            for(auto &node:mpiRanks)
            {
                std::cout<<"HSL Node ("<<n<<") \t"<<mpiComms[n]<<std::endl;
                ++n;
                for(auto &rank:node)
                    std::cout<<"MPI_COMM_WORLD rank : \t"<<rank<<std::endl;
            }
        }
        MPI_Barrier(world);
        if(isDiffusionNode)
        {
            size_t n=0;
            for(auto &node:mpiRanks)
            {
                if(n == whichHSL)
                {
                    std::cout<<"DIFFUSION LAYER:"<<std::endl;
                    std::cout<<"\t World: \t"<<world<<std::endl;
                    std::cout<<"\t Controller: \t"<<controllerComm<<std::endl;
                    std::cout<<"\t Workers: \t"<<workers<<std::endl;
                    std::cout<<"HSL Node ("<<n<<") \t"<<mpiComms[n]<<std::endl;
                    for(auto &rank:node)
                        std::cout<<"MPI_COMM_WORLD rank : \t"<<rank<<std::endl;
                }
                ++n;
                MPI_Barrier(workers);
                sleep(1);//sleep to synch output
            }
            std::cout<<std::endl;
        }
        MPI_Barrier(world);
    }
    //result: MPI nodes partitioned into controllerComm (root node), workers (diffusion solver nodes)

    resetTimers();

    MPI_Barrier(world);

    //compute number of nodes for data recording (uses different resolution of nodes/micron than for HSL signaling):
    nodesHighData = size_t(eQ::data::parameters["simulationTrapHeightMicrons"])*size_t(eQ::data::parameters["nodesPerMicronData"])+ 1;
    nodesWideData = size_t(eQ::data::parameters["simulationTrapWidthMicrons"])*size_t(eQ::data::parameters["nodesPerMicronData"])+ 1;

    if(isDiffusionNode)
    {
        mpiController.init(world, 0);
    }
    if(isControllerNode)
    {
        mpiHSL.assign(numHSLGrids, mpiController);//init vector by default copy constructing with (uninitialized) mpiController member
        for(size_t i(0); i<numHSLGrids; ++i)
        {
            int thisNode = int(i)+1;//iterate through worker nodes (1,...) but use world comm (0,...)
            mpiHSL[i].init(world, thisNode);
        }

        //=============================================================================
            //CREATE ABM CLASS:
        //=============================================================================
        //set the seed value for the ABM random number generator (set to same value for repeatability)

        paramsABM.zeroOne = params.zeroOne;
        //note:  these were first the dual-strain osc. parameters.
        //TODO:  define a class to pass rates as needed, or leave this to the Strain class
            init_kinetics_DSO();

        paramsABM.dataFiles = params.dataFiles;//a list of cell parameters to record using fenics
        ABM = std::make_shared<eQabm>(paramsABM);

    }
    MPI_Barrier(world);
}

void Simulation::computeGridParameters()
{
    fenicsParams.filePath           = params.fileIO->fbase;
    fenicsParams.dt                = double(eQ::data::parameters["dt"]);
    fenicsParams.trapWidthMicrons  = double(eQ::data::parameters["simulationTrapWidthMicrons"]);
    fenicsParams.trapHeightMicrons = double(eQ::data::parameters["simulationTrapHeightMicrons"]);

    globalNodesH = size_t(eQ::data::parameters["simulationTrapHeightMicrons"])*size_t(eQ::data::parameters["nodesPerMicronSignaling"]) + 1;
    globalNodesW = size_t(eQ::data::parameters["simulationTrapWidthMicrons"])*size_t(eQ::data::parameters["nodesPerMicronSignaling"]) + 1;
    globalNodes = globalNodesH * globalNodesW;


    auto uleft = double(eQ::data::parameters["simulationChannelLengthLeft"]);
    auto uright = double(eQ::data::parameters["simulationChannelLengthRight"]);
    nodesForChannels = unsigned(
                ceil( (double(eQ::data::parameters["simulationTrapWidthMicrons"]) + uleft + uright) * double(eQ::data::parameters["nodesPerMicronSignaling"]) ));
    nodesForChannels++;//always one more vertex for each dimension vs. # elements
}
void Simulation::create_HSLgrid()
{//all nodes (controller and workers) call this function

    MPI_Barrier(world);

    //DETERMINE NUMBER OF GRID POINTS:
    computeGridParameters();

    if(isControllerNode)
    {
        std::cout<<"Creating trap (hxw): "
                <<eQ::data::parameters["simulationTrapHeightMicrons"]<<" x "
               <<eQ::data::parameters["simulationTrapWidthMicrons"]
              <<" um^2 ==> "
             <<globalNodes<<" global grid nodes"
            <<std::endl;
    }

//=============================================================================
    //1.  DIFFUSION NODES CREATE DIFFUSION CLASS AND TRANSFER MAPPINGS
//=============================================================================
    if(isDiffusionNode)
    {//only workers compute fenics:
    //**************************************************************
    //              CREATE FENICS CLASS:
    //**************************************************************
        diffusionSolver = std::make_shared<fenicsInterface>();//create instance

//        if(bool(eQ::data::parameters["PETSC_SIMULATION"]))
//        {
//            diffusionSolver2 = std::make_shared<diffusionPETSc>();//create instance
//            diffusionSolver2->topBoundaryValue.assign(globalNodesW, 0.0);
//            diffusionSolver2->bottomBoundaryValue.assign(globalNodesW, 0.0);
//        }
//        std::cout<<"\nDiffusion solver created...\n";


        //TODO: should populate this perhaps in main, but left here for now
        std::vector<std::string> hslFilePaths;
        std::vector<std::string> topChannelFilePaths;
        std::vector<std::string> bottomChannelFilePaths;
        for(size_t i(0); i<numHSLGrids; i++)
        {
            auto thisFileString = params.fileIO->fbase + "H" + std::to_string(i+1) + ".pvd";
            hslFilePaths.push_back(thisFileString);

            thisFileString = params.fileIO->fbase + "/channels/ctH" + std::to_string(i+1) + ".pvd";
            topChannelFilePaths.push_back(thisFileString);
            thisFileString = params.fileIO->fbase + "/channels/cbH" + std::to_string(i+1) + ".pvd";
            bottomChannelFilePaths.push_back(thisFileString);
        }

        auto vecD = std::vector<double>(eQ::data::parameters["D_HSL"].get<std::vector<double>>());

        fenicsParams.uniqueID  = whichHSL;
        fenicsParams.comm      = mpiComms[whichHSL];
        fenicsParams.D_HSL     = vecD[whichHSL];
        fenicsParams.filePath  = hslFilePaths[whichHSL];
        fenicsParams.nodesPerMicron    = double(eQ::data::parameters["nodesPerMicronSignaling"]);
        fenicsParams.filePathTopChannel  = topChannelFilePaths[whichHSL];
        fenicsParams.filePathBottomChannel  = bottomChannelFilePaths[whichHSL];


        diffusionSolver->initDiffusion(fenicsParams);

/*
        if(bool(eQ::data::parameters["PETSC_SIMULATION"]))
        {
            fenicsParams.argc = argc;
            fenicsParams.argv = argv;
            fenicsParams.filePath = params.fileIO->fbase;

            diffusionSolver2->initDiffusion(fenicsParams);

            std::vector<double> thisData;

            if("MICROFLUIDIC_TRAP" == eQ::data::parameters["boundaryType"])
            {//USE CHANNEL SOLUTION FOR TOP/BOTTOM:
                //top channels are dirichlet (set N to 0)
                diffusionSolver2->gridData->topNeumannCoefficient = 0.0;
                diffusionSolver2->gridData->bottomNeumannCoefficient = 0.0;
                diffusionSolver2->gridData->topDirichletCoefficient = 1.0;
                diffusionSolver2->gridData->bottomDirichletCoefficient = 1.0;
                diffusionSolver2->topBoundaryValue = diffusionSolver->topChannelData;
                diffusionSolver2->bottomBoundaryValue = diffusionSolver->bottomChannelData;
            }
            else
            {//USE FIXED BOUNDARY CONDITION FOR TOP/BOTTOM
                thisData = std::vector<double>(eQ::data::parameters["boundaries"]["top"][1].get<std::vector<double>>());
                    diffusionSolver2->gridData->topNeumannCoefficient = thisData[0];
                    diffusionSolver2->gridData->topDirichletCoefficient = thisData[1];
                    diffusionSolver2->topBoundaryValue.assign(globalNodesW, thisData[2]);
                thisData = std::vector<double>(eQ::data::parameters["boundaries"]["bottom"][1].get<std::vector<double>>());
                    diffusionSolver2->gridData->bottomNeumannCoefficient = thisData[0];
                    diffusionSolver2->gridData->bottomDirichletCoefficient = thisData[1];
                    diffusionSolver2->bottomBoundaryValue.assign(globalNodesW, thisData[2]);
            }

            //LEFT/RIGHT BOUNDARY CONDITIONS:
            thisData = eQ::data::parameters["boundaries"]["left"][1].get<std::vector<double>>();
                    diffusionSolver2->gridData->leftNeumannCoefficient = thisData[0];
                    diffusionSolver2->gridData->leftDirichletCoefficient = thisData[1];
                    diffusionSolver2->gridData->leftBoundaryValue = thisData[2];
            thisData = eQ::data::parameters["boundaries"]["right"][1].get<std::vector<double>>();
                    diffusionSolver2->gridData->rightNeumannCoefficient = thisData[0];
                    diffusionSolver2->gridData->rightDirichletCoefficient = thisData[1];
                    diffusionSolver2->gridData->rightBoundaryValue = thisData[2];

            //set petsc boundary conditions here:
            diffusionSolver2->gridData->topDirichletCoefficient = 0.0;


            std::cout<<"PETSc solver, path: "<<fenicsParams.filePath<<" initialized..."<<std::endl;
        }
*/

        //VERIFY:
        if(diffusionSolver->shell->mesh->num_vertices() != globalNodes)
        {
            std::cout<<"ERROR VERTICES VS. numNodes!"<<std::endl;
        }

        //  INITIAL MPI DATA TRANSFER
        MPI_Barrier(world);
        std::cout<<"Diffusion node "<<whichHSL<<":  Sending dof/coords data via MPI_Send..."<<std::endl;
        mpiController << eQ::mpi::method::SEND;
        mpiController << diffusionSolver->shell->mesh_coords
                      << diffusionSolver->shell->dof_from_vertex;

        mpiController >> eQ::mpi::method::RECV;
        mpiController >> diffusionSolver->D11
                        >> diffusionSolver->D22
                        >> diffusionSolver->D12;
        std::cout<<"Diffusion node "<<whichHSL<<":  Received diffusion Tensor data via MPI_Recv..."<<std::endl;
        MPI_Barrier(world);

        //homogeneous boundary implementation (superseded by channel flow implementation, but maybe used)
        initBoundaryWell();

    }

//=============================================================================
    //2.  MPI RECEIVE HSL DATA, SEND ANISOTROPIC DIFFUSION TENSOR DATA:
//=============================================================================
    if(isControllerNode)
    {
        std::cout<<"\n\tController:  Transfering initial data from/to HSL nodes via MPI RECV/SEND..."<<std::endl<<std::endl;

        //set diffusion tensor to identity matrix (will be updated every time step, so just default values here):
        //init the grid to isotropic diffusion: (these grids copy directly to the .ufl file tensor)
        //note:  these are 1D vectors using standard translation for (i,j) entries using the eQ helper functions
        D11Grid.assign(globalNodes, 1.0);//diagonal scaling
        D22Grid.assign(globalNodes, 1.0);//diagonal scaling
        D12Grid.assign(globalNodes, 0.0);//symmetric, off-diagonal scaling


        //added an hsl array to record to json file:
        hslVector.assign(numHSLGrids, std::vector<double>(globalNodes, 0.0));
        hslLookup.assign(numHSLGrids, std::vector<eQ::nodeType>(globalNodes, 0));

        MPI_Barrier(world);
        for(auto &node : mpiHSL)
        {
            gridCoords.push_back(std::vector<double>(2*globalNodes, 0.0));//initialized vector to hold result
            gridDofs.push_back(std::vector<int>(globalNodes, 0.0));

            node >> eQ::mpi::method::RECV;
            node >> gridCoords.back()
                    >> gridDofs.back();
            std::cout<<"\n\tController:  Received dof/coords data via MPI_Recv..."<<std::endl;
        }
        for(auto &node : mpiHSL)
        {
            std::cout<<"\n\tController:  Sending diffusion tensor data via MPI_Send to node: "<<node.index()<<std::endl;
            node << eQ::mpi::method::SEND;
            node << D11Grid
                    <<D22Grid
                    <<D12Grid;
        }
        MPI_Barrier(world);
    }

//================================================================
    //3.  BUILD THE HSL LAYER DATA STRUCTURES TO PASS TO THE ABM:
//================================================================
    if(isControllerNode)
    {
        double npm = double(eQ::data::parameters["nodesPerMicronSignaling"]);

        for(size_t grid(0); grid<numHSLGrids; grid++)
        {
            //this will be the solution vector buffer, in dof order (for fenics)
            ABM->createDataVectors(globalNodesH, globalNodesW);

            //gridCoords is a vector (one entry per diffusion grid) of a vector of x,y pairs that map to the data vector for fenics
            //thus we step through every other entry to get x,y values iteratively (total size 2*globalNodes size)
            for (size_t i(0); i < globalNodes; ++i)
            {//mesh_coords is a 1D array of (x,y) pairs in vertex order (not in dof order):
                auto x = gridCoords[grid][2*i];
                auto y = gridCoords[grid][2*i+1];
                auto gridPoint = eQ::data::ij_from_xy(x,y,npm);
                //the point of this all is to have this mapping from physical node # to the fenics dof #,
                //allows to look that up the dof directly for read/write in the main acquisition loop
                //WRITE THE DOF LOOKUP TABLE:
                ABM->writeLookupTable(grid, gridPoint, eQ::nodeType(gridDofs[grid][i]));
            }

            ABM->dofLookupTable[grid]->linearize2Dgrid(hslLookup[grid]);
        }
        //to verify the dof are as expected (for small grids to test):
//        printData();
        std::cout<<"\n\tController:  createdHSLgrid = true...\n\n";
    }

    createdHSLgrid = true;
    MPI_Barrier(world);
}
void Simulation::init_ABM(int numSeedCells, std::vector<std::shared_ptr<Strain>> &strains)
{
    std::vector<std::shared_ptr<Ecoli>> cells;
    //ONLY THE CONTROLLER OWNS THE ABM MODEL:
    if(isControllerNode)
    {
        if("RANDOM" == eQ::data::parameters["cellInitType"])
        {
            ABM->initCells(eQabm::initType::RANDOM, numSeedCells, strains);
        }
        else
        {
            ABM->initCells(eQabm::initType::BANDED, numSeedCells, strains);
        }
        simulateABM = true;
    }
    else
        simulateABM = false;
}
void Simulation::stepSimulation(double simTime)
{
    MPI_Barrier(world);
    if(isControllerNode)
    {
        double tare = MPI_Wtime();
        // STEP THE PHYSICS ENGINE        
            ABM->stepSimulation();

        physicsTimer += (MPI_Wtime() - tare);

        if(HSL_signaling)
        {
            //get the HSL grid data from the fenics worker nodes:
            for(auto &node : mpiHSL)
            {
                node >> eQ::mpi::method::IRECV;
                node >> ABM->hslSolutionVector[node.index()];
            }
        }

        //UPDATES CELL POSITIONS, DIVIDES CELLS, AND REMOVES CELLS OUTSIDE THE TRAP
        ABM->updateCellData(simTime);//single threaded; divides and/or removes cells

        //WAIT FOR MPI DATA TRANSFER TO COMPLETE (HSL DATA TO CONTROLLER NODE)
        if(HSL_signaling)
        {
            waitMPI();//wait for HSL data transfer to complete (all nodes)
            //buffer the original HSL data
            for(auto &node : mpiHSL)
            {
                *ABM->hslSolutionBuffer[node.index()] = *ABM->hslSolutionVector[node.index()];//vector deep copy
            }

        }
        tare = MPI_Wtime();

        //SETUP THE CELL LIST AND DIVIDE TO THREADS:
//        ABM->updateCellModels(size_t(npes));//multi-threaded, passes world pe count (total)
        ABM->updateCellModels();//single-threaded

        //convert to a 1D array to ship via MPI (TODO: tie this to the grid function directly)
        ABM->D11grid->linearize2Dgrid(D11Grid);
        ABM->D22grid->linearize2Dgrid(D22Grid);
        ABM->D12grid->linearize2Dgrid(D12Grid);

        computeTimer += (MPI_Wtime() - tare);
        //non-controller nodes will wait here (but their cpu cores will be used in shared memory above)
        MPI_Barrier(world);

        if(HSL_signaling)
        {
            for(auto &node : mpiHSL)
            {
                node << eQ::mpi::method::ISEND;
                node << ABM->hslSolutionVector[node.index()];

                node << D11Grid
                        << D22Grid
                        << D12Grid;
            }
        }
    }
    if(isDiffusionNode)
    {
        if(HSL_signaling)
        {
            //==========================
            // SOLVE HSL DIFFUSION:
            //==========================
            //diffusion template methods:
            double tare = MPI_Wtime();
            diffusionSolver->stepDiffusion();
            diffusionTimer += (MPI_Wtime() - tare);

//            if(bool(eQ::data::parameters["PETSC_SIMULATION"]))
//            {
//                if("MICROFLUIDIC_TRAP" == eQ::data::parameters["boundaryType"])
//                {
//                    diffusionSolver2->topBoundaryValue = diffusionSolver->topChannelData;
//                    diffusionSolver2->bottomBoundaryValue = diffusionSolver->bottomChannelData;
//                }

//                double ptare = MPI_Wtime();
//                diffusionSolver2->stepDiffusion();
//                petscTimer += (MPI_Wtime() - ptare);
//            }

            mpiController << eQ::mpi::method::ISEND;
            mpiController << diffusionSolver->solution_vector;

            computeBoundaryWell();

            waitMPI();//wait for HSL data transfer to complete (all nodes)
            //non-controller nodes will wait here (but their cpu cores will be used in shared memory above)
            MPI_Barrier(world);

            mpiController >> eQ::mpi::method::IRECV;
            mpiController >> diffusionSolver->solution_vector;

            mpiController >> diffusionSolver->D11
                            >> diffusionSolver->D22
                            >> diffusionSolver->D12;
        }
    }
    simTime += double(eQ::data::parameters["dt"]);
    //NOTE:  DATA XFER TO HSL WORKER NODES IS STILL OPEN HERE;  ALL NODES CONTINUE WITHOUT A BARRIER...
}
void Simulation::stepFinalize()
{
    if(HSL_signaling)
    {
        waitMPI();//wait for HSL data transfer to complete (all nodes)\n"
    }
}
void Simulation::waitMPI()
{
//    std::cout<<"Entering waitMPI"<<std::endl;
    double tare = MPI_Wtime();
    for(auto &node : mpiHSL)
    {
        node >> eQ::mpi::method::WAIT_FOR_ACKNOWLEDGE;
    }
    mpiController >> eQ::mpi::method::WAIT_FOR_ACKNOWLEDGE;
    waitTimer += (MPI_Wtime() - tare);
//    std::cout<<"Exiting waitMPI"<<std::endl;
}

void Simulation::writeHSLFiles(double simTime)
{
    if((isDiffusionNode) && createdHSLgrid)
    {
//        std::cout<<"writing hsl files..."<<std::endl;
        diffusionSolver->writeDiffusionFiles(simTime);

//        if(bool(eQ::data::parameters["PETSC_SIMULATION"]))
//        {
//            diffusionSolver2->writeDiffusionFiles(simTime);
//        }

        //print the l2 norm and boundary flux of the HSL:
//        auto thisNorm = diffusionSolver->u->vector()->norm("l2");
//        auto flux = diffusionSolver->getBoundaryFlux();
        std::cout<<"\t boundaryFlux assemble result: ("<<my_PE_num<<") "<<diffusionSolver->totalBoundaryFlux
//                    <<"  u norm: "<<thisNorm
                   <<" well conc.: "<<boundaryWellConcentration
//                  <<" boundary underflow counter: "<<boundaryUnderFlow
//                 <<" boundary finite counter: "<<boundaryFinite
                  <<std::endl;

    }
}
//TODO:  GENERALIZE THIS FOR EACH TYPE OF GENE CIRCUIT
void Simulation::create_DataGrid()
{//for using Fenics for data recording only (populating grid functions)
    if(isControllerNode)//only controller makes data grid:
    {
        computeGridParameters();
        //create the fenics interface class:
        //**************************************************************
        //              CREATE FENICS CLASS:
        //**************************************************************
        fenicsParams.comm       = controllerComm;//single mpi node for the controller
        fenicsParams.dataFiles = params.dataFiles;
        fenicsParams.nodesPerMicron = double(eQ::data::parameters["nodesPerMicronData"]);
        //uses non-template version of diffusion class for fenics
        diffusionSolver = std::make_shared<fenicsInterface>(fenicsParams);
        std::cout<<"Fenics for data recording Class created..."<<std::endl;
    }
}
void Simulation::writeDataFiles(double simTime)
{//when using fenics to record cell data over the grid (not for HSL data here)
    if(params.dataFiles->empty()) return;
    if(isControllerNode)
    {
        diffusionSolver->writeDataFiles(simTime);
    }
}
void Simulation::computeBoundaryWell()
{
    double flux    = diffusionSolver->totalBoundaryFlux;

    //COMPUTE TRAP CHANNEL CONCENTRATION:
    boundaryWellConcentration
            += flux/wellScaling;//scale #c from flux integral
    boundaryWellConcentration
            -= double(eQ::data::parameters["dt"])*boundaryDecayRate*boundaryWellConcentration;

    if(boundaryWellConcentration < 0.0)
    {
        boundaryWellConcentration = 0.0;
        boundaryUnderFlow++;
    }
    if("DIRICHLET_UPDATE" == eQ::data::parameters["boundaryType"])
    {
        //note:  set "allBoundaries" to false if not using all boundaries with same value!
        eQ::data::parametersType bvals;
//        bvals["allBoundaries"] = bool(true);
//        bvals["boundaryValue"] = boundaryWellConcentration;
        diffusionSolver->setBoundaryValues(boundaryWellConcentration);
    }
}
void Simulation::setBoundaryRate()
{
    boundaryDecayRate =  //units: min^-1; = linearFlowRate/trapLength = 100/2000=0.05sec^-1
            double(eQ::data::parameters["simulationFlowRate"])
            / double(eQ::data::parameters["simulationTrapWidthMicrons"]);

    if((double(eQ::data::parameters["dt"]) * boundaryDecayRate) >= 1.0)
        std::cout<<"ERROR!  timestep too large for decay of boundary well concentration: "
                <<(double(eQ::data::parameters["dt"]) * boundaryDecayRate)<<std::endl;
}
void Simulation::initBoundaryWell()
{
    //BOUNDARY WELL MODEL INITIALIZATION
    //initialization w.r.t. the trap boundary (if used)
    boundaryWellConcentration = 0.0;
    boundaryUnderFlow = 0;

    //compute volume of flow channels (+ left/right sides in first approximation)
    wellScaling = 2.0 * 10.0 * (15.0/double(eQ::data::parameters["lengthScaling"])) //10um z-height, 15um y-height
            * (double(eQ::data::parameters["simulationTrapWidthMicrons"]) + double(eQ::data::parameters["simulationTrapHeightMicrons"]));//x channelThickness x channelHeight x #

    setBoundaryRate();
}
void Simulation::mpiAssignCommunicators()
{
    //ITERATE THROUGH EACH DIFFUSION LAYER RANK (1...N) FOR N+1 TOTAL RANKS (INCLUDING CONTROLLER RANK=0)
    int thisRank=0;//start with controller
    for(size_t i(0); i<numHSLGrids; ++i)
    {
        mpiRanks.push_back(std::vector<int>());//creates mpiRanks[i] vector
        for(size_t j(0); j<mpiNodesPerDiffusionLayer; ++j)
        {
            mpiRanks[i].push_back(++thisRank);//starts at rank=1 for grid=0
            if(my_PE_num == thisRank)
            {
                whichHSL = i;//index from 0; defines which HSL layer for each mpi process
//                std::cout<<"Node "<<my_PE_num<<" belongs to HSL layer "<<whichHSL<<std::endl;
            }
        }
    }

    //Generate workers MPI group and communicator
    std::vector<int> ranks = {0};//list initialized vector (0 = controller node)

    //create two MPI comms: one with server only (the controller), the other with all but server (the workers):
    MPI_Comm_group(world, &world_group);                                            //extract the group from the comm
    MPI_Group_excl(world_group, ranks.size(), ranks.data(), &worker_group);         //create a new group, which excludes the controller
    MPI_Group_incl(world_group, ranks.size(), ranks.data(), &controller_group);     //create a new group, only the controller

    MPI_Comm_create(world, worker_group, &workers);                                 //create a new comm from new group
    MPI_Comm_create(world, controller_group, &controllerComm);                      //create a new comm from new group

    //now form sub-comms for each HSL grid to pass to the diffusion solver
    for(size_t i(0); i<numHSLGrids; ++i)
    {
        MPI_Group   thisGroup;
        MPI_Comm    thisComm;
        MPI_Group_incl(world_group, mpiRanks[i].size(), mpiRanks[i].data(), &thisGroup);//create a new group with HSL layer[i] world ranks
        MPI_Comm_create(world, thisGroup, &thisComm);                                   //create a new comm from new group (uses world ranks)
        mpiComms.push_back(thisComm);
        MPI_Group_free(&thisGroup);
    }

    MPI_Group_free(&worker_group);//free the groups
    MPI_Group_free(&controller_group);
    MPI_Group_free(&world_group);
}

Simulation::~Simulation()
{
    std::cout<<"BEGIN Simulation::~Simulation()"<<std::endl;

//**************************************************************************//
//                  SHUTDOWN MPI SEQUENCE:
//**************************************************************************//
    MPI_Barrier(world);
    if( (bool(eQ::data::parameters["hslSignaling"])) && !isControllerNode)
    {
        diffusionSolver->finalize();
        if(bool(eQ::data::parameters["PETSC_SIMULATION"]))
        {
            diffusionSolver2->finalize();
        }
    }

    diffusionSolver2.reset();
    diffusionSolver.reset();

    MPI_Barrier(world);

    //since we initialize MPI ourselves, we must call this for fenics:
    dolfin::SubSystemsManager::singleton().finalize();

    if(world != workers)
        if(MPI_COMM_NULL != workers)
            MPI_Comm_free(&workers);
    if(world != controllerComm)
        if(MPI_COMM_NULL != controllerComm)
            MPI_Comm_free(&controllerComm);


    std::cout<<"END Simulation::~Simulation()"<<std::endl;
    MPI_Barrier(world);
}
void Simulation::printData()
{
    for(size_t rank=0; rank<numHSLGrids; rank++)
    {
        std::cout<<"GRID #"<<rank<<std::endl;
        std::cout<<"gridCoords"<<rank<<std::endl;
        for(auto &coord : gridCoords[rank])
            std::cout<< std::setw(3)<<std::left
                <<coord
                  <<" ";
        std::cout<<std::endl;
        std::cout<<"gridDofs"<<rank<<std::endl;
        for(auto &dof : gridDofs[rank])
            std::cout<< std::setw(3)<<std::left
                <<dof
                  <<" ";
        std::cout<<std::endl;

        std::cout<<"ABM->dofLookupTable"<<rank<<std::endl;
        for(int ii=int(globalNodesH)-1; ii>=0; ii--)//descend for pretty-printing
        {
            size_t i = size_t(ii);
            for(size_t j=0; j<globalNodesW; j++)
            {
                auto point = eQ::nodePoint{i,j};
                std::cout<< std::setw(3)<<std::left
                    <<ABM->dofLookupTable[rank]->operator[](point)
                      <<" ";
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
    loadParams(eQ::Cell::strainType::ACTIVATOR,pA,rates);
    loadParams(eQ::Cell::strainType::REPRESSOR,pR,rates);

    if("DUALSTRAIN_OSCILLATOR" == eQ::data::parameters["simType"])
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
/*
void Simulation::create_DSOgrid()
{
    //only the controller node should call this:
    //note:  without ABM, the controller node will also be an hsl node.
    computeGridParameters();
    init_kinetics_DSO();

    std::cout<<"Creating dsoGrid of size H, W (dt):"<<globalNodesH<<", "<<globalNodesW
            <<" ("<<eQ::data::parameters["dt"]<<")"<<std::endl;

    //CREATE 2D GRID OF STRAIN POINTERS AND INSTANTIATE WITH STRAIN TYPE:
    //for DSO, use two grids (one activator, one repressor at each lattice point)
    dsoGrid.push_back(std::make_shared<eQ::gridFunction<std::shared_ptr<Strain>>>(globalNodesH, globalNodesW));
    dsoGrid.push_back(std::make_shared<eQ::gridFunction<std::shared_ptr<Strain>>>(globalNodesH, globalNodesW));
    for(size_t i(0); i<globalNodesH; i++)
    {
        for(size_t j(0); j<globalNodesW; j++)
        {
            dsoGrid[0]->grid[i][j] =
                        std::make_shared<Strain>(eQ::Cell::strainType::ACTIVATOR, &pA,
                                                 eQ::data::parameters["dt"], eQ::data::parameters["nodesPerMicronSignaling"]);
            dsoGrid[1]->grid[i][j] =
                        std::make_shared<Strain>(eQ::Cell::strainType::REPRESSOR, &pR,
                                                 eQ::data::parameters["dt"], eQ::data::parameters["nodesPerMicronSignaling"]);
        }
    }

    simulateABM = false;
}
*/
