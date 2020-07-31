#include <dolfin.h>

#include <random>
#include <cmath>
#include <fstream>
#include <csignal>
#include <iostream>
#include <chrono>

//#include <mshr.h>
//using namespace dolfin;
//using namespace mshr;

#include "fHSL.h"

//public: (this is constructor with no arguments; => is an hsl fenics node)
fenicsInterface::fenicsInterface() :isDataRecordingNode(false)
{
    //switch on the diffusion node type desired:
    //TODO:  check parameter value
    //hard-code for the latest .ufl file uses:
    //shell is the pointer to the base class for access to the fenics-commom methods and variables.
    shell =  std::make_shared<fenicsBaseClass<hslD::FunctionSpace, hslD::Form_L, hslD::Form_a>>();
}
//for data recording only: (this is constructor with non-template initParms argument)
//calls the common init function with params set already (non-hsl node)
fenicsInterface::fenicsInterface(const eQ::diffusionSolver::params &initParams)
        : myParams(initParams), isDataRecordingNode(true)
{
    //the shell class consists of the mesh, functions, solvers etc...
    shell = std::make_shared<fenicsData>();
    //the init function will distinguish signaling and data recording via boolean isDataRecordingNode
    fenicsClassInit();
}

//called from Simulation::create_HSLgrid (for hsl nodes only)
void fenicsInterface::initDiffusion(eQ::diffusionSolver::params &initParams)
{
    myParams = initParams;
    //we convert the template init call into the internal (private) fenics init.
    fenicsClassInit();
    initHSLFiles();

    //NOTE:  h is in simulation units (not physical units)
    h = 1.0/double(eQ::data::parameters["nodesPerMicronSignaling"]);
    //compute volume of flow channel node: 10um z-height, 25um y-height, h wide
    wellScaling = 10.0 * (25.0/double(eQ::data::parameters["lengthScaling"])) * h;

    //buffers used to store total flux into channels per trap timestep:
    fluxTopChannel.assign(nodesW, 0.0);
    fluxBottomChannel.assign(nodesW, 0.0);

}
void fenicsInterface::computeBoundaryFlux()
{
    //TODO: implement each wall separately
    //TODO: check over-ride "Dirichlet_0" is not set
//    if(bool(bvals["allBoundaries"]) == true);

    size_t dofc, doft1, doft2;
    double ds = 1.0 * h;//surface are through which boundary flux passes (1um high x grid width)
    double gradc;

    //lambda function for flux computation:
    auto computeFlux = [&]()
    {
        return ds*(myParams.dt * myDiffusionConstant * gradc)/wellScaling;
    };
    size_t j;
    using point = std::pair<eQ::nodeType, eQ::nodeType &>;
    point point0{0,j}, point1{1,j}, pointN1{nodesH-1,j}, pointN2{nodesH-2,j};
    for(j=0; j<nodesW; ++j)
    {
        //LOWER BOUNDARY FLUX:
        dofc = bottomChannel->boundaryDofChannel[j];
        //finite-difference over 'h':
//        doft1 = shell->dofLookupTable->grid[1][j];
//        doft2 = shell->dofLookupTable->grid[0][j];
        doft1 = shell->dofLookupTable->operator[](point1);
        doft2 = shell->dofLookupTable->operator[](point0);
        gradc = (solution_vector[doft1] - solution_vector[doft2])/h;

        fluxBottomChannel[j] = computeFlux();

        //UPPER BOUNDARY FLUX:
            dofc = topChannel->boundaryDofChannel[j];
            //finite-difference over 'h':
//            doft1 = shell->dofLookupTable->grid[nodesH-2][j];
//            doft2 = shell->dofLookupTable->grid[nodesH-1][j];
            doft1 = shell->dofLookupTable->operator[](pointN2);
            doft2 = shell->dofLookupTable->operator[](pointN1);
            gradc = (solution_vector[doft1] - solution_vector[doft2])/h;

            fluxTopChannel[j] = computeFlux();
    }
}

void fenicsInterface::stepDiffusion()
{
    //note: trap top and bottom boundary conditions are in the solution to the prevous advection/diffusion solution of the channel


    //copies modified solution_vector to u0; solution_vector updated by ABM via mpi transfer    /*    /*
    shell->u0->vector()->set_local(solution_vector);//sets u0 to contents of solution_vector
    //solve the system (rate-limiting step)
    shell->LVS->solve();
    //sets solution_vector to contents of the new solution u (solution_vector is sent to, and then updated by, the ABM layer)
    shell->u->vector()->get_local(solution_vector);

    if( ("MICROFLUIDIC_TRAP" == eQ::data::parameters["boundaryType"]) && ("H_TRAP" != eQ::data::parameters["trapType"]) )
    {
        //now compute new flux into the channel boundaries and scale to channel volume element:
        computeBoundaryFlux();

        //copy the updated values into u0 and solve the next advection-diffusion timestep for both channels:

        size_t numIterations = size_t(eQ::data::parameters["channelSolverNumberIterations"]);
        auto dtx = std::make_shared<dolfin::Constant>(myParams.dt/double(numIterations));
        topChannel->L->dt = dtx;
        topChannel->a->dt = dtx;
        bottomChannel->L->dt = dtx;
        bottomChannel->a->dt = dtx;

        for(size_t i(0); i<numIterations; ++i)
        {
            for(size_t j(0); j<nodesW; ++j)
            {
                auto dofc = topChannel->boundaryDofChannel[j];
                solution_vectorTopChannel[dofc] += fluxTopChannel[j]/double(numIterations);
            }
            for(size_t j(0); j<nodesW; ++j)
            {
                auto dofc = bottomChannel->boundaryDofChannel[j];
                solution_vectorBottomChannel[dofc] += fluxBottomChannel[j]/double(numIterations);
            }

            topChannel->u0->vector()->set_local(solution_vectorTopChannel);//sets u0 to contents of solution_vector
            topChannel->LVS->solve();
            topChannel->u->vector()->get_local(solution_vectorTopChannel);
                bottomChannel->u0->vector()->set_local(solution_vectorBottomChannel);//sets u0 to contents of solution_vector
                bottomChannel->LVS->solve();
                bottomChannel->u->vector()->get_local(solution_vectorBottomChannel);
        }
        //copy vector to petsc exportable vector (grid ordering)
        for(size_t j(0); j<nodesW; ++j)
        {
            auto dofc = topChannel->boundaryDofChannel[j];
            topChannelData[j] = solution_vectorTopChannel[dofc];
            dofc = bottomChannel->boundaryDofChannel[j];
            bottomChannelData[j] = solution_vectorBottomChannel[dofc];
        }
    }

    //old boundary method:
    //compute boundary flux and scale to channel volume
    ub->interpolate(*shell->u);
    flux->assemble(*boundaryFlux, *V0);
    //add new fluxed-out value scaled by D, dt, and channel volume
    //new value = D*dt*flux/well Volume
    totalBoundaryFlux = (myDiffusionConstant * myParams.dt *  boundaryFlux->get_scalar_value());
}

//private:
void fenicsInterface::createMesh(MPI_Comm comm)
{
    //Note: mesh geometry will be provided by simulation class
    //to sync with trap geometry in chipmunk;
    //switch on simulation type here to create mesh geometry:
    auto p0 = Point(0.0, 0.0);
    auto p1 = Point(myParams.trapWidthMicrons, myParams.trapHeightMicrons);
    shell->mesh = std::make_shared<RectangleMesh>
            (comm, p0, p1, nodesW, nodesH, "right");

//    meshCreated = true;

//    double rlimit = 4.0;
//    double hlimit = 1.0;
//    double hby2 = hlimit/2.0;

//    auto domain = mshr::Rectangle(dolfin::Point(0., 0.), dolfin::Point(rlimit, hlimit));
//    auto cdomain = mshr::Circle(Point(0.0, hby2), hby2);
//    auto pdomain = mshr::Polygon(std::vector<dolfin::Point>{
//                                     Point(rlimit,0.0), Point(rlimit, hlimit), Point(rlimit-hby2, hby2)
//                                 });
//    auto combo = (domain + cdomain) - pdomain;
//    mesh = mshr::generate_mesh(combo, 50.);

//    auto myMeshFunc = MeshFunction<size_t>(mesh, mesh->topology().dim()-1, mesh->domains());
//    std::cout<<myMeshFunc.str(true);

//    dolfin::info(dolfin::parameters, true);
}
//private:
//void fenicsInterface::fenicsClassInit(const fenicsInterface::params &initParams)
void fenicsInterface::fenicsClassInit()
{
    std::cout<<"fenicsInterfaceInit with myParams.comm = "<<myParams.comm<<std::endl;
    //this class will be created by worker nodes for HSL, and by the controller node for data recording
    //we init the PETSC comm ourselves and set our local rank
    //this class creates and interfaces the fenics objects via the "private" interface

    //NOTE: MPI COMM PASSED IN IS FOR A SINGLE DIFFUSION NODE
    //the size of the comm determines number nodes/layer

    //we must set this explicitly for MPI sub-communicators when using fenics:
    PETSC_COMM_WORLD = myParams.comm;
    MPI_Comm_size(myParams.comm, &numFenicsPEs);
    MPI_Comm_rank(myParams.comm, &myRankMPI);
    isFenicsRootNode = (0 == myRankMPI);
    std::cout<<"I am "<<myRankMPI+1<<" of "<<numFenicsPEs
            <<" recording data to file: "
            <<myParams.filePath
           <<std::endl;

    //for MPI:
    std::string old_ghost_mode = dolfin::parameters["ghost_mode"];
     dolfin::parameters["ghost_mode"] = "shared_facet";
//    dolfin::parameters["ghost_mode"] = "shared_vertex";

//================================================================
//              //MESH GENERATION
//================================================================

    if(isDataRecordingNode)//this is a data recording version launched by the controller node
    {
        if(1 != numFenicsPEs)
            std::cout<<"Error: data recorder asserted but 1 != numFenicsPEs..."
                    <<numFenicsPEs<<" != "<<1<<std::endl;
        //Note: nodesPerMicron is changed to that for data recording in simulation::create_DataGrid()
        nodesH = unsigned(ceil(myParams.trapHeightMicrons)*myParams.nodesPerMicron);
        nodesW = unsigned(ceil(myParams.trapWidthMicrons)*myParams.nodesPerMicron);

        //CREATE THE LOCAL MESH:
        createMesh(myParams.comm);//data recording node creates its own mesh

    }
    else
    {//SUCCESS ON FENICS DIFFUSION NODE INITIALIZATION:
        MPI_Barrier(myParams.comm);

        //CREATE THE LOCAL MESH:
        nodesH = unsigned(ceil(myParams.trapHeightMicrons * myParams.nodesPerMicron));
        nodesW = unsigned(ceil(myParams.trapWidthMicrons * myParams.nodesPerMicron));

        myDiffusionConstant = myParams.D_HSL;//use a single-entry vector to pass D

        //================================================================
        //              //TRAP MESH
        //================================================================
        //create trap mesh:
        createMesh(myParams.comm);//each mpi node creates its own mesh

        //================================================================
        //              //FLOW CHANNEL MESH
        //================================================================
        //pass the trap data structure for use by channels:
        topChannel      = std::make_shared<fenicsChannel>(shell->data);
        bottomChannel   = std::make_shared<fenicsChannel>(shell->data);

        //TOP/BOTTOM CHANNELS MESH GENERATION:
        //use Robin BC instead of computing entire channel:
        auto uleft = double(0.0);
        auto uright = double(0.0);
//        auto uleft = double(eQ::data::parameters["simulationChannelLengthLeft"]);
//        auto uright = double(eQ::data::parameters["simulationChannelLengthRight"]);
        auto uwidth = double(eQ::data::parameters["simulationTrapWidthMicrons"]);
        nodesChannel = unsigned(
                    ceil( (double(eQ::data::parameters["simulationTrapWidthMicrons"]) + uleft + uright) * double(eQ::data::parameters["nodesPerMicronSignaling"]) ));

        //create channel mesh
        topChannel->mesh = std::make_shared<IntervalMesh>
                (myParams.comm, nodesChannel, -uleft, uwidth + uright);
        bottomChannel->mesh = std::make_shared<IntervalMesh>
                (myParams.comm, nodesChannel, -uleft, uwidth + uright);

    }

    MPI_Barrier(myParams.comm);

    //actual number of nodes created is +1 in each dimension; save these as actual #nodes:
    nodesW++;
    nodesH++;
    nodesChannel++;

//================================================================
//              //CREATE HSL OBJECTS
//================================================================
    if(isDataRecordingNode) //no hsl, check for data
    {
        if(!myParams.dataFiles->empty())
            std::cout<<"Diffusion is empty, creating data writing grid only..."<<std::endl;

        for(auto &file : *myParams.dataFiles)
        {
            std::cout<<"Creating data file: "<<file.fileName<<std::endl;//uses operator<< in
            //create a scalar data function space and output file for each:
//            std::string fileName = eQ::dataStrings[file.second];
            //create the fenics File class using MPI comm and filename passed in dataFiles vector
            outputFiles.push_back(std::make_shared<dolfin::File>(
                                      myParams.comm,  std::string(myParams.filePath + file.fileName), "compressed")
                                  );
        }
        VscalarData = std::make_shared<data::FunctionSpace>(shell->mesh);
        shell->u   = std::make_shared<dolfin::Function>(VscalarData);  shell->u->set_allow_extrapolation(true);
        //vector data:
        VvectorData = std::make_shared<vdata::FunctionSpace>(shell->mesh);
        uvec = std::make_shared<dolfin::Function>(VvectorData); uvec->set_allow_extrapolation(true);

    }
    else
    {

        D11 = std::make_shared<std::vector<double>>();
        D22 = std::make_shared<std::vector<double>>();
        D12 = std::make_shared<std::vector<double>>();

        createHSL();//private fenics interface

        //Dxx are the (diagonal, off-diagonal) for the symmetric tensor "D", which can in general be anisotropic
        //init to 1,0 means diagonal is 1, off-diagonal is 0: isotropic D
        D11->assign(solution_vector.size(), 1.0);
        D22->assign(solution_vector.size(), 1.0);
        D12->assign(solution_vector.size(), 0.0);

        //thin the output data:
        dolfin::set_log_level(dolfin::WARNING);
    }
}

//public:
void fenicsInterface::setRobinBoundaryConditions()
{
    fenicsVariables &data = shell->data;

//      ======= MEDIA FLOW RATE SETTING =======
        channelFlowVelocity = double(eQ::data::parameters["simulationFlowRate"]);
        data.v     = std::make_shared<dolfin::Constant>(channelFlowVelocity);//units: um/min


        double channelLengthLeft    = double(eQ::data::parameters["simulationChannelLengthLeft"]);//channel length is in simulation units
        double channelLengthRight    = double(eQ::data::parameters["simulationChannelLengthRight"]);//channel length is in simulation units
        double lvdl        = (channelLengthLeft * channelFlowVelocity)/myDiffusionConstant;
        double lvdr        = (channelLengthRight * channelFlowVelocity)/myDiffusionConstant;



        //Robin boundary condition rates for left,right of trap (or channel)
        //see my Overleaf writeup for this derivation:
        if(channelFlowVelocity > 1.0e-6)
        {
            leftRate = channelFlowVelocity * (1.0/(1.0 - exp(-lvdl)));
            rightRate = channelFlowVelocity * (1.0/(exp(lvdr) - 1.0));
        }
        else
        {//set to a linear slope to the channel left/right boundary:
            leftRate = myDiffusionConstant/channelLengthLeft;
            rightRate = myDiffusionConstant/channelLengthRight;
        }
        data.s_left     = std::make_shared<dolfin::Constant>(0.0);
        data.s_right    = std::make_shared<dolfin::Constant>(0.0);
        data.r_left     = std::make_shared<dolfin::Constant>(leftRate);
        data.r_right    = std::make_shared<dolfin::Constant>(rightRate);

}

//private:
void fenicsInterface::createHSL()
{
    fenicsVariables &data = shell->data;

    //PERIMETER BOUNDARY BC:
    //SUBDOMAINS TO MARK LEFT/RIGHT BOUNDARIES:
    auto leftWall       = std::make_shared<DirichletBoundary_TrapEdge>(myParams.trapHeightMicrons, myParams.trapWidthMicrons, DirichletBoundary_TrapEdge::LEFT);
    auto rightWall      = std::make_shared<DirichletBoundary_TrapEdge>(myParams.trapHeightMicrons, myParams.trapWidthMicrons, DirichletBoundary_TrapEdge::RIGHT);
    auto topWall        = std::make_shared<DirichletBoundary_TrapEdge>(myParams.trapHeightMicrons, myParams.trapWidthMicrons, DirichletBoundary_TrapEdge::TOP);
    auto bottomWall     = std::make_shared<DirichletBoundary_TrapEdge>(myParams.trapHeightMicrons, myParams.trapWidthMicrons, DirichletBoundary_TrapEdge::BOTTOM);

    //================================================================
    //              //CLASS VARIABLES
    //================================================================

    data.dt     = std::make_shared<dolfin::Constant>(myParams.dt);
    data.D      = std::make_shared<dolfin::Constant>(myDiffusionConstant);
    data.zero   = std::make_shared<dolfin::Constant>(0.0);
    data.one    = std::make_shared<dolfin::Constant>(1.0);

    data.tensorD11 = std::make_shared<AnisotropicDiffusionTensor>(
                D11, size_t(myParams.nodesPerMicron), nodesH, nodesW, 1.0);
    data.tensorD22 = std::make_shared<AnisotropicDiffusionTensor>(
                D22, size_t(myParams.nodesPerMicron), nodesH, nodesW, 1.0);
    data.tensorD12 = std::make_shared<AnisotropicDiffusionTensor>(
                D12, size_t(myParams.nodesPerMicron), nodesH, nodesW, 0.0);


    setRobinBoundaryConditions();

//================================================================
//             FUNCTION SPACES, FORMS, FUNCTIONS:
//================================================================
    //call the derived-class specific space and form creation methods

    shell->createSpaceFormsFunctions();
        topChannel->createSpaceFormsFunctions();
        bottomChannel->createSpaceFormsFunctions();

//================================================================
//              //FLOW CHANNELS
//================================================================
    //CREATE CHANNEL ROBIN BOUNDARY CONDITION OBJECTS:

    //MESH FUNCTION FOR TOP/BOTTOM CHANNELS: MARK LEFT=1, RIGHT=2, for use in .ufl file ds(1), ds(2)
    data.meshFunctionChannel = std::make_shared<MeshFunction<size_t>>(topChannel->mesh, topChannel->mesh->topology().dim()-1, 0);
    leftWall->mark(*data.meshFunctionChannel, 1);
    rightWall->mark(*data.meshFunctionChannel, 2);
        topChannel->createLinearProblem();//copy data struct from shell

    //repeat for bottom channel
    data.meshFunctionChannel = std::make_shared<MeshFunction<size_t>>(bottomChannel->mesh, bottomChannel->mesh->topology().dim()-1, 0);
    leftWall->mark(*data.meshFunctionChannel, 1);
    rightWall->mark(*data.meshFunctionChannel, 2);
        bottomChannel->createLinearProblem();//copy data struct from shell

    //actual #dofs are in Function sizes (excludes ghost cells)
    topChannel->u->vector()->get_local(solution_vectorTopChannel);
    solution_vectorTopChannel.assign(solution_vectorTopChannel.size(), 0.0);//initialize the vector to xfer data
        bottomChannel->u->vector()->get_local(solution_vectorBottomChannel);
        solution_vectorBottomChannel.assign(solution_vectorBottomChannel.size(), 0.0);//initialize the vector to xfer data
    //set additional vector for grid-ordering of channels (for petsc use)
    topChannelData.assign(solution_vectorTopChannel.size(), 0.0);
    bottomChannelData.assign(solution_vectorBottomChannel.size(), 0.0);

//================================================================
//                    BOUNDARY CONDITIONS
//================================================================

    //SUBDOMAIN CLASS INSTANTIATIONS:
    auto dbc_openWalled     = std::make_shared<DirichletBoundary_openWalls>();
    auto dbc_threeWalled    = std::make_shared<DirichletBoundary_threeWalls>();
    auto dbc_twoWalled      = std::make_shared<DirichletBoundary_twoWalls>(myParams.trapHeightMicrons, myParams.trapWidthMicrons);
    auto dbc_oneWall        = std::make_shared<DirichletBoundary_oneWall>(myParams.trapHeightMicrons, myParams.trapWidthMicrons);
    auto dbc_leftWall       = std::make_shared<DirichletBoundary_leftWall>(myParams.trapHeightMicrons, myParams.trapWidthMicrons);


    //UPDATING DIRICHLET BOUNDARY:  switch to update expression class or zero

    if("MICROFLUIDIC_TRAP" == eQ::data::parameters["boundaryType"])
    {
        if("H_TRAP" == eQ::data::parameters["trapType"])
        {   //reset the rates to channel velocity for H-trap (left/right boundaries are narrow, transverse flow)
            //note: should be set with Robin BC set in "boundaries" parameters
            leftRate = channelFlowVelocity;
            rightRate = channelFlowVelocity;
        }

        data.meshFunction = std::make_shared<MeshFunction<size_t>>(shell->mesh, shell->mesh->topology().dim()-1, 0);
        leftWall->mark(*data.meshFunction, 1);
        rightWall->mark(*data.meshFunction, 2);

//not sure that this works:
//        leftWall->mark_facets(*shell->mesh, 0);//translates to ds(0) in .ufl code
//        rightWall->mark_facets(*shell->mesh, 1);//translates to ds(1) in .ufl

        shell->dbc.clear();

        //DECODE THE BOUNDARY CONDITIONS:

        //LEFT WALL:
        std::vector<double> thisData;
        thisData = eQ::data::parameters["boundaries"]["left"][1].get<std::vector<double>>();
        if(0.0 == thisData[0])
        {//DIRICHLET, SET TO BV:
            auto leftData = std::make_shared<dolfin::Constant>(thisData[2]);
            shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, leftData, leftWall));
        }
        else if(0.0 == thisData[1])
        {//NEUMANN, SET TO BV:
            //ONLY SET HOMOGENEOUS NEUMANN...
            data.s_left     = std::make_shared<dolfin::Constant>(0.0);
            data.r_left     = std::make_shared<dolfin::Constant>(0.0);
        }
        else
        {//ROBIN, only set to rate computed, above
            data.s_left     = std::make_shared<dolfin::Constant>(0.0);
            data.r_left     = std::make_shared<dolfin::Constant>(leftRate);
        }

        //RIGHT WALL:
        thisData = eQ::data::parameters["boundaries"]["right"][1].get<std::vector<double>>();
        if(0.0 == thisData[0])
        {//DIRICHLET, SET TO BV:
            auto rightData = std::make_shared<dolfin::Constant>(thisData[2]);
            shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, rightData, rightWall));//set to zero always
        }
        else if(0.0 == thisData[1])
        {//NEUMANN, SET TO BV:
            //ONLY SET HOMOGENEOUS NEUMANN...
            data.s_right     = std::make_shared<dolfin::Constant>(0.0);
            data.r_right     = std::make_shared<dolfin::Constant>(0.0);
        }
        else
        {//ROBIN:
            data.s_right     = std::make_shared<dolfin::Constant>(0.0);
            data.r_right     = std::make_shared<dolfin::Constant>(rightRate);
        }


        //TOP WALL
        thisData = eQ::data::parameters["boundaries"]["top"][1].get<std::vector<double>>();
        if(0.0 == thisData[0])//DIRICHLET
        {
            if(-1.0 == thisData[2])
            {
                //TOP/BOTTOM WALL (SET TO FLOW CHANNEL):
                shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, topChannel->u, topWall));
            }
            else
            {
                auto topData = std::make_shared<dolfin::Constant>(thisData[2]);
                shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, topData, topWall));
            }
        }
        //NEUMANN:  do nothing...

        //BOTTOM WALL
        thisData = eQ::data::parameters["boundaries"]["bottom"][1].get<std::vector<double>>();
        if(0.0 == thisData[0])//DIRICHLET
        {
            if(-1.0 == thisData[2])
            {
                //TOP/BOTTOM WALL (SET TO FLOW CHANNEL):
                shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, bottomChannel->u, bottomWall));
            }
            else
            {
                auto bottomData = std::make_shared<dolfin::Constant>(thisData[2]);
                shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, bottomData, bottomWall));
            }
        }
        //NEUMANN:  do nothing...
    }
    else
    {
        std::cout<<"\n\t Initializing boundary conditions via old method...should depracate\n"<<std::endl;

        if("DIRICHLET_UPDATE" == eQ::data::parameters["boundaryType"])
        {
                shell->dbc.clear();

                //COMPARTMENT IMPLEMENTATION:
                boundaryCompartment =  std::make_shared<updatingDirchletBoundary>(0.0);
                shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, boundaryCompartment, leftWall));//update value each timestep
                shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, boundaryCompartment, rightWall));
                shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, boundaryCompartment, topWall));//update value each timestep
                shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, boundaryCompartment, bottomWall));
        }
        //switch on simulation/trap type here:  uses the expression subclass instance "bValues" as defined above
        else if("DIRICHLET_0" == eQ::data::parameters["boundaryType"])
        {
            //GENERATE THE DIRICHLET BOUNDARY CONDITION OBJECT (space, values, subdomain)
            std::shared_ptr<SubDomain> boundaryDomain;
            if("NOWALLED" == eQ::data::parameters["trapType"]) boundaryDomain = dbc_openWalled;
            if("THREEWALLED" == eQ::data::parameters["trapType"]) boundaryDomain = dbc_threeWalled;
            if("TWOWALLED" == eQ::data::parameters["trapType"]) boundaryDomain = dbc_twoWalled;
            if("ONEWALLED" == eQ::data::parameters["trapType"]) boundaryDomain = dbc_oneWall;

            walls = std::make_shared<dolfin::DirichletBC>(shell->V, data.zero, boundaryDomain);

            shell->dbc.push_back(walls);
        }
        else if("NEUMANN_3WALLED_TEST" == eQ::data::parameters["boundaryType"])
            shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, data.zero, dbc_oneWall));
        else
            shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, data.zero, dbc_openWalled));
    }

////////////////////////////////////////////////////////////////////////////////
//             FORMS, VARIATIONAL PROBLEMS AND SOLVERS::
////////////////////////////////////////////////////////////////////////////////
    //call the derived-class method to set parameters into the .ufl file
    shell->createLinearProblem();

    //actual #dofs are in Function sizes (excludes ghost cells)
    shell->u->vector()->get_local(solution_vector);
    solution_vector.assign(solution_vector.size(), 0.0);//initialize the vector to xfer data

    if(isFenicsRootNode)
        std::cout<<"solution_vector.size() = "<<solution_vector.size()<<std::endl<<std::endl;


    //for computation of boundary flux for updating BC:
    boundaryFlux    = std::make_shared<dolfin::Scalar>(myParams.comm);
    V0      = std::make_shared<boundary::Functional>(shell->mesh);
    VU      = std::make_shared<boundary::CoefficientSpace_u>(shell->mesh);
    ub      = std::make_shared<Function>(VU); ub->set_allow_extrapolation(true);
    V0->u = ub;
    flux      = std::make_shared<dolfin::Assembler>();


}

void fenicsInterface::setBoundaryValues(const double bval)
{
    boundaryCompartment->updateBoundary(bval);
}
//private:
void fenicsInterface::initHSLFiles()
{
    //input data is: filePath (from template, passed through to parameters) and mpiRank (for filename lookup)
    hslWriter.push_back(std::make_pair(//pair is (file, function)
                            std::make_shared<dolfin::File>(//a file needs passed mpi comm, file path+name string, and here type "compressed"
                               myParams.comm,
                                myParams.filePath, "compressed"),
                            shell->u));//"u" is the solution vector for hsl

    if("MICROFLUIDIC_TRAP" == eQ::data::parameters["boundaryType"])
    {
       hslWriter.push_back(std::make_pair(//pair is (file, function)
                                std::make_shared<dolfin::File>(//a file needs passed mpi comm, file path+name string, and here type "compressed"
                                   myParams.comm,
                                    myParams.filePathTopChannel, "compressed"),
                                topChannel->u));//"u" is the solution vector for hsl
        hslWriter.push_back(std::make_pair(//pair is (file, function)
                                std::make_shared<dolfin::File>(//a file needs passed mpi comm, file path+name string, and here type "compressed"
                                   myParams.comm,
                                    myParams.filePathBottomChannel, "compressed"),
                                bottomChannel->u));//"u" is the solution vector for hsl
    }
}
//public:
void fenicsInterface::writeDiffusionFiles(double timeStamp)
{
    for(auto &file : hslWriter)
    {
        *file.first << std::pair<const Function*, double>(file.second.get(), timeStamp);
    }
}
void fenicsInterface::writeDataFiles(double dt)
{
    for(size_t i(0); i<outputFiles.size(); i++)
    {
        if(eQ::data::tensor::rank::VECTOR == myParams.dataFiles->at(i).data->getRank())
        {//VECTOR DATA
            auto expression = vectorDataExpression(myParams.dataFiles->at(i).data);
            uvec->interpolate(expression);
            *outputFiles[i] << std::pair<const Function*, double>(uvec.get(), dt);
        }
        else if(eQ::data::tensor::rank::SCALAR == myParams.dataFiles->at(i).data->getRank())
        {//SCALAR DATA
            auto expression = scalarDataExpression(myParams.dataFiles->at(i).data);
            shell->u->interpolate(expression);
            *outputFiles[i] << std::pair<const Function*, double>(shell->u.get(), dt);
        }
    }
}
void fenicsInterface::finalize()
{

}
fenicsInterface::~fenicsInterface()
{
}
