#include <dolfin.h>

#include <random>
#include <cmath>
#include <fstream>
#include <csignal>
#include <iostream>
#include <chrono>

#include <mshr.h>
using namespace dolfin;
using namespace mshr;

#include "fHSL.h"

//public: (this is constructor with no arguments; => is an hsl fenics node)
fenicsInterface::fenicsInterface() :isDataRecordingNode(false)
{
    //switch on the diffusion node type desired:
    //TODO:  check parameter value
    //hard-code for the latest .ufl file uses:
    //shell is the pointer to the base class for access to the fenics-commom methods and variables.
    hslD =  std::make_shared<fenicsBaseClass<hslD::FunctionSpace, hslD::Form_L, hslD::Form_a>>(shell);
}
//for data recording only: (this is constructor with non-template initParms argument)
//calls the common init function with params set already (non-hsl node)
fenicsInterface::fenicsInterface(const eQ::diffusionSolver::params &initParams)
        : myParams(initParams), isDataRecordingNode(true)
{
    //the shell class consists of the mesh, functions, solvers etc...
    shell = std::make_shared<fenicsShell>();
    //the init function will distinguish signaling and data recording via boolean isDataRecordingNode
    fenicsClassInit();
}

//called from Simulation::create_HSLgrid (for hsl nodes only)
//void fenicsInterface::initDiffusion(MPI_Comm comm, std::vector<std::string> filePaths, int argc, char* argv[])
//void fenicsInterface::initDiffusion(size_t id, MPI_Comm comm, std::string filePath, double D, double dt, int argc, char* argv[])
void fenicsInterface::initDiffusion(eQ::diffusionSolver::params &initParams)
{
    myParams = initParams;
    //we convert the template init call into the internal (private) fenics init.
    std::cout<<"Template initDiffusion() called with mpi comm: "<<myParams.comm<<std::endl;
    fenicsClassInit();
    initHSLFiles();

    h = 1.0/double(eQ::parameters["nodesPerMicronSignaling"]);
    //compute volume of flow channel node: 10um z-height, 25um y-height, h wide
    wellScaling = 10.0 * (25.0/double(eQ::parameters["lengthScaling"])) * h;

    //buffers used to store total flux into channels per trap timestep:
    fluxTopChannel.assign(nodesW, 0.0);
    fluxBottomChannel.assign(nodesW, 0.0);

}
void fenicsInterface::setBoundaryValues(const eQ::parametersType &bvals)
{}

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
    for(size_t j(0); j<nodesW; ++j)
    {
        //LOWER BOUNDARY FLUX:
        dofc = bottomChannel->boundaryDofChannel[j];
        //finite-difference over 'h':
        doft1 = shell->dofLookupTable->grid[1][j];
        doft2 = shell->dofLookupTable->grid[0][j];
        gradc = (solution_vector[doft1] - solution_vector[doft2])/h;

        fluxBottomChannel[j] = computeFlux();

        //UPPER BOUNDARY FLUX:
            dofc = topChannel->boundaryDofChannel[j];
            //finite-difference over 'h':
            doft1 = shell->dofLookupTable->grid[nodesH-2][j];
            doft2 = shell->dofLookupTable->grid[nodesH-1][j];
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

    if( ("MICROFLUIDIC_TRAP" == eQ::parameters["boundaryType"]) && ("H_TRAP" != eQ::parameters["trapType"]) )
    {
        //now compute new flux into the channel boundaries and scale to channel volume element:
        computeBoundaryFlux();

        //copy the updated values into u0 and solve the next advection-diffusion timestep for both channels:

        size_t numIterations = size_t(eQ::parameters["channelSolverNumberIterations"]);
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
eQ::parametersType fenicsInterface::getBoundaryFlux(void)
{
    //TODO: implement flux computation out of each wall separately
    eQ::parametersType fluxData;
    fluxData["totalFlux"] = totalBoundaryFlux;
    return fluxData;
}
//private:
void fenicsInterface::createMesh(MPI_Comm comm)
{
    //Note: mesh geometry will be provided by simulation class
    //to sync with trap geometry in chipmunk;
    //switch on simulation type here to create mesh geometry:
    auto p0 = Point(0.0, 0.0);
    auto p1 = Point(myParams.trapWidthMicrons, myParams.trapHeightMicrons);
//    auto p1 = Point(eQ::parameters["simulationTrapWidthMicrons"], eQ::parameters["simulationTrapHeightMicrons"]);
//    auto p1 = Point(myParams.simData.trapWidthMicrons, myParams.simData.trapHeightMicrons);
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
    std::cout<<"I am "<<myRankMPI+1<<" of "<<myParams.filePath<<" fenics worker"<<std::endl;

    //for MPI:
    std::string old_ghost_mode = dolfin::parameters["ghost_mode"];
     dolfin::parameters["ghost_mode"] = "shared_facet";
//    dolfin::parameters["ghost_mode"] = "shared_vertex";

    ////////////////////////////////////////////////////////////////////////////////
    //                          MESH
    ////////////////////////////////////////////////////////////////////////////////

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
        //create new communicators, one per mpi node:
//        int zeroKey=0;
//        MPI_Comm_split(myParams.comm, myRankMPI, zeroKey, &mySingletonComm);
//        std::cout<<"Created new MPI comm! mySingletonComm = "<<mySingletonComm<<std::endl;
        MPI_Barrier(myParams.comm);

        //CREATE THE LOCAL MESH:
        nodesH = unsigned(ceil(myParams.trapHeightMicrons * myParams.nodesPerMicron));
        nodesW = unsigned(ceil(myParams.trapWidthMicrons * myParams.nodesPerMicron));

        //create trap mesh:
        createMesh(myParams.comm);//each mpi node creates its own mesh

        topChannel = std::make_shared<fenicsBaseClass<AdvectionDiffusion::FunctionSpace, AdvectionDiffusion::Form_L, AdvectionDiffusion::Form_a>>();
        bottomChannel = std::make_shared<fenicsBaseClass<AdvectionDiffusion::FunctionSpace, AdvectionDiffusion::Form_L, AdvectionDiffusion::Form_a>>();

        //TOP/BOTTOM CHANNELS MESH GENERATION:
        //use Robin BC instead of computing entire channel:
        auto uleft = double(0.0);
        auto uright = double(0.0);
//        auto uleft = double(eQ::parameters["channelLengthMicronsLeft"]);
//        auto uright = double(eQ::parameters["channelLengthMicronsRight"]);
        auto uwidth = double(eQ::parameters["simulationTrapWidthMicrons"]);
        nodesChannel = unsigned(
                    ceil( (double(eQ::parameters["simulationTrapWidthMicrons"]) + uleft + uright) * double(eQ::parameters["nodesPerMicronSignaling"]) ));

        //create channel mesh
        topChannel->mesh = std::make_shared<IntervalMesh>
                (myParams.comm, nodesChannel, -uleft, uwidth + uright);
        bottomChannel->mesh = std::make_shared<IntervalMesh>
                (myParams.comm, nodesChannel, -uleft, uwidth + uright);

        myDiffusionConstant = myParams.D_HSL;//use a single-entry vector to pass D
    }

    MPI_Barrier(myParams.comm);

    //actual number of nodes created is +1 in each dimension; save these as actual #nodes:
    nodesW++;
    nodesH++;
    nodesChannel++;

    if(isFenicsRootNode)
    {
        std::cout
                <<std::endl<<std::endl
                <<"Grid total nodes H,W: "<<nodesH<<" x "<<nodesW
                <<" = "<<nodesH * nodesW << " total nodes."
               <<" Channel nodes: "<<nodesChannel
               <<std::endl<<std::endl;
    }
    if(isDataRecordingNode) //no hsl, check for data
    {
        if(!myParams.dataFiles.empty())
            std::cout<<"Diffusion is empty, creating data writing grid only..."<<std::endl;

        //dataFiles is a vector of pair<dataSource, fileName>
        for(auto &file : myParams.dataFiles)
        {
            std::cout<<"Creating data file: "<<file.second<<std::endl;//uses operator<< in
            //create a scalar data function space and output file for each:
            std::string fileName = eQ::dataStrings[file.second];
            //create the fenics File class using MPI comm and filename passed in dataFiles vector
            outputFiles.push_back(std::make_shared<dolfin::File>(
                                      myParams.comm,  std::string(myParams.filePath + fileName), "compressed")
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
//        tensorDx->printSize();
//        tensorDy->printSize();

        //thin the output data:
        dolfin::set_log_level(dolfin::WARNING);
    }
}
//private:
void fenicsInterface::createHSL()
{
////////////////////////////////////////////////////////////////////////////////
//                      CLASS VARIABLES
////////////////////////////////////////////////////////////////////////////////
    fenicsVariable::data data;

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


//                              ======= MEDIA FLOW RATE SETTING =======
    double channelFlowVelocity = double(eQ::parameters["trapChannelLinearFlowRate"])/double(eQ::parameters["lengthScaling"]);
    data.v     = std::make_shared<dolfin::Constant>(channelFlowVelocity);//units: um/min


    double channelLengthLeft    = double(eQ::parameters["channelLengthMicronsLeft"]);//channel length is in simulation units
    double channelLengthRight    = double(eQ::parameters["channelLengthMicronsRight"]);//channel length is in simulation units
    double lvdl        = (channelLengthLeft * channelFlowVelocity)/myDiffusionConstant;
    double lvdr        = (channelLengthRight * channelFlowVelocity)/myDiffusionConstant;

    double rightRate = 0.0;
    double leftRate = 0.0;


    if(channelFlowVelocity > 1.0e-6)
    {
        leftRate = channelFlowVelocity * (1.0/(1.0 - exp(-lvdl)));
        rightRate = channelFlowVelocity * (1.0/(exp(lvdr) - 1.0));
    }
    else
    {//set to a linear slope to the channel left/right boundary:
        rightRate = myDiffusionConstant/channelLengthLeft;
        leftRate = myDiffusionConstant/channelLengthRight;
    }
    data.s_left     = std::make_shared<dolfin::Constant>(0.0);
    data.s_right     = std::make_shared<dolfin::Constant>(0.0);
    data.r_left     = std::make_shared<dolfin::Constant>(leftRate);
    data.r_right     = std::make_shared<dolfin::Constant>(rightRate);

////////////////////////////////////////////////////////////////////////////////
//                      FUNCTION SPACES, FORMS, FUNCTIONS:
//////////////////////////////////////////////////////////////////////////////
    //call the derived-class specific space and form creation methods

    createSpaceFormsFunctions();
        topChannel->createSpaceFormsFunctions();
        bottomChannel->createSpaceFormsFunctions();


    //CREATE CHANNEL ROBIN BOUNDARY CONDITION OBJECTS:
    //SUBDOMAINS TO MARK LEFT/RIGHT BOUNDARIES:
    auto leftWall       = std::make_shared<DirichletBoundary_TrapEdges>
            (myParams.trapHeightMicrons, myParams.trapWidthMicrons, DirichletBoundary_TrapEdges::edge::LEFT);
    auto rightWall      = std::make_shared<DirichletBoundary_TrapEdges>
            (myParams.trapHeightMicrons, myParams.trapWidthMicrons, DirichletBoundary_TrapEdges::edge::RIGHT);

    //MESH FUNCTION FOR TOP/BOTTOM CHANNELS: MARK LEFT=1, RIGHT=2, for use in .ufl file ds(1), ds(2)
    data.meshFunctionChannel = std::make_shared<MeshFunction<size_t>>(topChannel->mesh, topChannel->mesh->topology().dim()-1, 0);
    leftWall->mark(*data.meshFunctionChannel, 1);
    rightWall->mark(*data.meshFunctionChannel, 2);
    //set the .ufl form parameters from data:
        topChannel->setFormParameters(data);
        topChannel->createLinearVariationalSolver();
        topChannel->createCoordinatesToDofMappingChannel();

    //repeat for bottom channel
    data.meshFunctionChannel = std::make_shared<MeshFunction<size_t>>(bottomChannel->mesh, bottomChannel->mesh->topology().dim()-1, 0);
    leftWall->mark(*data.meshFunctionChannel, 1);
    rightWall->mark(*data.meshFunctionChannel, 2);
        bottomChannel->setFormParameters(data);
        bottomChannel->createLinearVariationalSolver();
        bottomChannel->createCoordinatesToDofMappingChannel();

    //actual #dofs are in Function sizes (excludes ghost cells)
    topChannel->u->vector()->get_local(solution_vectorTopChannel);
    solution_vectorTopChannel.assign(solution_vectorTopChannel.size(), 0.0);//initialize the vector to xfer data
        bottomChannel->u->vector()->get_local(solution_vectorBottomChannel);
        solution_vectorBottomChannel.assign(solution_vectorBottomChannel.size(), 0.0);//initialize the vector to xfer data
    //set additional vector for grid-ordering of channels (for petsc use)
    topChannelData.assign(solution_vectorTopChannel.size(), 0.0);
    bottomChannelData.assign(solution_vectorBottomChannel.size(), 0.0);

////////////////////////////////////////////////////////////////////////////////
//                    BOUNDARY CONDITIONS
////////////////////////////////////////////////////////////////////////////////

    //SUBDOMAIN CLASS INSTANTIATIONS:
    auto dbc_openWalled     = std::make_shared<DirichletBoundary_openWalls>();
    auto dbc_threeWalled    = std::make_shared<DirichletBoundary_threeWalls>();
    auto dbc_twoWalled      = std::make_shared<DirichletBoundary_twoWalls>(myParams.trapHeightMicrons, myParams.trapWidthMicrons);
    auto dbc_oneWall        = std::make_shared<DirichletBoundary_oneWall>(myParams.trapHeightMicrons, myParams.trapWidthMicrons);
    auto dbc_leftWall       = std::make_shared<DirichletBoundary_leftWall>(myParams.trapHeightMicrons, myParams.trapWidthMicrons);

    //PERIMETER BOUNDARY BC:

    //UPDATING DIRICHLET BOUNDARY:  switch to update expression class or zero

    if("MICROFLUIDIC_TRAP" == eQ::parameters["boundaryType"])
    {
        if("H_TRAP" == eQ::parameters["trapType"])
        {   //reset the rates to channel velocity for H-trap (left/right boundaries are narrow, transverse flow)
            //note: should be set with Robin BC set in "boundaries" parameters
            leftRate = channelFlowVelocity;
            rightRate = channelFlowVelocity;
        }
        auto leftWall       = std::make_shared<DirichletBoundary_TrapEdges>(myParams.trapHeightMicrons, myParams.trapWidthMicrons, DirichletBoundary_TrapEdges::edge::LEFT);
        auto rightWall      = std::make_shared<DirichletBoundary_TrapEdges>(myParams.trapHeightMicrons, myParams.trapWidthMicrons, DirichletBoundary_TrapEdges::edge::RIGHT);
        auto topWall        = std::make_shared<DirichletBoundary_TrapEdges>(myParams.trapHeightMicrons, myParams.trapWidthMicrons, DirichletBoundary_TrapEdges::edge::TOP);
        auto bottomWall     = std::make_shared<DirichletBoundary_TrapEdges>(myParams.trapHeightMicrons, myParams.trapWidthMicrons, DirichletBoundary_TrapEdges::edge::BOTTOM);

//        boundaryChannelLeft =  std::make_shared<updatingDirchletBoundary>(0.0);
//        boundaryChannelRight =  std::make_shared<updatingDirchletBoundary>(0.0);

        data.meshFunction = std::make_shared<MeshFunction<size_t>>(shell->mesh, shell->mesh->topology().dim()-1, 0);
        leftWall->mark(*data.meshFunction, 1);
        rightWall->mark(*data.meshFunction, 2);

//not sure that this works:
//        leftWall->mark_facets(*shell->mesh, 0);//translates to ds(0) in .ufl code
//        rightWall->mark_facets(*shell->mesh, 1);//translates to ds(1) in .ufl

        //GENERATE THE DIRICHLET BOUNDARY CONDITION OBJECT (space, values, subdomain)
        //use the channel solutions to populate the BC for the trap
        shell->dbc.clear();

        //DECODE THE BOUNDARY CONDITIONS:
        //LEFT/RIGHT WALL:
        std::vector<double> thisData;
        thisData = eQ::parameters["boundaries"]["left"][1].get<std::vector<double>>();
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
        thisData = eQ::parameters["boundaries"]["right"][1].get<std::vector<double>>();
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


        //TOP/BOTTOM WALL (SET TO FLOW CHANNEL):
        shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, topChannel->u, topWall));
        shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, bottomChannel->u, bottomWall));
        //SET TO 0
//        shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, data.zero, topWall));
//        shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, data.zero, bottomWall));

        if("H_TRAP" == eQ::parameters["trapType"])
        {  //reset the top,bottom boundaries to reflecting (no Dirichlet)
            shell->dbc.clear();
        }
    }
    else
    {
        //GENERATE THE DIRICHLET BOUNDARY CONDITION OBJECT (space, values, subdomain)
        std::shared_ptr<SubDomain> boundaryDomain;

        //switch on simulation/trap type here:  uses the expression subclass instance "bValues" as defined above
        if("NOWALLED" == eQ::parameters["trapType"]) boundaryDomain = dbc_openWalled;
        if("THREEWALLED" == eQ::parameters["trapType"]) boundaryDomain = dbc_threeWalled;
        if("TWOWALLED" == eQ::parameters["trapType"]) boundaryDomain = dbc_twoWalled;
        if("ONEWALLED" == eQ::parameters["trapType"]) boundaryDomain = dbc_oneWall;

        walls = std::make_shared<dolfin::DirichletBC>(shell->V, data.zero, boundaryDomain);

        if("DIRICHLET_0" == eQ::parameters["boundaryType"])
            shell->dbc.push_back(walls);
        else if("NEUMANN_3WALLED_TEST" == eQ::parameters["boundaryType"])
            shell->dbc.push_back(std::make_shared<dolfin::DirichletBC>(shell->V, data.zero, dbc_oneWall));
        else
        //Robin b.c.: (delete dirichlet BC):
            shell->dbc.clear();
    }

////////////////////////////////////////////////////////////////////////////////
//             FORMS, VARIATIONAL PROBLEMS AND SOLVERS::
////////////////////////////////////////////////////////////////////////////////
    //call the derived-class method to set parameters into the .ufl file
    setFormParameters(data);
    shell->createLinearVariationalSolver();
    shell->createGridCoordinatesToDofMapping();

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
//private:
void fenicsInterface::initHSLFiles()
{
    //input data is: filePath (from template, passed through to parameters) and mpiRank (for filename lookup)
    hslWriter.push_back(std::make_pair(//pair is (file, function)
                            std::make_shared<dolfin::File>(//a file needs passed mpi comm, file path+name string, and here type "compressed"
//                               mySingletonComm,
//                                std::string(myParams.filePaths[size_t(myRankMPI)]), "compressed"),
                               myParams.comm,
                                myParams.filePath, "compressed"),
                            shell->u));//"u" is the solution vector for hsl

    auto channelPath = std::string("./images/channel/")
            + std::to_string(size_t(myParams.comm)) + std::string("_") + std::to_string(size_t(myParams.uniqueID)) + std::string("_")
            + std::string("channelTop.pvd");
    hslWriter.push_back(std::make_pair(//pair is (file, function)
                            std::make_shared<dolfin::File>(//a file needs passed mpi comm, file path+name string, and here type "compressed"
                               myParams.comm,
                                channelPath, "compressed"),
                            topChannel->u));//"u" is the solution vector for hsl

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
        if(1 == myParams.dataFiles[i].first->getRank())
        {//VECTOR DATA
            auto expression = vectorDataExpression(myParams.dataFiles[i].first);
            uvec->interpolate(expression);
            *outputFiles[i] << std::pair<const Function*, double>(uvec.get(), dt);
        }
        else if(0 == myParams.dataFiles[i].first->getRank())
        {//SCALAR DATA
            auto expression = scalarDataExpression(myParams.dataFiles[i].first);
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
    //these must be cleared for MPI in fenics:
    for(auto &file : outputFiles)
        file.reset();

    for(auto &file : hslWriter)
        file.first.reset();

    boundaryFlux.reset();
    shell->mesh.reset();
    topChannel->mesh.reset();
    bottomChannel->mesh.reset();

//        shell.reset();
//        topChannel.reset();
//        bottomChannel.reset();
}
