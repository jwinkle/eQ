#include "diffuclass.h"
#include <iostream>
#include <math.h>
#include <boost/filesystem.hpp>

//------------------------------------------------------------------------------
//                           Declare functions
//------------------------------------------------------------------------------

extern PetscErrorCode ComputeMatrix(KSP, Mat, Mat, void*);

PetscErrorCode MyMatMult(Mat,Vec,Vec);


//------------------------------------------------------------------------------
//							  Public Methods
//------------------------------------------------------------------------------
//TODO: include all boundary condition elements and individually parallelized grids
void diffusionPETSc::initDiffusion(eQ::diffusionSolver::params &initParams)
//void diffusionPETSc::initDiffusion(MPI_Comm comm, std::vector<std::string> filePaths, int argc, char* argv[])
{
    auto argc = initParams.argc;
    auto argv = initParams.argv;
    int 	 i, j;

    /*
    //Split input communicator into subcommunicators for each desired diffusion grid
    MPI_Comm subComm;
    int 	 myRankMPI, commSize, i, j;

    MPI_Comm_size(comm, &commSize);
    MPI_Comm_rank(comm, &myRankMPI);

    MPI_Comm_split(comm, myRankMPI, 0, &subComm);

    initData.subCommunicator = subComm;
    //Record vector of diffusion constants and choose the appropriate constant for each processor
    std::vector<double>		d_vector;
    d_vector = std::vector<double>(eQ::parameters["D_HSL"].get<std::vector<double>>());
    initData.diffusionConstant = d_vector.at(myRankMPI % d_vector.size());

    //Set domain size based on parameters
    initData.xLengthMicrons = int(eQ::parameters["simulationTrapWidthMicrons"]);
    initData.yLengthMicrons = int(eQ::parameters["simulationTrapHeightMicrons"]);

    //Set time step
    initData.dt = eQ::parameters["dt"];

    //Set grid spacing in domain
    initData.h = 1 / double(eQ::parameters["nodesPerMicronSignaling"]);

    //Set the appropriate filepath to each grid for data recording
    initData.directoryName = filePaths.at(myRankMPI % d_vector.size());
*/
    initData.subCommunicator = initParams.comm;
    initData.diffusionConstant = initParams.D_HSL;//direct read from init params
    initData.xLengthMicrons = initParams.trapWidthMicrons;
    initData.yLengthMicrons = initParams.trapHeightMicrons;
    initData.dt = initParams.dt;
    initData.h = 1.0 / initParams.nodesPerMicron;

    initData.directoryName = initParams.filePath + "petsc";
     boost::filesystem::path dstFolder = initData.directoryName;
     boost::filesystem::create_directory(dstFolder);

     initData.objectName = "grid";

    //Set initial boundary conditions
    if("DIRICHLET_0" == eQ::parameters["boundaryType"]){
        initData.homogeneousDirichlet = PETSC_TRUE;

        initData.topDirichletCoefficient = 1;
        initData.bottomDirichletCoefficient = 1;
        initData.leftDirichletCoefficient = 1;
        initData.rightDirichletCoefficient = 1;

        initData.topNeumannCoefficient = 0;
        initData.bottomNeumannCoefficient = 0;
        initData.leftNeumannCoefficient = 0;
        initData.rightNeumannCoefficient = 0;

        initData.topBoundaryValue = 0;
        initData.bottomBoundaryValue = 0;
        initData.leftBoundaryValue = 0;
        initData.rightBoundaryValue = 0;
    }

    //Initialize PETSc element of diffusion
    InitializeDiffusion(&initData, argc, argv);

    //Resize and populate vectors for data transfer
    solution_vector.resize(gridNodesX * gridNodesY);
    allXCoordinates.resize(gridNodesX * gridNodesY);
    allYCoordinates.resize(gridNodesX * gridNodesY);

    for (j=0; j < gridNodesY; j++){
        for (i=0; i < gridNodesX; i++){
            allXCoordinates.at(i + j * gridNodesX) = i*initData.h;
            allYCoordinates.at(i + j * gridNodesX) = j*initData.h;
            solution_vector.at(i + j * gridNodesX) = 0;
        }
    }

    //Copy initial conditions to solution vector
    ReadGridValues(allXCoordinates, allYCoordinates, &solution_vector);
}

void diffusionPETSc::stepDiffusion(void)
{
    //Copy modified solution vector to u0
    //solution_vector is our public vector of all solution data - copy that over.

    WriteGridValues(allXCoordinates, allYCoordinates, solution_vector);

    TimeStep();

    ReadGridValues(allXCoordinates, allYCoordinates, &solution_vector);
}

//TODO: implement each wall separately
void diffusionPETSc::setBoundaryValues(const eQ::parametersType &bvals)
{
    if(bool(bvals["allBoundaries"]) == true){
        gridData->topBoundaryValue = double(bvals["boundaryValue"]);
        gridData->bottomBoundaryValue = double(bvals["boundaryValue"]);
        gridData->leftBoundaryValue = double(bvals["boundaryValue"]);
        gridData->rightBoundaryValue = double(bvals["boundaryValue"]);
    }
}

//TODO: verify accuracy, check about choosing walls
eQ::parametersType diffusionPETSc::getBoundaryFlux(void)
{
    //Compute and return boundary flux across all walls
    eQ::parametersType 	fluxData;
    double 			   	totalBoundaryFlux, boundarySlope;
    int 				i, j;

    for (i = 0; i < gridNodesX; i++){
        for (j = 0; j < gridNodesY; j++){
            if (i == 0){
                boundarySlope = boundarySlope + (solution_vector.at(2 + (gridNodesX) * j) - \
                        solution_vector.at((gridNodesX)*j))/(2*gridData->h);
            }
            else if (i == gridNodesX - 1){
                boundarySlope = boundarySlope + (solution_vector.at(i - 2 + (gridNodesX-1)*j) - \
                        solution_vector.at(i + (gridNodesX-1)*j))/(2*gridData->h);
            }
            else if (j == 0){
                boundarySlope = boundarySlope + (solution_vector.at(i + (gridNodesX-1) * 2) - \
                        solution_vector.at(i))/(2*gridData->h);
            }
            else if (j == gridNodesY - 1){
                boundarySlope = boundarySlope + (solution_vector.at(i + (gridNodesX-1) * (j - 2)) - \
                        solution_vector.at(i + (gridNodesX-1)*j))/(2*gridData->h);
            }
        }
    }

    totalBoundaryFlux = gridData->diffusionConstant * gridData->dt * boundarySlope;
    fluxData["totalFlux"] = totalBoundaryFlux;
    return fluxData;
}

void diffusionPETSc::writeDiffusionFiles(double timestamp)
{
    //Work in timestamp
    RecordData();
}

/*
void diffusionPETSc::writeDataFiles(double dt)
{
    //Need to write full function
}
*/

void diffusionPETSc::finalize(void)
{
    DiffusionFinalize(PETSC_TRUE);
}


//------------------------------------------------------------------------------
//                            Private Methods
//------------------------------------------------------------------------------

PetscErrorCode diffusionPETSc::ApplyBoundaryConditions()
{
  //Applys Neumann/Robin Boundary Conditions to RHS vector
  PetscInt       i,j,xm,ym,xs,ys;
  PetscScalar    **RHSarray;

  //Get subgrid information for each processor and get arrays from vector
  ierr = DMDAGetCorners(distributedArray,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(distributedArray, globalVector, &RHSarray);CHKERRQ(ierr);

  //JW NOTE:  each node can check its range whether it has boundary points
  //  then, traverse the boundary explicitly, rather than traversing the entire i,j space
    double twoFh = 2 * gridData->fourierNumber * gridData->h;//
    double topD = gridData->topDirichletCoefficient;
    double bottomD = gridData->bottomDirichletCoefficient;
    double leftD = gridData->leftDirichletCoefficient;
    double rightD = gridData->rightDirichletCoefficient;
        double topN = gridData->topNeumannCoefficient;
        double bottomN = gridData->bottomNeumannCoefficient;
        double leftN = gridData->leftNeumannCoefficient;
        double rightN = gridData->rightNeumannCoefficient;

  //Alter values in arrays to match Boundary Conditions
  for (j = ys; j < ys+ym; j++){
    for (i = xs; i < xs+xm; i++){
//      if (i == 0 || i == gridNodesX-1 || j == 0 || j == gridNodesY-1){
          if(j == gridNodesY-1)
          {
              double topBV = gridData->topBoundaryValue;
              //Set non-Dirichlet for top wall
              if (topN != 0)
              {//if TL (TR) corner, only set to N if left (right) wall is N
                  if ((i != 0 || leftN != 0) && (i != gridNodesX-1 || rightN != 0))
                        RHSarray[j][i] += (twoFh * topBV) / topN;
              }
              else
                  //JW: is topD ever 0??
                  RHSarray[j][i] = topBV / topD;
          }
          else if(j == 0)
          {
              double bottomBV = gridData->bottomBoundaryValue;
              //Set non-Dirichlet for bottom wall
              if (bottomN != 0)
              {
                  if ((i != 0 || leftN != 0) && (i != gridNodesX-1 || rightN != 0))
                        RHSarray[j][i] +=  (twoFh * bottomBV) / bottomN;
              }
              else
                  RHSarray[j][i] = bottomBV / bottomD;
          }
          if(i == gridNodesX-1)
          {
              double rightBV = gridData->rightBoundaryValue;
              //Set non-Dirichlet for right wall
              if (rightN != 0)
              {
                  if ((j != 0 || bottomN != 0) && (j != gridNodesY-1 || topN != 0))
                      RHSarray[j][i] += (twoFh * rightBV) / rightN;
              }
              else
                  RHSarray[j][i] = rightBV / rightD;
          }
          else if(i == 0)
          {
              double leftBV = gridData->leftBoundaryValue;
              //Set non-Dirichlet for left wall
              if (leftN != 0)
              {
                  if ((j != 0 || bottomN != 0) && (j != gridNodesY-1 || topN != 0))
                      RHSarray[j][i] += (twoFh * leftBV) / leftN;
              }
              else
                  RHSarray[j][i] = leftBV / leftD;
          }
//      }
    }
  }
/*
  //Alter values in arrays to match Boundary Conditions
  for (j = ys; j < ys+ym; j++){
    for (i = xs; i < xs+xm; i++){
      if (i == 0 || i == gridNodesX-1 || j == 0 || j == gridNodesY-1){
        //Set non-Dirichlet for top wall
        if (gridData->topNeumannCoefficient != 0 && j == gridNodesY-1){
          //Only change corners if they are also Neumann Boundary Conditions
          if ((i != 0 || gridData->leftNeumannCoefficient != 0) &&
              (i != gridNodesX-1 || gridData->rightNeumannCoefficient != 0)){
                RHSarray[j][i] = RHSarray[j][i] + (2 * gridData->fourierNumber * gridData->h * \
                  gridData->topBoundaryValue) / gridData->topNeumannCoefficient;
          }
        }
        //Set non-Dirichlet for bottom wall
        if (gridData->bottomNeumannCoefficient != 0 && j == 0){
          //Only change corners if they are also Neumann Boundary Conditions
          if ((i != 0 || gridData->leftNeumannCoefficient != 0) &&
              (i != gridNodesX-1 || gridData->rightNeumannCoefficient != 0)){
            RHSarray[j][i] = RHSarray[j][i] + (2 * gridData->fourierNumber * gridData->h * \
              gridData->bottomBoundaryValue) / gridData->bottomNeumannCoefficient;
          }
        }
        //Set non-Dirichlet for left wall
        if (gridData->leftNeumannCoefficient != 0 && i == 0){
          //Only change corners if they are also Neumann Boundary Conditions
          if ((j != 0 || gridData->bottomNeumannCoefficient != 0) &&
              (j != gridNodesY-1 || gridData->topNeumannCoefficient != 0)){
            RHSarray[j][i] = RHSarray[j][i] + (2 * gridData->fourierNumber * gridData->h * \
              gridData->leftBoundaryValue) / gridData->leftNeumannCoefficient;
          }
        }
        //Set non-Dirichlet for right wall
        if (gridData->rightNeumannCoefficient != 0 && i == gridNodesX-1){
          //Only change corners if they are also Neumann Boundary Conditions
          if ((j != 0 || gridData->bottomNeumannCoefficient != 0) &&
              (j != gridNodesY-1 || gridData->topNeumannCoefficient != 0)){
            RHSarray[j][i] = RHSarray[j][i] + (2 * gridData->fourierNumber * gridData->h * \
              gridData->rightBoundaryValue) / gridData->rightNeumannCoefficient;
          }
        }

        //Otherwise, set Dirichlet conditions
        if (gridData->topNeumannCoefficient == 0 && j == gridNodesY-1){
          RHSarray[j][i] = gridData->topBoundaryValue / gridData->topDirichletCoefficient;
        }
        else if (gridData->bottomNeumannCoefficient == 0 && j == 0){
          RHSarray[j][i] = gridData->bottomBoundaryValue / gridData->bottomDirichletCoefficient;
        }
        else if (gridData->leftNeumannCoefficient == 0 && i == 0){
          RHSarray[j][i] = gridData->leftBoundaryValue / gridData->leftDirichletCoefficient;
        }
        else if (gridData->rightNeumannCoefficient == 0 && i == gridNodesX-1){
          RHSarray[j][i] = gridData->rightBoundaryValue / gridData->rightDirichletCoefficient;
        }
      }
    }
  }
*/
  ierr = DMDAVecRestoreArray(distributedArray, globalVector, &RHSarray);CHKERRQ(ierr);

  return ierr;
}

//Figure out argc and argv
PetscErrorCode diffusionPETSc::InitializeDiffusion (DiffusionData *dataStruct, int argc, char **argv)
{
  //Initializes sets data values, inializes PETSc, and otherwise takes care of
  //all preparation needed in order to iterate. It should be noted that this
  //function does NOT set any initial conditions.

  gridData = dataStruct;
  //Initialize PETSc and assign passed communicator to DIFFU_COMM
  DIFFU_COMM = gridData->subCommunicator;
  ierr = PetscInitialize(&argc,&argv,(char*)0,NULL);
  MPI_Comm_size(DIFFU_COMM, &commSize);

  //Set up viewer for print statements
  printViewer = PETSC_VIEWER_STDOUT_(DIFFU_COMM);

  //Calculate crucial values for grid
  gridNodesX = (gridData->xLengthMicrons/gridData->h) + 1;
  gridNodesY = (gridData->yLengthMicrons/gridData->h) + 1;
  gridData->totalNodes = gridNodesX * gridNodesY;
  gridData->fourierNumber = (gridData->diffusionConstant*gridData->dt) / (gridData->h*gridData->h);

  //Set boundary and stencil types for distributed array
  boundary = DM_BOUNDARY_NONE;
  stencilType = DMDA_STENCIL_STAR;

  //Create Krylov Subspace method(KSP) context and Distributed Array(DMDA)
  ierr = KSPCreate(DIFFU_COMM, &krylovSolver);CHKERRQ(ierr);
  ierr = DMDACreate2d(DIFFU_COMM, boundary, boundary, stencilType, gridNodesX, gridNodesY, PETSC_DECIDE, \
    PETSC_DECIDE, 1, 1, NULL, NULL, &distributedArray);CHKERRQ(ierr);

  //Configure DMDA
  ierr = DMSetFromOptions(distributedArray);CHKERRQ(ierr);
  ierr = DMSetUp(distributedArray);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(distributedArray, 0, gridData->xLengthMicrons, 0, \
    gridData->yLengthMicrons, 0, 0);CHKERRQ(ierr);
  ierr = DMView(distributedArray, printViewer);CHKERRQ(ierr);

  //Create solution vector and local helper vector for scattering
  ierr = DMCreateGlobalVector(distributedArray, &globalVector);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(distributedArray, &gridData->localVector);CHKERRQ(ierr);

  //Obtain additional needed values for data structure
  ierr = VecGetOwnershipRange(globalVector, &gridData->localMinimum, &gridData->localMaximum);CHKERRQ(ierr);

  //Configure KSP
  ierr = KSPSetDM(krylovSolver, distributedArray);CHKERRQ(ierr);
  ierr = KSPSetComputeOperators(krylovSolver, ComputeMatrix, gridData);CHKERRQ(ierr);
  ierr = KSPGetPC(krylovSolver, &preconditioner);CHKERRQ(ierr);
  ierr = PCSetType(preconditioner, PCNONE);CHKERRQ(ierr);
  ierr = KSPSetType(krylovSolver, KSPFBCGSR);
  ierr = KSPSetFromOptions(krylovSolver);CHKERRQ(ierr);
  ierr = KSPSetUp(krylovSolver);CHKERRQ(ierr);

  //Prepare Application Ordering (AO)
  ierr = DMDAGetAO(distributedArray, &appOrder);CHKERRQ(ierr);

  ierr = PetscPrintf(DIFFU_COMM, "\nPETSC CLASS InitializeDiffusion(): %s \n", gridData->directoryName.c_str());CHKERRQ(ierr);

  //Set initial timestep number
  step = 0;
  return ierr;
}

PetscErrorCode diffusionPETSc::TimeStep()
{
  //Takes a step forward in time in the simulation
  step++;
  ierr = ApplyBoundaryConditions();CHKERRQ(ierr);
  ierr = KSPSolve(krylovSolver, globalVector, globalVector);CHKERRQ(ierr);
  ierr = PetscPrintf(DIFFU_COMM,"\r%d steps complete", step);CHKERRQ(ierr);
  return ierr;
}

//Ensure vector references work as desired
PetscErrorCode diffusionPETSc::WriteGridValues(vector<double> xCoordinates, vector<double> \
  yCoordinates, vector<double> values)
{
  //Takes values and their coordinates in vector form, converts coordinates to
  //a discretized form, and then inserts the values into the grid based on these
  //coordinates.

  //Declare temporary variables for use in method
  PetscInt      xm, ym, xs, ys, xindex, yindex;
  unsigned int  size, i;
  PetscScalar   **globalarray;

  //Record size variable and check for size consistency across vectors
  size = xCoordinates.size();
  if (size != yCoordinates.size() || size != values.size()){
    ierr = PetscPrintf(PETSC_COMM_SELF, \
        "Error: vector sizes for WriteGridValues must match\n");CHKERRQ(ierr);
  }

  //Get coordinate range and array for all data on processor
  ierr = DMDAGetCorners(distributedArray,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(distributedArray, globalVector, &globalarray);CHKERRQ(ierr);

  //Iterate through all coordinates, filter out those not on processor, and
  //insert values into grid
  for (i=0; i<size; i++){
    xindex = round(xCoordinates.at(i)/gridData->h);
    yindex = round(yCoordinates.at(i)/gridData->h);
    if (xindex >= xs && xindex < xs+xm && yindex >= ys && yindex < ys+ym){
        globalarray[yindex][xindex] = values.at(i);
    }
  }

  //Restore array to vector, finalizing the insert of these new values.
  ierr = DMDAVecRestoreArray(distributedArray, globalVector, &globalarray);CHKERRQ(ierr);

  return ierr;
}

//TODO: Ensure vector references work as desired
//Extensively test commSize > 1 case - need MPI to properly share data?
PetscErrorCode diffusionPETSc::ReadGridValues(vector<double> xCoordinates, vector<double> \
  yCoordinates, vector<double> *values)
{
  //Given some coordinates and a vector for values, discretizes the coordinates
  //and finds the values associated with them in the grid. Puts these values
  //into the given vector. Note that this only does so for on-process values.

  //Record size variable and check for size consistency across vectors
  unsigned int  readSize;
  readSize = xCoordinates.size();
  if (readSize != yCoordinates.size() || readSize != values->size()){
    ierr = PetscPrintf(PETSC_COMM_SELF, \
    "Error: vector sizes for ReadGridValues must match\n");CHKERRQ(ierr);
    return ierr;
  }

  //Read values from grid - use alternative algorithm if grid is parallelized
  if (commSize == 1){
    //Declare temporary variables for use in method
    PetscInt            xm, ym, xs, ys, xindex, yindex;
    unsigned int        i;
    PetscScalar         **readarray;


    //Get coordinate range and array for all data on processor
    ierr = DMDAGetCorners(distributedArray,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(distributedArray, globalVector, &readarray);CHKERRQ(ierr);

    //Get list of natural grid coordinates
    for (i=0; i < readSize; i++){
        xindex = round(xCoordinates.at(i)/gridData->h);
        yindex = round(yCoordinates.at(i)/gridData->h);
        if (xindex >= xs && xindex < xs+xm && yindex >= ys && yindex < ys+ym){
            values->at(i) = readarray[yindex][xindex];
        }
    }

    //Restore read array to vector
    ierr = DMDAVecRestoreArrayRead(distributedArray, globalVector, &readarray);CHKERRQ(ierr);
  }

  else{
    //Declare temporary variables for use in method
    PetscInt 		   	xindex, yindex, globalIndex[readSize],
                        localIndex[readSize];
    unsigned int		i;
    const PetscScalar	*readarray;
    IS					from, to;
    Vec					sequentialVector;
    VecScatter			scatter;

    //Create index vector for later PETSc vector scattering
    for (i=0; i < readSize; i++){
        localIndex[i] = i;
    }

    //Get list of natural grid coordinates
    for (i=0; i < readSize; i++){
        xindex = round(xCoordinates.at(i)/gridData->h);
        yindex = round(yCoordinates.at(i)/gridData->h);
        globalIndex[i] = yindex*gridNodesX + xindex;
    }

    //Replace natural grid coordinates with PETSc grid coordinates
    ierr = AOApplicationToPetsc(appOrder, readSize, globalIndex);CHKERRQ(ierr);

    //Prepares index sets for PETSc vector scattering
    ierr = ISCreateGeneral(DIFFU_COMM, readSize, globalIndex, PETSC_COPY_VALUES, &from);
        CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_SELF, readSize, localIndex, PETSC_COPY_VALUES, &to);
        CHKERRQ(ierr);

    //Create PETSc vector to be scattered to
    ierr = VecCreateSeq(PETSC_COMM_SELF, readSize, &sequentialVector);CHKERRQ(ierr);

    //Create PETSc vector scatter and execute it
    ierr = VecScatterCreate(globalVector, from, sequentialVector, to, &scatter);CHKERRQ(ierr);
    ierr = VecScatterBegin(scatter, globalVector, sequentialVector, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ(ierr);
    ierr = VecScatterEnd(scatter, globalVector, sequentialVector, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ(ierr);

    //Obtain data array from scattered PETSc vector
    ierr = VecGetArrayRead(sequentialVector, &readarray);CHKERRQ(ierr);

    //Copy data array into vector (non-PETSc vector)
    for (i=0; i < readSize; i++){
        values->at(i) = readarray[i];
    }

    //Restore array to vector
    ierr = VecRestoreArrayRead(sequentialVector, &readarray);CHKERRQ(ierr);

    //Destroy unnecessary objects
    ierr = ISDestroy(&from);CHKERRQ(ierr);
    ierr = ISDestroy(&to);CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter);CHKERRQ(ierr);
    ierr = VecDestroy(&sequentialVector);CHKERRQ(ierr);
  }

  return ierr;
}

PetscErrorCode diffusionPETSc::RecordData()
{
  //Writes grid data into a .vtr file to be read by ParaView

  //Declare temporary variables for use in method
  char filename[100];
  PetscViewer vtrviewer;

  //Get name of file to be written
  sprintf(filename, "%s/%s%04d.vtr", gridData->directoryName.c_str(), \
          gridData->objectName.c_str(), step);

    //    sprintf(filename, "%s/%s%d.vtr", gridData->directoryName.c_str(),

  //Open viewer and write data to file
  ierr = PetscViewerVTKOpen(DIFFU_COMM, filename, FILE_MODE_WRITE, \
    &vtrviewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(vtrviewer, PETSC_VIEWER_VTK_VTR);CHKERRQ(ierr);
  ierr = VecView(globalVector, vtrviewer);CHKERRQ(ierr);

  //Destroy viewer and return ierr
  ierr = PetscViewerDestroy(&vtrviewer);CHKERRQ(ierr);
  return ierr;
}

PetscErrorCode diffusionPETSc::DiffusionFinalize(PetscBool last)
{
    //Destroys all remaining objects in class and finalizes PETSc if "last" is set
    //to true.

    //Destroy objects created for an instance of this class
    ierr = DMDestroy(&distributedArray);CHKERRQ(ierr);
    ierr = KSPDestroy(&krylovSolver);CHKERRQ(ierr);
    ierr = VecDestroy(&globalVector);CHKERRQ(ierr);
    ierr = VecDestroy(&gridData->localVector);CHKERRQ(ierr);

    ierr = PetscPrintf(DIFFU_COMM, "\nDone\n");CHKERRQ(ierr);

    //If this is the last instance of this class to be finished, finalize PETSc.
    if (last){
        ierr = PetscFinalize();
    }

    return ierr;
}

//------------------------------------------------------------------------------
//                    Helper functions (external to class)
//------------------------------------------------------------------------------

PetscErrorCode ComputeMatrix(KSP krylovSolver, Mat A, Mat jac, void *user){
    //Establishes shell matrix for KSP operation

    //Initialize function---------------------------------------------------------
    //Declare variables
    PetscErrorCode ierr;
    DiffusionData  *gridData = (DiffusionData*)user;

    //Begin function
    PetscFunctionBegin;

    //Create matrices-------------------------------------------------------------
    //Form matrix for KSP process
    ierr = MatSetSizes(A, gridData->localMaximum-gridData->localMinimum, gridData->localMaximum-\
        gridData->localMinimum, gridData->totalNodes, gridData->totalNodes);CHKERRQ(ierr);
    ierr = MatSetType(A,MATSHELL);CHKERRQ(ierr);
    ierr = MatShellSetContext(A,gridData);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);

    //Set multiplication function for KSP process matrix
    ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void))MyMatMult);
        CHKERRQ(ierr);

    //Return error code-----------------------------------------------------------
    PetscFunctionReturn(ierr);
}

PetscErrorCode MyMatMult(Mat A, Vec X, Vec Y){
    //Defines the multiplication function for shell matrix

    //Initialize function---------------------------------------------------------
    //Declare variables
    PetscErrorCode ierr;
    void           *ptr;
    DiffusionData  *user;
    PetscScalar    **yarray, **xarray;
    DM             distributedArray;
    PetscInt       i,j,xm,ym,xs,ys,gxm,gym,gxs,gys;

    //Begin function
    PetscFunctionBegin;

    //Retrieve context and DM
    ierr = MatShellGetContext(A,&ptr);CHKERRQ(ierr);
    ierr = VecGetDM(X,&distributedArray);CHKERRQ(ierr);
    user = (DiffusionData*)ptr;

    //Prepare for computation-----------------------------------------------------
    //Begin scattering
    ierr = DMGlobalToLocalBegin(distributedArray, X, INSERT_VALUES, user->localVector);CHKERRQ(ierr);

    //Get local domain in grid
    ierr = DMDAGetCorners(distributedArray,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
    ierr = DMDAGetGhostCorners(distributedArray,&gxs,&gys,0,&gxm,&gym,0);CHKERRQ(ierr);

    //Get local array for product vector
    ierr = DMDAVecGetArray(distributedArray, Y, &yarray);CHKERRQ(ierr);

    //Finish scattering and get local array from scattered vector
    ierr = DMGlobalToLocalEnd(distributedArray, X, INSERT_VALUES, user->localVector);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(distributedArray, user->localVector, &xarray);CHKERRQ(ierr);

    //Compute product-------------------------------------------------------------
    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {
            //Sets value for edge points
            if (i == gxs || j == gys || i == gxs+gxm-1 || j == gys+gym-1){
                //Evaluate in case where top wall is Neumann/Robin
                if (user->topNeumannCoefficient != 0 && j == gys+gym-1){
                    //Evaluate left corner when left wall is Neumann/Robin
                    if (i == gxs && user->leftNeumannCoefficient != 0){
                        yarray[j][i] = -2*user->fourierNumber*xarray[j][i+1]-2*user->fourierNumber*xarray[j-1][i]\
                            +(1+(4 + (2*user->h*user->topDirichletCoefficient/user->topNeumannCoefficient)\
                            +(2*user->h*user->leftDirichletCoefficient/user->leftNeumannCoefficient))*\
                            user->fourierNumber)*xarray[j][i];
                    }
                    //Evaluate right corner when right wall is Neumann/Robin
                    else if (i == gxs+gxm-1 && user->rightNeumannCoefficient != 0){
                        yarray[j][i] = -2*user->fourierNumber*xarray[j][i-1]-2*user->fourierNumber*xarray[j-1][i]\
                            +(1+(4 + (2*user->h*user->topDirichletCoefficient/user->topNeumannCoefficient)\
                            +(2*user->h*user->rightDirichletCoefficient/user->rightNeumannCoefficient))*\
                            user->fourierNumber)*xarray[j][i];
                    }
                    //Evaluate all other edge points on top wall
                    else{
                        yarray[j][i] = -2*user->fourierNumber*xarray[j-1][i]-user->fourierNumber*xarray[j][i-1]\
                            -user->fourierNumber*xarray[j][i+1]+(1 + (4 + 2*user->h*user->topDirichletCoefficient\
                            /user->topNeumannCoefficient)*user->fourierNumber)*xarray[j][i];
                    }
                }
                //Evaluate in case where bottom wall is Neumann/Robin
                else if (user->bottomNeumannCoefficient != 0 && j == gys){
                    //Evaluate left corner when left wall is Neumann/Robin
                    if (i == gxs && user->leftNeumannCoefficient != 0){
                        yarray[j][i] = -2*user->fourierNumber*xarray[j][i+1]-2*user->fourierNumber*xarray[j+1][i]\
                            +(1+(4 + (2*user->h*user->leftDirichletCoefficient/user->leftNeumannCoefficient) + \
                            (2*user->h*user->bottomDirichletCoefficient/user->bottomNeumannCoefficient))*\
                            user->fourierNumber)*xarray[j][i];
                    }
                    //Evaluate right corner when right wall is Neumann/Robin
                    else if (i == gxs+gxm-1 && user->rightNeumannCoefficient != 0){
                        yarray[j][i] = -2*user->fourierNumber*xarray[j][i-1]-2*user->fourierNumber*xarray[j+1][i]\
                            +(1+(4 + (2*user->h*user->rightDirichletCoefficient/user->rightNeumannCoefficient) + \
                            (2*user->h*user->bottomDirichletCoefficient/user->bottomNeumannCoefficient))*\
                            user->fourierNumber)*xarray[j][i];
                    }
                    //Evaluate all other edge points on bottom wall
                    else{
                        yarray[j][i] = -2*user->fourierNumber*xarray[j+1][i]-user->fourierNumber*xarray[j][i-1]\
                            -user->fourierNumber*xarray[j][i+1]+(1+(4 + (2*user->h*user->bottomDirichletCoefficient\
                            /user->bottomNeumannCoefficient))*user->fourierNumber)*xarray[j][i];
                    }
                }
                //Evaluate in case where left wall is Neumann/Robin
                else if ((user->leftNeumannCoefficient != 0 && i == gxs) && (j != gys && j != gys+gym-1)){
                    yarray[j][i] = -user->fourierNumber*xarray[j-1][i]-user->fourierNumber*xarray[j+1][i]-\
                        2*user->fourierNumber*xarray[j][i+1]+(1+(4 + (2*user->h*user->leftDirichletCoefficient\
                        /user->leftNeumannCoefficient))*user->fourierNumber)*xarray[j][i];
                }
                //Evaluate in case where right wall is Neumann/Robin
                else if ((user->rightNeumannCoefficient != 0 && i == gxs+gxm-1) && (j != gys && j != gys+gym-1)){
                    yarray[j][i] = -user->fourierNumber*xarray[j-1][i]-user->fourierNumber*xarray[j+1][i]-\
                        2*user->fourierNumber*xarray[j][i-1]+(1+(4 + (2*user->h*user->rightDirichletCoefficient\
                        /user->rightNeumannCoefficient))*user->fourierNumber)*xarray[j][i];
                }
                //Evaluate wall in case of Dirichlet condition
                else{
                yarray[j][i] = xarray[j][i];
                }
            }
            //Sets value for non-edge points
            else{
                yarray[j][i] = 1.0 - user->fourierNumber*(
                        xarray[j-1][i] + xarray[j+1][i] + xarray[j][i-1] + xarray[j][i+1] - 4*xarray[j][i]);
                //        yarray[j][i] = -user->fourierNumber*xarray[j-1][i]-user->fourierNumber*xarray[j+1][i]-\
                //            user->fourierNumber*xarray[j][i-1] -user->fourierNumber*xarray[j][i+1]+\
                //            (1+4*user->fourierNumber)*xarray[j][i];
            }
        }
    }
    //Finish function-------------------------------------------------------------
    //Restore arrays
    ierr = DMDAVecRestoreArray(distributedArray, Y, &yarray);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(distributedArray, user->localVector, &xarray);CHKERRQ(ierr);

    //Return error code
    PetscFunctionReturn(ierr);
}
