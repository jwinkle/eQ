#ifndef __DIFFUCLASS_H_INCLUDED 
#define __DIFFUCLASS_H_INCLUDED

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <vector>
#include <string>
//JJW mod:
//#include "eQ.h"
#include "./src/eQ.h"

using namespace std;


//Data structure for passing to various PETSc functions
typedef struct {
	//Unifying the interface - how to set non PETSc datatypes? 
	//Also, how to completely unify interface so as not to have 
	//PETSc-specific values in data?  
  PetscInt    xLengthMicrons, yLengthMicrons, localMinimum, localMaximum, totalNodes;

  PetscReal   diffusionConstant, h, dt, fourierNumber,

              topDirichletCoefficient, bottomDirichletCoefficient,
              leftDirichletCoefficient, rightDirichletCoefficient,

              topNeumannCoefficient, bottomNeumannCoefficient,
              leftNeumannCoefficient, rightNeumannCoefficient,

              topBoundaryValue, bottomBoundaryValue,
              leftBoundaryValue, rightBoundaryValue;

  PetscBool   homogeneousDirichlet;

  string      directoryName, objectName;

  MPI_Comm    subCommunicator;

  Vec         localVector;
}DiffusionData;

//Class for simulating 2d diffusion
class diffusionPETSc : public eQ::diffusionSolver
{
//    diffusionPETSc() {}

    PetscErrorCode  ierr;
    MPI_Comm        DIFFU_COMM;
    PetscInt        gridNodesX, gridNodesY, step;
	int				commSize;
    Vec             globalVector;
    DMBoundaryType  boundary;
    DMDAStencilType stencilType;
    DM              distributedArray;
	AO				appOrder;
    KSP             krylovSolver;
    PC              preconditioner;
    DiffusionData   initData, *gridData;
    PetscViewer     printViewer;

    PetscErrorCode  ApplyBoundaryConditions();
    PetscErrorCode  InitializeDiffusion(DiffusionData*, int, char **);
    PetscErrorCode  TimeStep();
    int             KSPIterationCountLastTimeStep();
    PetscErrorCode  WriteGridValues(vector<double>, vector<double>, vector<double>);
    PetscErrorCode  ReadGridValues(vector<double>, vector<double>, vector<double> *);
    PetscErrorCode  RecordData();
    PetscErrorCode  DiffusionFinalize(PetscBool);

	std::vector<double> allXCoordinates, allYCoordinates;

	public:

		diffusionPETSc() {}
		//Vector for storing vector information
		std::vector<double> solution_vector;

//		void initDiffusion(MPI_Comm, std::vector<std::string>, int, char**);
		void initDiffusion(eQ::diffusionSolver::params &initParams);
		void stepDiffusion();
		void setBoundaryValues(const eQ::parametersType &bvals);
		eQ::parametersType getBoundaryFlux(void);
		double getDiffusionConstant(void);
		void writeDiffusionFiles(double timestamp);
		void finalize(void);
};


#endif
