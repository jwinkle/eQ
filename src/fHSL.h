#ifndef FENICS_H
#define FENICS_H
#include <dolfin.h>

#include "eQ.h"
#include "./abm/eColi.h"

#include "../qtensor/pressure.h"
#include "../qtensor/velocity.h"
#include "../qtensor/qtensor.h"
#include "../qtensor/gradv.h"
#include "../qtensor/data.h"
#include "../qtensor/vdata.h"

#include "../fenics/hsl.h"
#include "../fenics/hslD.h"
#include "../fenics/hslRobin.h"
#include "../fenics/boundary.h"
#include "../fenics/AdvectionDiffusion.h"
#include "../fenics/AD1Dss.h"


#include <random>
#include <cmath>
#include <fstream>
#include <csignal>
#include <iostream>
#include <chrono>
#include <string>


#include "Expressions.h"




using namespace dolfin;
#include <mshr.h>
using namespace mshr;


/*
class fenicsParent
{
public:
    enum  {U0, DT, F, D, tD11, tD22, tD12, NUM_FENICSVARS};
    //Generic tuple of Function space, BiLinear form and Liner form:
    typedef std::tuple<std::shared_ptr<dolfin::FunctionSpace>,
                        std::shared_ptr<dolfin::Form>,
                        std::shared_ptr<dolfin::Form>>
            fenicsSpaceForms;
    //Generic vector of variables (constants and expressions):
    typedef std::vector<std::shared_ptr<dolfin::Variable>>
            fenicsVariables;

    //Create the fenicsSpaceForms tuple by passing in the mesh:
    virtual fenicsSpaceForms makeSpaceForms(std::shared_ptr<Mesh> mesh) =0;
//    virtual fenicsSpaceForms makeSpaceForms(std::shared_ptr<Mesh> mesh);
    //Set the form variables with vector input:
    virtual void setFormParameters(fenicsParent::fenicsVariables &, double) =0;
    virtual void setRobinBoundaryValue(double) =0;
    virtual void setRobinBoundaryRate(double) =0;
};


//Templated derived class from fenicsParent class
template <class FS, class LF, class BF>
class fenicsForms : public fenicsParent
{
public:
    fenicsSpaceForms makeSpaceForms(std::shared_ptr<Mesh> mesh)
    {
        functionSpace   = std::make_shared<FS>(mesh);
        linearForm      = std::make_shared<LF>(functionSpace);
        bilinearForm    = std::make_shared<BF>(functionSpace,functionSpace);
        return std::make_tuple(functionSpace, linearForm, bilinearForm);
    }
    ~fenicsForms()
    {
        bilinearForm.reset();
        linearForm.reset();
        functionSpace.reset();
    }
    void setRobinBoundaryValue(double value);
    void setRobinBoundaryRate(double rate);
    void setFormParameters(fenicsVariables &, double);
    std::shared_ptr<FS> functionSpace;
    std::shared_ptr<LF> linearForm;
    std::shared_ptr<BF> bilinearForm;
};
//for the general template, empty functions for setting boundary value, rate
template <class FS, class LF, class BF>
void fenicsForms<FS, LF, BF>::setRobinBoundaryValue(double value) {}//set the boundary well value:
template <class FS, class LF, class BF>
void fenicsForms<FS, LF, BF>::setRobinBoundaryRate(double rate) {}//set the boundary well value:

template <> inline  //<> used for specific types, as below
void fenicsForms<hsl::FunctionSpace, hsl::LinearForm, hsl::BilinearForm>::setFormParameters(fenicsVariables &fvars, double r)
{//ONLY INITIALIZE WHAT THE FORM ACTUALLY HAS:

    linearForm->u0 =  std::dynamic_pointer_cast<dolfin::Function>(fvars.at(fenicsForms::U0));
    linearForm->dt =  std::dynamic_pointer_cast<dolfin::Constant>(fvars.at(fenicsForms::DT));
    linearForm->f =  std::dynamic_pointer_cast<dolfin::Constant>(fvars.at(fenicsForms::F));
    bilinearForm->D =  std::dynamic_pointer_cast<dolfin::Constant>(fvars.at(fenicsForms::D));
    bilinearForm->dt =  std::dynamic_pointer_cast<dolfin::Constant>(fvars.at(fenicsForms::DT));
}
template <> inline  //<> used for specific types, as below
//void fenicsForms<hslD::FunctionSpace, hslD::LinearForm, hslD::BilinearForm>::setFormParameters(fenicsVariables &fvars, double r)
void fenicsForms<hslD::FunctionSpace, hslD::Form_L, hslD::Form_a>::setFormParameters(fenicsVariables &fvars, double r)
{//ONLY INITIALIZE WHAT THE FORM ACTUALLY HAS:
    linearForm->u0 =  std::dynamic_pointer_cast<dolfin::Function>(fvars.at(fenicsForms::U0));
    linearForm->dt =  std::dynamic_pointer_cast<dolfin::Constant>(fvars.at(fenicsForms::DT));
    linearForm->f =  std::dynamic_pointer_cast<dolfin::Constant>(fvars.at(fenicsForms::F));
    bilinearForm->D =  std::dynamic_pointer_cast<dolfin::Constant>(fvars.at(fenicsForms::D));
    bilinearForm->D11 =  std::dynamic_pointer_cast<dolfin::Expression>(fvars.at(fenicsForms::tD11));
    bilinearForm->D22 =  std::dynamic_pointer_cast<dolfin::Expression>(fvars.at(fenicsForms::tD22));
    bilinearForm->D12 =  std::dynamic_pointer_cast<dolfin::Expression>(fvars.at(fenicsForms::tD12));
    bilinearForm->dt =  std::dynamic_pointer_cast<dolfin::Constant>(fvars.at(fenicsForms::DT));
}


//  ROBIN BOUNDARY CONDITION TEMPLATE FORMS:

//for the Robin template, set the storage value (s) and the rate (r)
template <> inline  //<> used for specific types, as below
void fenicsForms<hslRobin::FunctionSpace, hslRobin::LinearForm, hslRobin::BilinearForm>::setRobinBoundaryValue(double value)
    {linearForm->s = std::make_shared<Constant>(value);}//set the boundary well value:
template <> inline  //<> used for specific types, as below
void fenicsForms<hslRobin::FunctionSpace, hslRobin::LinearForm, hslRobin::BilinearForm>::setRobinBoundaryRate(double rate)
    {linearForm->r = std::make_shared<Constant>(rate);}//set the boundary well value:
template <> inline  //<> used for specific types, as below
void fenicsForms<hslRobin::FunctionSpace, hslRobin::LinearForm, hslRobin::BilinearForm>::setFormParameters
                (fenicsVariables &fvars, double rate)
{
    linearForm->u0 =  std::dynamic_pointer_cast<dolfin::Function>(fvars.at(fenicsForms::U0));
    linearForm->dt =  std::dynamic_pointer_cast<dolfin::Constant>(fvars.at(fenicsForms::DT));
    linearForm->f =  std::dynamic_pointer_cast<dolfin::Constant>(fvars.at(fenicsForms::F));
    bilinearForm->D =  std::dynamic_pointer_cast<dolfin::Constant>(fvars.at(fenicsForms::D));
    bilinearForm->dt =  std::dynamic_pointer_cast<dolfin::Constant>(fvars.at(fenicsForms::DT));

    auto zero = std::make_shared<dolfin::Constant>(0.0);
    auto robinRate = std::make_shared<dolfin::Constant>(rate);
    linearForm->r = robinRate;
    bilinearForm->r = robinRate;
    linearForm->s = zero;
}
*/

class fenicsVariable
{//superset of all possible constants and expression classess for the diffusion solvers
public:
    struct data
    {
        std::shared_ptr<dolfin::Constant> zero;
        std::shared_ptr<dolfin::Constant> one;
        std::shared_ptr<dolfin::Constant> dt;
        std::shared_ptr<dolfin::Constant> D;
        std::shared_ptr<AnisotropicDiffusionTensor> tensorD11;
        std::shared_ptr<AnisotropicDiffusionTensor> tensorD22;
        std::shared_ptr<AnisotropicDiffusionTensor> tensorD12;

        std::shared_ptr<dolfin::Constant> v;//constant velocity

        std::shared_ptr<dolfin::Constant> s_left;//Robin boundary value
        std::shared_ptr<dolfin::Constant> r_left;//Robin "rate" (effective units=length^-1...1/distance from trap boundary to Robin boundary value location)
        std::shared_ptr<dolfin::Constant> s_right;//Robin boundary value
        std::shared_ptr<dolfin::Constant> r_right;//Robin "rate" (effective units=length^-1...1/distance from trap boundary to Robin boundary value location)

        std::shared_ptr<MeshFunction<size_t>> meshFunction;
        std::shared_ptr<MeshFunction<size_t>> meshFunctionChannel;

    };
};

class fenicsShell
{
public:
    typedef  std::shared_ptr<eQ::gridFunction<size_t>>  HSLgrid;


    std::shared_ptr<Mesh>                       mesh;
    std::shared_ptr<dolfin::Function>           u;
    std::shared_ptr<dolfin::Function>           u0;
    std::vector<std::shared_ptr<const dolfin::DirichletBC>> dbc;

    std::shared_ptr<dolfin::LinearVariationalProblem>   LVP;
    std::shared_ptr<dolfin::LinearVariationalSolver>    LVS;

    std::shared_ptr<dolfin::FunctionSpace>              V;

    //after derived classes create spaces and functions and populate forms, generate the problem+solver objects:
    void createLinearVariationalSolver()
    {
        LVP = std::make_shared<dolfin::LinearVariationalProblem>(_a, _L, u, dbc);
        LVS = std::make_shared<dolfin::LinearVariationalSolver>(LVP);

        //sizes here are local and include ghost nodes:
        dof_from_vertex 	= dolfin::vertex_to_dof_map(*V);
        vertex_from_dof 	= dolfin::dof_to_vertex_map(*V);
        mesh_coords 		= mesh->coordinates();
        numNodesLocal 		= mesh_coords.size()/2;//numNodes is for this process only
    }

    void createGridCoordinatesToDofMapping()
    {
        auto globalNodesH = size_t(eQ::parameters["simulationTrapHeightMicrons"])*size_t(eQ::parameters["nodesPerMicronSignaling"]) + 1;
        auto globalNodesW = size_t(eQ::parameters["simulationTrapWidthMicrons"])*size_t(eQ::parameters["nodesPerMicronSignaling"]) + 1;

        //CREATE THE LOOKUP TABLE FOR THIS DIFFUSION LAYER:
        dofLookupTable = std::make_shared<eQ::gridFunction<size_t>>(globalNodesH, globalNodesW);

        //we step through every other entry to get x,y values iteratively (total size 2*globalNodes size)
        for (size_t i(0); i < numNodesLocal; i++)
        {//mesh_coords is a 1D array of (x,y) pairs in vertex order (not in dof order):
            //convert to a 1D array of std::pair
            auto x = mesh_coords[2*i];
            auto y = mesh_coords[2*i+1];

            unsigned jx = unsigned(round(x * double(eQ::parameters["nodesPerMicronSignaling"])));
            unsigned iy = unsigned(round(y * double(eQ::parameters["nodesPerMicronSignaling"])));
            //the point of this all is to have this mapping from physical node # to the fenics dof #,
            //allows to look that up the dof directly for read/write in the main acquisition loop
//                dof_from_grid[grid]->grid[iy][jx] = gridDofs[grid][i];//note row,column = y,x

            //WRITE THE DOF LOOKUP TABLE:
            //write directly to the ABM parameters data structure:
            dofLookupTable->grid[iy][jx] = size_t(dof_from_vertex[i]);
        }

    }

    void createCoordinatesToDofMappingChannel()
    {
        //dimensions are simulation scale:
        auto uleft = double(0.0);
        auto uright = double(0.0);
//        auto uleft      = double(eQ::parameters["channelLengthMicronsLeft"]);
//        auto uright     = double(eQ::parameters["channelLengthMicronsRight"]);
        auto utrap      = double(eQ::parameters["simulationTrapWidthMicrons"]);

        nodesForChannels = size_t(
                    ceil( (utrap + uleft + uright) * double(eQ::parameters["nodesPerMicronSignaling"]) ));
        nodesForChannels++;//always one more vertex for each dimension vs. # elements

        zerovecChannel = std::vector<double>(nodesForChannels);

        //CREATE THE LOOKUP TABLE FOR THIS DIFFUSION LAYER:
        dofLookupTableChannel = std::vector<size_t>(nodesForChannels);

        for (size_t i(0); i < nodesForChannels; i++)
        {//for an interval mesh, mesh_coords is a 1D array of x locations in vertex order (not in dof order):
            //reset x position to left-most edge is x=0 for coverting to an index value from 0:
            auto x = mesh_coords[i] + uleft;

            size_t jx = size_t(round(x * double(eQ::parameters["nodesPerMicronSignaling"])));

            //WRITE THE DOF LOOKUP TABLE:
            //write directly to the ABM parameters data structure:
            dofLookupTableChannel[jx] = size_t(dof_from_vertex[i]);
        }
        //create a vector of dof's for the boundary interface to the trap:
        size_t trapLeft = size_t(round(uleft * double(eQ::parameters["nodesPerMicronSignaling"])));
        size_t trapRight = size_t(round((uleft + utrap) * double(eQ::parameters["nodesPerMicronSignaling"])));
        boundaryDofChannel = std::vector<size_t>(trapRight-trapLeft+1);

        size_t j=0;
        for(size_t i(trapLeft); i<trapRight+1; ++i)
        {
            boundaryDofChannel[j++] = (dofLookupTableChannel[i]);
        }

        std::cout<<"dofLookupTableChannel, size = "<<dofLookupTableChannel.size()<<std::endl;
//        for(auto dof : dofLookupTableChannel)
//            std::cout<<dof<<", ";
//        std::cout<<std::endl;

    }

    std::shared_ptr<dolfin::Form>           _L, _a;
    size_t                                  numNodesLocal, nodesForChannels;
    std::vector<double>                     mesh_coords;
    std::vector<size_t>                     vertex_from_dof;//dof_to_vertex_map(const FunctionSpace& space);
    std::vector<int>                        dof_from_vertex;//vertex_to_dof_map(const FunctionSpace& space);
    HSLgrid                                 dofLookupTable;
    std::vector<size_t>                     dofLookupTableChannel;
    std::vector<size_t>                     boundaryDofChannel;
    std::vector<double>                     zerovecChannel;

};

//Base class for working with Fenics .ufl files (PDE "forms") for linear problems:
template <class FS, class LF, class BF>//FunctionSpace, LinearForm, BiLinearForm
class fenicsBaseClass : public fenicsShell
{
public:
    //note: if defined outside class definition, add this to the signature:
    //template <class FS, class LF, class BF>
    fenicsBaseClass(){}
    fenicsBaseClass(std::shared_ptr<fenicsShell> &shell)
    {
        //populate the base class pointer with the derived class object called (allows access to base variables)
        shell = std::shared_ptr<fenicsShell>(this);
    }
    void createSpaceFormsFunctions()
    {
        V = std::make_shared<FS>(mesh);
        L = std::make_shared<LF>(V);
        a = std::make_shared<BF>(V, V);

        u   = std::make_shared<dolfin::Function>(V);    u->set_allow_extrapolation(true);
        u0  = std::make_shared<dolfin::Function>(V);    u0->set_allow_extrapolation(true);

        //copy to the shell's base class pointers
        _L = L;
        _a = a;
    }

    //necessary for each .ufl file implementation to define this for specific template class types:
    virtual void setFormParameters(fenicsVariable::data &data);

    std::shared_ptr<LF>                       L;
    std::shared_ptr<BF>                       a;

};
template <class FS, class LF, class BF>
void fenicsBaseClass<FS,LF,BF>::setFormParameters(fenicsVariable::data &data)
{
    std::cout<<"ERROR: fenicsBaseClass<FS,LF,BF>::setFormParameters NOT IMPLEMENTED FOR THIS TEMPLATE!"<<std::endl;
}
template<> inline
void fenicsBaseClass<hslD::FunctionSpace, hslD::LinearForm, hslD::BilinearForm>::setFormParameters(fenicsVariable::data &data)
{
    L->u0   =  u0;
    L->dt   =  data.dt;
    L->f    =  data.zero;
    L->r1    =  data.r_left;
    L->s1    =  data.s_left;
    L->r2    =  data.r_right;
    L->s2    =  data.s_right;

    a->D    =  data.D;
    a->D11  =  data.tensorD11;
    a->D22  =  data.tensorD22;
    a->D12  =  data.tensorD12;
    a->dt   =  data.dt;
    a->r1    =  data.r_left;
    a->r2    =  data.r_right;

    L->ds =     data.meshFunction;
    a->ds =     data.meshFunction;
}
template<> inline
void fenicsBaseClass<
AdvectionDiffusion::FunctionSpace, AdvectionDiffusion::LinearForm, AdvectionDiffusion::BilinearForm>
    ::setFormParameters(fenicsVariable::data &data)
{
    L->D = data.D;
    L->dt = data.dt;
    L->v  = data.v;
    L->u0 = u0;

    a->D = data.D;
    a->dt = data.dt;
    a->v  = data.v;

    L->r1    =  data.r_left;
    L->s1    =  data.s_left;
    L->r2    =  data.r_right;
    L->s2    =  data.s_right;
    a->r1    =  data.r_left;
    a->r2    =  data.r_right;

    L->ds =     data.meshFunctionChannel;
    a->ds =     data.meshFunctionChannel;

}

template<> inline
void fenicsBaseClass<AD1Dss::FunctionSpace, AD1Dss::LinearForm, AD1Dss::BilinearForm>::setFormParameters(fenicsVariable::data &data)
{
//    L->D = data.D;
//    L->dt = data.dt;
//    L->v  = data.v;
//    L->u0 = u0;
    L->f = data.zero;

    a->D = data.D;
//    a->dt = data.dt;
    a->v  = data.v;

}


class fenicsInterface : public eQ::diffusionSolver
{
public:


    virtual ~fenicsInterface();
    fenicsInterface(void);
    fenicsInterface(const eQ::diffusionSolver::params &);

    eQ::diffusionSolver::params                     myParams;

    //this pointer is returned by the member base class constructor as a parameter;
    //used for base class methods without having to know the template pointer type (avoids messy templating in this class)
    std::shared_ptr<fenicsShell>                   shell;

    //declare all possible pointers (avoids using a messy template class for the interface)
    std::shared_ptr<fenicsBaseClass
        <hslD::FunctionSpace, hslD::Form_L, hslD::Form_a>>                   hslD;
    std::shared_ptr<fenicsBaseClass
        <hslRobin::FunctionSpace, hslRobin::Form_L, hslRobin::Form_a>>       hslRobin;
    std::shared_ptr<fenicsBaseClass
        <hsl::FunctionSpace, hsl::Form_L, hsl::Form_a>>                      hsl;

    void createSpaceFormsFunctions()
    {
        std::cout<<"POINTER STATUS IS:"
                <<"\n\t hslD = "<<hslD
               <<"\n\t hslRobin = "<<hslRobin
              <<"\n\t hsl = "<<hsl
             <<std::endl;
        //wrapper method to init the form that is created (depends on only one being created!)
        if(nullptr != hslD) hslD->createSpaceFormsFunctions();
        else if(nullptr != hslRobin) hslRobin->createSpaceFormsFunctions();
        else if(nullptr != hsl) hsl->createSpaceFormsFunctions();
    }
    void setFormParameters(fenicsVariable::data &data)
    {
        //wrapper method to init the form that is created (depends on only one being created!)
        if(nullptr != hslD) hslD->setFormParameters(data);
        else if(nullptr != hslRobin) hslRobin->setFormParameters(data);
        else if(nullptr != hsl) hsl->setFormParameters(data);
    }
    void computeBoundaryFlux();

    std::shared_ptr<fenicsBaseClass
        <AdvectionDiffusion::FunctionSpace, AdvectionDiffusion::Form_L, AdvectionDiffusion::Form_a>>       topChannel;
    std::shared_ptr<fenicsBaseClass
        <AdvectionDiffusion::FunctionSpace, AdvectionDiffusion::Form_L, AdvectionDiffusion::Form_a>>       bottomChannel;


    bool                        isFenicsRootNode;
    size_t                      nodesChannel;
    double                      wellScaling;



    //DIFFUSION SOLVER TEMPLATE METHODS:
//    void initDiffusion(MPI_Comm comm, std::vector<std::string> filePaths, int argc, char* argv[]);
    void initDiffusion(size_t id, MPI_Comm comm, std::string filePath, double D, double dt, int argc, char* argv[]);
    void initDiffusion(eQ::diffusionSolver::params &);
    void stepDiffusion();
    void writeDiffusionFiles(double dt);
    void finalize(void);
    void setBoundaryValues(const eQ::parametersType &bvals);
    eQ::parametersType getBoundaryFlux(void);

    void updateChannels();
    void writeDataFiles(double dt);


    int myRankMPI, numFenicsPEs;
    size_t nodesH, nodesW;

    std::vector<std::shared_ptr<dolfin::File>> outputFiles;

    std::vector<
            std::pair<
                    std::shared_ptr<dolfin::File>,
                    std::shared_ptr<dolfin::Function>
                    >>
                hslWriter;


    std::vector<std::pair<double, double>>  coords;

    std::vector<double> solution_vector;
    std::vector<double> solution_vectorTopChannel;
    std::vector<double> solution_vectorBottomChannel;
    std::vector<double> fluxTopChannel;
    std::vector<double> fluxBottomChannel;

    //FUNCTIONS:
    std::shared_ptr<dolfin::Function> uvec;
    std::shared_ptr<dolfin::Function> u2;
    std::shared_ptr<dolfin::Scalar> boundaryFlux;

    std::shared_ptr<dolfin::Function> uCT, uCB;
    std::shared_ptr<dolfin::Function> u0CT, u0CB;


    std::shared_ptr<updatingDirchletBoundary> boundaryChannelLeft;
    std::shared_ptr<updatingDirchletBoundary> boundaryChannelRight;
    std::shared_ptr<updatingDirchletBoundary> boundaryCompartment;

    std::shared_ptr<std::vector<double>> D11, D22, D12;


protected:
//private:
    bool isDataRecordingNode;
    void fenicsClassInit();
    void initHSLFiles();

    void createMesh(MPI_Comm);
    void createHSL();
    void mapGridNodes();

    double myDiffusionConstant;
    double totalBoundaryFlux;

    std::shared_ptr<data::FunctionSpace> VscalarData;
    std::shared_ptr<vdata::FunctionSpace> VvectorData;

    //BOUNDARY CONDITIONS:
    //the boundary resolution in simulation units:
    double      h;

    std::vector<std::shared_ptr<const dolfin::DirichletBC>> dbce;
    std::shared_ptr<dolfin::DirichletBC> 	walls;




    std::shared_ptr<boundary::Functional> V0;
    std::shared_ptr<boundary::CoefficientSpace_u> VU;
    std::shared_ptr<Function> ub;
    std::shared_ptr<dolfin::Assembler> flux;
};



////////////////////////////////////////////////////////////////////////////////
//                          BOUNDARY CONDITIONS
////////////////////////////////////////////////////////////////////////////////
class DirichletBoundary_openWalls : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
        {return on_boundary;}
  };
class DirichletBoundary_threeWalls : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
      return on_boundary
      and (near(x[1], 0));  //bottom wall open => dbc is near y=0 (and on boundary)
    }
    //            or
    //            near(x[1], gridHeight)
  };
class DirichletBoundary_twoWalls : public SubDomain
{
public:
    DirichletBoundary_twoWalls(double gridHeight, double gridWidth)
        : trapH(gridHeight), trapW(gridWidth) {}
    bool inside(const Array<double>& x, bool on_boundary) const
    {
      return on_boundary
              and (
                  near(x[1], 0)
                or
                  near(x[1], trapH)
                  );  //bottom wall open => dbc is near y=0 (and on boundary)
    }
    double trapH, trapW;
};
class DirichletBoundary_oneWall: public SubDomain
{
public:
    DirichletBoundary_oneWall(double gridHeight, double gridWidth)
        : trapH(gridHeight), trapW(gridWidth) {}
    bool inside(const Array<double>& x, bool on_boundary) const
    {
      return on_boundary
                and (
                  near(x[0], 0)
                or
                  near(x[0], trapW)
                or
                  near(x[1], 0)
              );  //bottom wall open => dbc is near y=0 (and on boundary)
    }
    double trapH, trapW;
};
class DirichletBoundary_leftWall: public SubDomain
{
public:
    DirichletBoundary_leftWall(double gridHeight, double gridWidth)
        : trapH(gridHeight), trapW(gridWidth) {}
    bool inside(const Array<double>& x, bool on_boundary) const
    {
      return on_boundary
                and (
                  near(x[1], trapH)
                or
                  near(x[0], trapW)
                or
                  near(x[1], 0)
              );  //bottom wall open => dbc is near y=0 (and on boundary)
    }
    double trapH, trapW;
};
class LEGI_boundary: public SubDomain
{
public:
    LEGI_boundary() {}
    bool inside(const Array<double>& x, bool on_boundary) const
    {
      return on_boundary;
    }
};

class DirichletBoundary_Edges1D : public SubDomain
{
public:
    DirichletBoundary_Edges1D(double left, double right)
        : leftEdge(left), rightEdge(right) {}
    virtual bool inside(const Array<double>& x, bool on_boundary) const
    {
      return on_boundary
          and (
              near(x[0], leftEdge)
           or
              near(x[0], rightEdge)
          );
    }
    double leftEdge, rightEdge;
};

class DirichletBoundary_Edges1DL : public DirichletBoundary_Edges1D
{
public:
    virtual bool inside(const Array<double>& x, bool on_boundary) const
    {
      return on_boundary
          and (
              near(x[0], leftEdge)
//           or
//              near(x[0], rightEdge)
          );
    }
};
class DirichletBoundary_Edges1DR : public DirichletBoundary_Edges1D
{
public:
    virtual bool inside(const Array<double>& x, bool on_boundary) const
    {
      return on_boundary
          and (
//              near(x[0], leftEdge)
//           or
              near(x[0], rightEdge)
          );
    }
};


class DirichletBoundary_TrapEdges: public SubDomain
{
public:
    enum class edge
    {
        TOP, BOTTOM, LEFT, RIGHT,
        NUM_EDGES
    };

    DirichletBoundary_TrapEdges(double gridHeight, double gridWidth, edge which)
        : trapH(gridHeight), trapW(gridWidth), whichEdge(which) {}
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        bool nearThisEdge = false;
        if(edge::TOP == whichEdge) nearThisEdge = near(x[1], trapH);
        else if(edge::BOTTOM == whichEdge) nearThisEdge = near(x[1], 0.0);
        else if(edge::LEFT == whichEdge) nearThisEdge = near(x[0], 0.0);
        else if(edge::RIGHT == whichEdge) nearThisEdge = near(x[0], trapW);

      return (on_boundary && nearThisEdge);
    }
    double trapH, trapW;
    edge whichEdge;
};


#endif // FENICS_H
