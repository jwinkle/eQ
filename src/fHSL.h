#ifndef FENICS_H
#define FENICS_H
#include <dolfin.h>

#include "eQ.h"
#include "./abm/Ecoli.h"

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

#include <petscsys.h>

#include "Expressions.h"

//using namespace dolfin;
//#include <mshr.h>
//using namespace mshr;

struct fenicsVariables
{//superset of all possible constants and expression classess for the diffusion solvers
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

class fenicsShell
{
public:
    typedef  std::shared_ptr<eQ::gridFunction<size_t>>  HSLgrid;

    virtual ~fenicsShell()=default;
    virtual void createSpaceFormsFunctions()    {}
    virtual void setFormParameters()            {}
    virtual void updateSolver()                 {}
    virtual void createLinearProblem()          {}
    virtual void createLinearVariationalSolver(){}

    std::shared_ptr<Mesh>                       mesh;
    std::shared_ptr<dolfin::Function>           u;
    std::shared_ptr<dolfin::Function>           u0;
    std::vector<std::shared_ptr<const dolfin::DirichletBC>> dbc;

    std::shared_ptr<dolfin::LinearVariationalProblem>   LVP;
    std::shared_ptr<dolfin::LinearVariationalSolver>    LVS;

    std::shared_ptr<dolfin::FunctionSpace>              V;

    void getVertexMappings()
    {
        //sizes here are local and include ghost nodes:
        dof_from_vertex 	= dolfin::vertex_to_dof_map(*V);
        vertex_from_dof 	= dolfin::dof_to_vertex_map(*V);
        mesh_coords 		= mesh->coordinates();
    }
    virtual void createGridCoordinatesToDofMapping()
    {
        getVertexMappings();

        auto globalNodesH = size_t(eQ::data::parameters["simulationTrapHeightMicrons"])*size_t(eQ::data::parameters["nodesPerMicronSignaling"]) + 1;
        auto globalNodesW = size_t(eQ::data::parameters["simulationTrapWidthMicrons"])*size_t(eQ::data::parameters["nodesPerMicronSignaling"]) + 1;

        //CREATE THE LOOKUP TABLE FOR THIS DIFFUSION LAYER:
        dofLookupTable = std::make_shared<eQ::gridFunction<size_t>>(globalNodesH, globalNodesW);

        //default grid is 2D rectangular:
        //we step through every other entry to get x,y values iteratively (total size 2*globalNodes size)
        for (size_t i(0); i < mesh_coords.size()/2; i++)
        {//mesh_coords is a 1D array of (x,y) pairs in vertex order (not in dof order):
            //convert to a 1D array of std::pair
            auto x = mesh_coords[2*i];
            auto y = mesh_coords[2*i+1];

            size_t jx = size_t(round(x * double(eQ::data::parameters["nodesPerMicronSignaling"])));
            size_t iy = size_t(round(y * double(eQ::data::parameters["nodesPerMicronSignaling"])));
            //the point of this all is to have this mapping from physical node # to the fenics dof #,
            //allows to look that up the dof directly for read/write in the main acquisition loop
            auto point = std::pair<size_t, size_t>{iy,jx};
            dofLookupTable->operator[](point) = size_t(dof_from_vertex[i]);
        }
    }

    std::vector<double>                     mesh_coords;
    std::vector<size_t>                     vertex_from_dof;//dof_to_vertex_map(const FunctionSpace& space);
    std::vector<int>                        dof_from_vertex;//vertex_to_dof_map(const FunctionSpace& space);
    HSLgrid                                 dofLookupTable;
    fenicsVariables                         data;

};


class fenicsData : public fenicsShell
{//for future expansion:
public:
};


//Base class for working with Fenics .ufl files (PDE "forms") for linear problems:
template <class FS, class LF, class BF>//FunctionSpace, LinearForm, BiLinearForm
class fenicsBaseClass : public fenicsShell
{
public:
    //note: if defined outside class definition, add this to the signature:
    //template <class FS, class LF, class BF>
    //necessary for each .ufl file implementation to define this for specific template class types:
    void setFormParameters() override;
    void updateSolver() override;
    void createSpaceFormsFunctions() override
    {
        V = std::make_shared<FS>(mesh);
        L = std::make_shared<LF>(V);
        a = std::make_shared<BF>(V, V);

        u   = std::make_shared<dolfin::Function>(V);    u->set_allow_extrapolation(true);
        u0  = std::make_shared<dolfin::Function>(V);    u0->set_allow_extrapolation(true);
    }
    void createLinearProblem() override
    {
        setFormParameters();
        createLinearVariationalSolver();
        createGridCoordinatesToDofMapping();
    }
    void createLinearVariationalSolver() override
    {
        LVP = std::make_shared<dolfin::LinearVariationalProblem>(a, L, u, dbc);
        LVS = std::make_shared<dolfin::LinearVariationalSolver>(LVP);
    }
    //we need template versions of these to set ufl form parameters for each unique .ufl file:
    std::shared_ptr<LF>                       L;
    std::shared_ptr<BF>                       a;
};

template <class FS, class LF, class BF>
void fenicsBaseClass<FS,LF,BF>::
setFormParameters()
{
    std::cout<<"ERROR: fenicsBaseClass<FS,LF,BF>::setFormParameters NOT IMPLEMENTED FOR THIS TEMPLATE!"<<std::endl;
}
template <class FS, class LF, class BF>
void fenicsBaseClass<FS,LF,BF>::
updateSolver()
{
    std::cout<<"ERROR: fenicsBaseClass<FS,LF,BF>::updateSolver NOT IMPLEMENTED FOR THIS TEMPLATE!"<<std::endl;
}


template<> inline
void fenicsBaseClass<hslD::FunctionSpace, hslD::LinearForm, hslD::BilinearForm>::
updateSolver()
{
    L->f    =  data.zero;
    L->r1    =  data.r_left;
    L->r2    =  data.r_right;

    a->D    =  data.D;
    a->r1    =  data.r_left;
    a->r2    =  data.r_right;
}

template<> inline
void fenicsBaseClass<hslD::FunctionSpace, hslD::LinearForm, hslD::BilinearForm>::
setFormParameters()
{
    std::cout<<"called template<> inline setFormParameters() override"<<std::endl;

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
void fenicsBaseClass<AD1Dss::FunctionSpace, AD1Dss::LinearForm, AD1Dss::BilinearForm>::
setFormParameters()
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

using baseAdvDiff = fenicsBaseClass<
                AdvectionDiffusion::FunctionSpace,
                AdvectionDiffusion::Form_L,
                AdvectionDiffusion::Form_a>;

class fenicsChannel : public baseAdvDiff
{
public:
    fenicsChannel(fenicsVariables &dataRef) : cdata(dataRef) {} //tie channel data dirctly to the trap data
    void updateSolver() override
    {
        L->v  = cdata.v;
        a->v  = cdata.v;
        L->r1    =  cdata.r_left;
        L->r2    =  cdata.r_right;
        a->r1    =  cdata.r_left;
        a->r2    =  cdata.r_right;
    }
    void setFormParameters() override
    {
        std::cout<<"called channel setFormParameters() override"<<std::endl;
        L->D = cdata.D;
        L->dt = cdata.dt;
        L->v  = cdata.v;
        L->u0 = u0;

        a->D = cdata.D;
        a->dt = cdata.dt;
        a->v  = cdata.v;

        L->r1    =  cdata.r_left;
        L->s1    =  cdata.s_left;
        L->r2    =  cdata.r_right;
        L->s2    =  cdata.s_right;
        a->r1    =  cdata.r_left;
        a->r2    =  cdata.r_right;

        L->ds =     cdata.meshFunctionChannel;
        a->ds =     cdata.meshFunctionChannel;

    }

    void createGridCoordinatesToDofMapping() override
    {
        getVertexMappings();

        //dimensions are simulation scale:
        auto uleft = double(0.0);
        auto uright = double(0.0);
//        auto uleft      = double(eQ::data::parameters["simulationChannelLengthLeft"]);
//        auto uright     = double(eQ::data::parameters["simulationChannelLengthRight"]);
        auto utrap      = double(eQ::data::parameters["simulationTrapWidthMicrons"]);

        nodesForChannels = size_t(
                    ceil( (utrap + uleft + uright) * double(eQ::data::parameters["nodesPerMicronSignaling"]) ));
        nodesForChannels++;//always one more vertex for each dimension vs. # elements

        zerovecChannel = std::vector<double>(nodesForChannels);

        //CREATE THE LOOKUP TABLE FOR THIS DIFFUSION LAYER:
        dofLookupTableChannel = std::vector<size_t>(nodesForChannels);

        for (size_t i(0); i < nodesForChannels; i++)
        {//for an interval mesh, mesh_coords is a 1D array of x locations in vertex order (not in dof order):
            //reset x position to left-most edge is x=0 for coverting to an index value from 0:
            auto x = mesh_coords[i] + uleft;

            size_t jx = size_t(round(x * double(eQ::data::parameters["nodesPerMicronSignaling"])));

            //WRITE THE DOF LOOKUP TABLE:
            //write directly to the ABM parameters data structure:
            dofLookupTableChannel[jx] = size_t(dof_from_vertex[i]);
        }
        //create a vector of dof's for the boundary interface to the trap:
        size_t trapLeft = size_t(round(uleft * double(eQ::data::parameters["nodesPerMicronSignaling"])));
        size_t trapRight = size_t(round((uleft + utrap) * double(eQ::data::parameters["nodesPerMicronSignaling"])));
        boundaryDofChannel = std::vector<size_t>(trapRight-trapLeft+1);

        size_t j=0;
        for(size_t i(trapLeft); i<trapRight+1; ++i)
        {
            boundaryDofChannel[j++] = (dofLookupTableChannel[i]);
        }

//        std::cout<<"boundaryDofChannel, size = "<<boundaryDofChannel.size()<<std::endl;
//        for(auto dof : boundaryDofChannel)
//            std::cout<<dof<<", ";
//        std::cout<<std::endl;

    }

    std::vector<size_t>                     dofLookupTableChannel;
    std::vector<size_t>                     boundaryDofChannel;
    std::vector<double>                     zerovecChannel;
    size_t                                  nodesForChannels;
private:
    fenicsVariables             &cdata;
};
class fenicsInterface : public eQ::diffusionSolver
{
public:
    ~fenicsInterface();
    fenicsInterface(void);
    fenicsInterface(const eQ::diffusionSolver::params &);

    eQ::diffusionSolver::params                     myParams;

    //this pointer is returned by the member base class constructor as a parameter;
    //used for base class methods without having to know the template pointer type (avoids messy templating in this class)
    std::shared_ptr<fenicsShell>                   shell;

    void computeBoundaryFlux();

    std::shared_ptr<fenicsChannel>      topChannel;
    std::shared_ptr<fenicsChannel>      bottomChannel;


    bool                        isFenicsRootNode;
    size_t                      nodesChannel;
    double                      wellScaling;

    //DIFFUSION SOLVER TEMPLATE METHODS:
    void initDiffusion(size_t id, MPI_Comm comm, std::string filePath, double D, double dt, int argc, char* argv[]);
    void initDiffusion(eQ::diffusionSolver::params &) override;
    void stepDiffusion() override;
    void writeDiffusionFiles(double dt) override;
    void finalize(void) override;
    void setBoundaryValues(const double);
    void setRobinBoundaryConditions();

    void updateRobinParameters()
    {
        topChannel->updateSolver();
        topChannel->createLinearVariationalSolver();
            bottomChannel->updateSolver();
            bottomChannel->createLinearVariationalSolver();
            if(false == lateralNeumann)
            {
                shell->updateSolver();
                shell->createLinearVariationalSolver();
            }
    }

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
    std::vector<double> solution_vectorModified;
    std::vector<double> solution_vectorTopChannel;
    std::vector<double> solution_vectorBottomChannel;
    std::vector<double> topChannelData;
    std::vector<double> bottomChannelData;
    std::vector<double> fluxTopChannel;
    std::vector<double> fluxBottomChannel;

    double totalBoundaryFlux;

    //FUNCTIONS:
    std::shared_ptr<dolfin::Function> uvec;
    std::shared_ptr<dolfin::Function> u2;
    std::shared_ptr<dolfin::Scalar> boundaryFlux;

    std::shared_ptr<dolfin::Function> uCT, uCB;
    std::shared_ptr<dolfin::Function> u0CT, u0CB;


//    std::shared_ptr<updatingDirchletBoundary> boundaryChannelLeft;
//    std::shared_ptr<updatingDirchletBoundary> boundaryChannelRight;
    std::shared_ptr<updatingDirchletBoundary> boundaryCompartment;

    std::shared_ptr<std::vector<double>> D11, D22, D12;

private:
    bool isDataRecordingNode;
    void fenicsClassInit();
    void initHSLFiles();

    void createMesh(MPI_Comm);
    void createHSL();
    void mapGridNodes();

    double myDiffusionConstant;

    std::shared_ptr<data::FunctionSpace> VscalarData;
    std::shared_ptr<vdata::FunctionSpace> VvectorData;

    double rightRate = 0.0;
    double leftRate = 0.0;

    bool    lateralNeumann = false;

    double channelFlowVelocity = 0.0;

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


class DirichletBoundary_TrapEdge: public SubDomain
{
public:
    enum edge
    {
        TOP, BOTTOM, LEFT, RIGHT,
        NUM_EDGES
    };

    DirichletBoundary_TrapEdge(double gridHeight, double gridWidth, edge which)
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
