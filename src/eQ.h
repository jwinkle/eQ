#ifndef EQ_H
#define EQ_H

#include <iostream>
#include "mpi.h"


//https://github.com/nlohmann/json
#include "../nlohmann/json.hpp"
// for convenience
using json = nlohmann::json;

namespace eQ {

using nodeType = size_t;

//=======================================================================================
enum class HSLType
{
    C4HSL,
    C14HSL,
    NUM_HSLTYPES
};
enum class dataParameterType
{
    HSL,
    PROTEIN,
    VECTOR,

    NUM_DATAPARAMETERTYPES
};
//=======================================================================================
template<class T>
class gridFunction
{
public:
    gridFunction(size_t nh, size_t nw)
                : nodesHigh(nh), nodesWide(nw)
    {
        for(size_t i(0); i<nodesHigh; i++)
            grid.push_back(std::vector<T>(nodesWide, T(0)));
    }
    void clear()
    {
        if(grid.empty())
        {
            std::cout<<"grid empty error!"<<std::endl;
            return;
        }
        for(size_t i(0); i<nodesHigh; i++)
            for(size_t j(0); j<nodesWide; j++)
                grid[i][j] = T(0);
    }
    void assign(T value)
    {
        if(grid.empty())
        {
            std::cout<<"grid empty error!"<<std::endl;
            return;
        }
        for(size_t i(0); i<nodesHigh; i++)
            for(size_t j(0); j<nodesWide; j++)
                grid[i][j] = T(value);
    }
    bool isValidIndex(std::pair<size_t, size_t> point)
    {
        return (point.first < nodesHigh) && (point.second < nodesWide);
    }
    void linearize2Dgrid(std::vector<T> &dataArray)
    {
        auto nrows = grid.size();
        auto ncols = grid[0].size();
        size_t index=0;
        for(size_t i(0); i<nrows; i++)
            for(size_t j(0); j<ncols; j++)
                dataArray.at(index++) = grid[i][j];
    }
    std::vector<std::vector<T>> grid;
private:
    size_t nodesHigh, nodesWide;
};
//=======================================================================================

//=======================================================================================
//=======================================================================================
class data
{
public:
    class tensor;
    struct record
    {
        eQ::dataParameterType                   type;
        size_t                                  index;
        std::string                             fileName;
        std::shared_ptr<eQ::data::tensor>       data;
    };
    using parametersType    = json;
    using stringsType       = std::string;
    using files_t           = std::vector<data::record>;

    static bool							isControllerNode;
    static bool							initializedParameters;
    static eQ::data::parametersType     parameters;

    //HSL type=key, {D,d} (diffusion constants in media D and membrane rate d)
    //note: initialized in main
    static std::map<std::string, std::vector<double>> physicalDiffusionRates;

    //=======================================================================================
        static inline std::pair<double, double> xy_from_ij(size_t i, size_t j, double n)
        {//Note: j is horizontal (x) direction, i is vertical (y)
            return std::make_pair(double(j)/n,double(i)/n);
        }
        static inline std::pair<size_t, size_t> ij_from_xy(double x, double y, double n)
        {//Note: j is horizontal (x) direction, i is vertical (y)
            size_t j = size_t(round(x * n));
            size_t i = size_t(round(y * n));
            return std::make_pair(i,j);
        }
        static inline std::pair<size_t, size_t> ij_from_xy(std::pair<double, double> xy, double n)
        {//Note: j is horizontal (x) direction, i is vertical (y)
            size_t j = size_t(round(xy.first * n));
            size_t i = size_t(round(xy.second * n));
            return std::make_pair(i,j);
        }
        static inline size_t index_from_ij(size_t i, size_t j, size_t nh, size_t nw)
        {//Note: j is horizontal (x) direction, i is vertical (y)
            return ((i < nh) && (j < nw)) ? j+i*nw : 0;
        }
    //=======================================================================================

    class tensor
        {
        public:
            enum rank
            {
                SCALAR = 0,
                VECTOR = 1,
                NUM_RANKS
            };
            tensor(size_t n, size_t y, size_t x, size_t thisRank)
                : nodesPerMicron(n),//use the Data version of per-micron
                  rank(thisRank)//scalar data source default
            {
                auto nh = y*n + 1;
                auto nw = x*n + 1;
                for(size_t i(0); i<=rank; ++i)//NOTE: rank=0 is scalar, loop once (uses <=)
                {
                    grid.push_back(std::make_shared<eQ::gridFunction<double>>(nh,nw));
                }
            }

            tensor(const tensor &)
                {std::cout<<"tensorDataSource copyConstructor not defined!\n";}

            size_t getRank(){return rank;}
            double eval(double x, double y)
            {
                auto ij = ij_from_xy(x,y,nodesPerMicron);
                return grid[0]->grid[ij.first][ij.second];
            }
            std::pair<double,double> evalVector(double x, double y)
            {
                auto ij = ij_from_xy(x,y,nodesPerMicron);
                return std::make_pair(
                    grid[0]->grid[ij.first][ij.second],
                    grid[1]->grid[ij.first][ij.second]);
            }
            std::vector<std::shared_ptr<eQ::gridFunction<double>>> grid;

        private:
            double nodesPerMicron;
            size_t rank;
        };
};

//=======================================================================================
//=======================================================================================
class simulationTiming
    {
    public:
        simulationTiming()=default;
        void dt(double dt)
        {
            _dt = dt;
            stepsPerMin = size_t(round(1.0/_dt));
            stepsPerHour = 60*stepsPerMin;
            _timer = 0.0;
            _timeSteps = 0;
            flags.clear();
        }
        void    setSimulationTimeMinutes(size_t mins)   { _simulationTimeMinutes = mins; }
        size_t  getSimulationTimeMinutes()              { return _simulationTimeMinutes; }
        double  dt()            { return _dt; }
        double  getTime()       { return MPI_Wtime(); }
        double  tare()          { _timeTare = getTime(); return _timeTare; }
        double  timer()         { return getTime() - _timeTare; }
        double  simTime()       { return _timer; }
        size_t  steps()         { return _timeSteps; }
        bool    stepTimer()
        {
            ++_timeSteps;
            _timer += _dt;
            for(auto &flag : flags)
            {
                flag.second.thrown = (_timer >= flag.second.when);
            }
            return (_timer <= _simulationTimeMinutes);
        }

        bool    periodicTimeMinutes(size_t mins)
        {
            return (_timeSteps%(mins*stepsPerMin) == 0);
        }
        void    setTimerFlag(std::string key, double when)
        {
            flags.insert({key, {when, false, false}});
        }
        bool    flagThrown(std::string key)
        {
            return (flags.at(key).thrown && !flags.at(key).ignore);
        }
        void    flagIgnore(std::string key)
        {
            flags.at(key).ignore = true;
        }
//        operator bool() { return _timer < _simulationTimeMinutes; }

        size_t stepsPerMin, stepsPerHour;

    private:
        struct alarm { double when; bool thrown; bool ignore; };
        double  _dt;
        size_t  _simulationTimeMinutes;
        size_t  _timeSteps=0;
        double  _timeTare;
        double  _timer;
        std::map<std::string, alarm> flags;
    };
//=======================================================================================
class diffusionSolver
{
public:
    struct params
    {
        int argc;
        char** argv;
        size_t                                      uniqueID;
        MPI_Comm                                    comm;
        double                                      dt;
        double                                      D_HSL;
        std::string                                 filePath;
        std::shared_ptr<eQ::data::files_t>            dataFiles;
        double                                      trapHeightMicrons;
        double                                      trapWidthMicrons;
        double                                      nodesPerMicron;
        double                                      trapChannelVelocity;
    };
//        virtual void initDiffusion(MPI_Comm comm, std::vector<std::string> filePaths, int argc, char* argv[]) =0;// {}
//        virtual void initDiffusion(size_t id, MPI_Comm comm, std::string filePath, double D, double dt, int argc, char* argv[]) =0;// {}
    virtual void initDiffusion(eQ::diffusionSolver::params &) =0;// {}
    virtual void stepDiffusion() {}
    virtual void setBoundaryValues(const eQ::data::parametersType &bvals) =0;
    virtual eQ::data::parametersType getBoundaryFlux(void) =0; //must at least define value of "totalFlux"
    virtual void writeDiffusionFiles(double timestamp) =0;
    virtual void finalize(void) {}
};

class boundaryCondition
{
public:
    enum class wall
    {
        LEFT, RIGHT, TOP, BOTTOM,
        ALL,
        NUM_WALLS
    };
    enum class type
    {
        DIRICHLET, NEUMANN, ROBIN,
        CHANNELUPDATE,
        DIRICHLET_0, NEUMANN_0,
        NUM_TYPES
    };

    boundaryCondition & operator<<(boundaryCondition::wall wall)
    {
        _wall = wall;
        if(type::CHANNELUPDATE == _type)
        {
            if( (wall::TOP == _wall) || (wall::BOTTOM == _wall) )
                writeDirichlet(_wall, -1.0);
        }
        else if(type::DIRICHLET_0 == _type)
        {
            if(wall::ALL == wall)
            {
                writeDirichlet(wall::TOP, 0.0);
                writeDirichlet(wall::BOTTOM, 0.0);
                writeDirichlet(wall::LEFT, 0.0);
                writeDirichlet(wall::RIGHT, 0.0);
            }
            else
                writeDirichlet(_wall, 0.0);
        }
        else if(type::NEUMANN_0 == _type)
            writeNeumann(_wall, 0.0);
        else if(type::DIRICHLET == _type)
            writeDirichlet(_wall, 0.0);
        else if(type::NEUMANN == _type)
            writeNeumann(_wall, 0.0);
        else if(type::ROBIN == _type)
            writeRobin(_wall, 1.0);

        return (*this);
    }

    boundaryCondition & operator<<(boundaryCondition::type type)
    {
        _type = type;
        return (*this);
    }

protected:
    eQ::data::parametersType      _bcs;
    boundaryCondition::wall _wall;
    boundaryCondition::type _type;

private:
    int     _robinState =0;
    double  value1      =0.0;

    std::vector<std::string> walls = {"left", "right", "top", "bottom"};
    std::vector<std::string> typess = {"Dirichlet", "Neumann", "Robin"};

    void writeDirichlet(wall w, double v)
    {
        size_t type = static_cast<size_t>(type::DIRICHLET);
        size_t wall = static_cast<size_t>(w);
        _bcs[walls[wall].c_str()] =
            {typess[type].c_str() , {0.0, 1.0, v}};
    }
    void writeNeumann(wall w, double v)
    {
        size_t type = static_cast<size_t>(type::NEUMANN);
        size_t wall = static_cast<size_t>(w);
        _bcs[walls[wall].c_str()] =
            {typess[type].c_str() , {1.0, 0.0, v}};
    }
    void writeRobin(wall w, double v)
    {
        size_t type = static_cast<size_t>(type::ROBIN);
        size_t wall = static_cast<size_t>(w);
        _bcs[walls[wall].c_str()] =
            {typess[type].c_str() , {v, 1.0, 0.0}};
    }

friend eQ::data::parametersType &operator<<(eQ::data::parametersType &os, const eQ::boundaryCondition &bc);
};


}//end namespace eQ::




#endif // EQ_H
