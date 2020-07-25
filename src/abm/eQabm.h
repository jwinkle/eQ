#ifndef EQABM_H
#define EQABM_H

#include "../eQ.h"
#include "../eQcell.h"

#include "../Strain.h"

#include "../inputOutput.h"
#include "Ecoli.h"
#include "cpmEcoli.h"
#include "cpmTrap.h"
#include "cpmHabitat.h"



class eQBaseABM
{
public:

    //virtual destructor needed in base class so that derived classes can be deleted properly:
    virtual ~eQBaseABM()=default;

//    virtual void initCells(int)=0;
    virtual void stepSimulation()=0;
    virtual void updateCellData(double simTime)=0;
    virtual void updateCellModels()=0;
//    virtual void updateCellModels(size_t threads)=0;

    virtual void createDataVectors(eQ::nodeType, eQ::nodeType)=0;

    virtual void writeLookupTable(eQ::nodeType grid, eQ::nodePoint point, eQ::nodeType dof)=0;

};

class eQabm : public eQBaseABM
{
public:
    using HSLdata               = std::vector<double>;
    using HSLgridNodes          = eQ::gridFunction<eQ::nodeType>;

    using HSLgrid               = std::shared_ptr<HSLgridNodes>;
    using HSLsolution           = std::shared_ptr<HSLdata>;
    using channelHSLsolution    = std::shared_ptr<HSLdata>;

    struct Params
	{
        std::shared_ptr<inputOutput>                fileIO;
        std::shared_ptr<eQ::data::files_t>          dataFiles;
        std::shared_ptr<eQ::uniformRandomNumber>    zeroOne;
    };
    Params params;

//    using eQ2DRectangularABM::eQ2DRectangularABM;
    eQabm(const eQabm::Params &);
    ~eQabm() override;

    enum class initType
    {
        RANDOM,
        BANDED,
        NUM_INITTYPES
    };

    void initCells(eQabm::initType howToInit, int numCellsToInit, std::vector<std::shared_ptr<Strain> > &strains);

    class cellContainer
    {
    public:
        virtual ~cellContainer()=default;

        size_t  cellCount(){return _cellCount;}

        void init()
        {
            if(!_cellList.empty())
                std::cout<<"cellList not empty!"<<std::endl;
            clear();
        }
        void clear()
        {
            _cellCount = 0;
            _cellList.clear();
            _strainCounts.clear();
            _eraseCounter=0;
        }
        eQabm::cellContainer &operator<<(std::shared_ptr<Ecoli> cell)
        {//adds cell to list, increment count
            _cellList.push_front(cell);
            ++_cellCount;//must count these manually using a forward_list
            ++_strainCounts[cell->params.strain->getStrainType()];
            return (*this);
        }
        eQabm::cellContainer &operator>>(std::shared_ptr<Ecoli> &cell)
        {//fetches cell only (does not remove it)
            cell = bool(*this) ? (*_icells) : nullptr;
            return (*this);
        }
        eQabm::cellContainer &operator++()
        {//increment the list pointer for next access
            //check if we're at the end of list (in case >> ended)
            if( bool(*this) )
            {//necessary for forward_list to have 2 pointers:
                _icellsPrev = _icells;
                ++_icells;
            }
            return (*this);
        }
        eQabm::cellContainer &operator--()
        {//removes cell from list, decrement count
            --_strainCounts[(*_icells)->params.strain->getStrainType()];
            --_cellCount;
            ++_eraseCounter;
            _icells = _cellList.erase_after(_icellsPrev);
            return (*this);
        }
        operator bool()
        {
            return _icells != _cellList.end();
        }
        void average(double &data)
        {
            if(_cellCount > 0)
                data /= double(_cellCount);
        }
        bool operator>>(double &size)
        {
            size = double(_cellCount);
            return (_cellCount > 0);
        }
        void beginIteration()
        {
            _icells     = _cellList.begin();
            _icellsPrev = _cellList.before_begin();//needed for forward_list C++ data structure
        }
        bool strainFixation()
        {
            bool fixation(true);
            for(strainCount_t & strain : _strainCounts)
            {//first occurence of between (0,1) will clear the flag to 'false'
                fixation &= ((0 == strain.second) || (strain.second == _cellCount));
            }
            return fixation;
        }
        std::vector<double> strainFractions()
        {
            std::vector<double> fracs;
            if(_cellCount > 0)
                for(const strainCount_t strain : _strainCounts)
                        fracs.push_back(double(strain.second)/double(_cellCount));
            return fracs;
        }
        size_t eraseCounter() { return _eraseCounter; }

    protected:
        std::forward_list<std::shared_ptr<Ecoli>>   _cellList;
        size_t                                      _cellCount=0;//must count these manually using a forward_list
        size_t                                      _eraseCounter=0;
        std::map<eQ::Cell::strainType, size_t>      _strainCounts;
        using strainCount_t = std::pair<const eQ::Cell::strainType, size_t>;

        std::forward_list<std::shared_ptr<Ecoli>>::iterator _icells, _icellsPrev;


    };

    cellContainer                                   cellList;
    std::map<std::pair<size_t,size_t>, size_t>      overWrites;

    void printOverWrites()
    {
        for(auto where : overWrites)
        {
            std::cout<<std::endl;
            std::cout<<"overwrite of cell grid pointer: (x,y) = "
                    <<where.second<<"X at: "
                   <<where.first.second<<", "<<where.first.first //(i,j) order => (y,x)
                  <<std::endl;
            std::cout<<std::endl;
        }
        overWrites.clear();
    }



    //2D grid at signaling resolution (2D array of simulation lattice points)
    //dofLookupTable grid is used to lookup the Fenics DOF# from the lattice point location.
    //its size is determined at runtime and it is instantiated after data transfer of the DOF mappings via MPI
    std::vector<HSLgrid>                dofLookupTable;
    std::vector<HSLsolution>            hslSolutionVector;
    std::vector<channelHSLsolution>     channelSolutionVector;
    std::vector<double>                 membraneDiffusionRates;
    std::vector<HSLgrid>                petscLookupTable;
    std::vector<HSLsolution>            petscSolutionVector;

    void createDataVectors(eQ::nodeType globalNodesH, eQ::nodeType globalNodesW) override
    {
        hslSolutionVector.push_back(std::make_shared<HSLdata>(globalNodesH*globalNodesW));
        //this will map x,y position directly to dof
        //CREATE THE LOOKUP TABLE FOR THIS DIFFUSION LAYER:
        dofLookupTable.push_back(std::make_shared<HSLgridNodes>(globalNodesH, globalNodesW));

//        if(bool(eQ::parameters["PETSC_SIMULATION"]))
//        {
//            petscSolutionVector.push_back(std::make_shared<HSLdata>(globalNodesH*globalNodesW));
//            petscLookupTable.push_back(std::make_shared<HSLgridNodes>(globalNodesH, globalNodesW));
//            //populate the grid lookup table here (standard i,j vertex assignment):
//            for (size_t i(0); i < globalNodesH*globalNodesW; ++i)
//            {
//                unsigned jx = i%globalNodesW;
//                unsigned iy = i/globalNodesW;
//                petscLookupTable.back()->grid[iy][jx] = i;
//            }
//        }
    }

    void writeLookupTable(eQ::nodeType grid, eQ::nodePoint point, eQ::nodeType dof) override
    {
        dofLookupTable[grid]->operator[](point) = dof;
    }

    void assignDefaultParameters(Ecoli::Params &cellParams)
    {

        cellParams.space = habitat->get_cpSpace();

        //set default initial cell parameters
        cellParams.mass                        = eQ::Cell::DEFAULT_CELL_MASS;
        cellParams.moment                      = eQ::Cell::DEFAULT_CELL_MOMENT;
        cellParams.width                       = eQ::Cell::DEFAULT_CELL_WIDTH_MICRONS;
        cellParams.meanDivisionLength          = eQ::Cell::DEFAULT_DIVISION_LENGTH_MICRONS;//default, possibly reset below:
        cellParams.divisionLength              = eQ::Cell::DEFAULT_DIVISION_LENGTH_MICRONS;//default, possibly reset below:
        cellParams.doublingPeriodMinutes       = eQ::Cell::DEFAULT_CELL_DOUBLING_PERIOD_MINUTES;
        cellParams.x          = 0.0;
        cellParams.y          = 0.0;
        cellParams.angle      = 0.0;
        cellParams.vx          = 0.0;
        cellParams.vy          = 0.0;
        cellParams.stabilityScaling     = stabilityScaling;

    }

    void stepSimulation()                   override;
    void updateCellData(double simTime)     override;
    void updateCellModels()                 override;
//    void updateCellModels(size_t threads)   override {updateCellModels();}

    void finalizeDataRecording(std::string fpath);


    bool pressureInductionFlag = false;

    double averagePointsPerCell;

    std::shared_ptr<eQ::gridFunction<std::shared_ptr<Ecoli>>>   cellPointers;
    std::shared_ptr<eQ::gridFunction<size_t>>                   gridDataCounter;

    std::shared_ptr<eQ::gridFunction<double>> D11grid;
    std::shared_ptr<eQ::gridFunction<double>> D22grid;
    std::shared_ptr<eQ::gridFunction<double>> D12grid;

    double maxSpringCompression;


    std::vector<std::tuple<double,long,double,long,double>>     divisionList;

private:

    std::shared_ptr<cpmHabitat>                 habitat;
    std::shared_ptr<cpmTrap>                    trap;
    long                                        cellID_factory=1;
    double                                      stabilityScaling;
    bool                                        hslSignaling=false;


    void        recordDivisionEvent(double, std::shared_ptr<Ecoli>, std::shared_ptr<Ecoli>);
    void        updateCells();
    double      rn();


    double dt;
    size_t nodesHighData, nodesWideData;
    size_t nodesHigh, nodesWide;
    double trapWidthMicrons, trapHeightMicrons;
    double nodesPerMicron, nodesPerMicronData;

    size_t nodesToEdge;
};

#endif // EQABM_H
