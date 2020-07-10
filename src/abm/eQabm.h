#ifndef EQABM_H
#define EQABM_H

#include "../eQ.h"
#include "../eQcell.h"

#include "../Strain.h"

#include "../inputOutput.h"
#include "eColi.h"
#include "cpmEColi.h"
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
    virtual void updateCellModels(size_t threads)=0;

    virtual void createDataVectors(eQ::nodeType, eQ::nodeType)=0;
    virtual void writeLookupTable(eQ::nodeType grid, eQ::nodeType iy, eQ::nodeType jx, eQ::nodeType dof)=0;

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
        std::shared_ptr<eQ::data::files_t>            dataFiles;
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

    std::vector<double> initCells(eQabm::initType howToInit, int numCellsToInit, std::vector<std::shared_ptr<Strain> > &strains);

    class cellContainer
    {
    public:
        cellContainer()=default;
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
        }
        eQabm::cellContainer &operator<<(std::shared_ptr<eColi> cell)
        {
            _cellList.push_front(cell);
            ++_cellCount;//must count these manually using a forward_list
            return (*this);
        }

        eQabm::cellContainer &operator>>(std::shared_ptr<eColi> &cell)
        {
            cell = bool(*this) ? (*_icells) : nullptr;
            return (*this);
        }
        eQabm::cellContainer &operator>>(size_t &size)
        {
            size = _cellCount;
            return (*this);
        }
        eQabm::cellContainer &operator++()
        {//necessary for forward_list to have 2 pointers:
            _icellsPrev = _icells;
            ++_icells;
            return (*this);
        }
        eQabm::cellContainer &operator--()
        {
            _icells = _cellList.erase_after(_icellsPrev);
            --_cellCount;
            return (*this);
        }
        operator bool()
        {
            return _icells != _cellList.end();
        }
        void beginIteration()
        {
            _icells     = _cellList.begin();
            _icellsPrev = _cellList.before_begin();//needed for forward_list C++ data structure
        }
    protected:
        std::forward_list<std::shared_ptr<eColi>>   _cellList;
        size_t                                      _cellCount=0;//must count these manually using a forward_list

        std::forward_list<std::shared_ptr<eColi>>::iterator _icells, _icellsPrev;

    };

    cellContainer                                   cellList;




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

    void writeLookupTable(eQ::nodeType grid, eQ::nodeType iy, eQ::nodeType jx, eQ::nodeType dof) override
    {
        dofLookupTable[grid]->grid[iy][jx] = dof;
    }

    void assignDefaultParameters(eColi::Params &cellParams)
    {

        cellParams.space = habitat->get_cpSpace();

        //set default initial cell parameters
        cellParams.baseData.mass                        = eQ::Cell::DEFAULT_CELL_MASS;
        cellParams.baseData.moment                      = eQ::Cell::DEFAULT_CELL_MOMENT;
        cellParams.baseData.width                       = eQ::Cell::DEFAULT_CELL_WIDTH_MICRONS;
        cellParams.baseData.meanDivisionLength          = eQ::Cell::DEFAULT_DIVISION_LENGTH_MICRONS;//default, possibly reset below:
        cellParams.baseData.divisionLength              = eQ::Cell::DEFAULT_DIVISION_LENGTH_MICRONS;//default, possibly reset below:
        cellParams.baseData.doublingPeriodMinutes       = eQ::Cell::DEFAULT_CELL_DOUBLING_PERIOD_MINUTES;
        cellParams.baseData.x          = 0.0;
        cellParams.baseData.y          = 0.0;
        cellParams.baseData.angle      = 0.0;
        cellParams.baseData.vx          = 0.0;
        cellParams.baseData.vy          = 0.0;
        cellParams.stabilityScaling     = stabilityScaling;

    }

//    void initCells(int) override;
    void stepSimulation() override;
    void updateCellData(double simTime) override;
    void updateCellModels() override;
    void updateCellModels(size_t threads)override {updateCellModels();}

    void finalizeDataRecording(std::string fpath);


    bool pressureInductionFlag = false;

    int eraseCounter;
    double strainRatio;

    double averagePointsPerCell;

    std::shared_ptr<eQ::gridFunction<std::shared_ptr<eColi>>> cellPointers;
    std::shared_ptr<eQ::gridFunction<size_t>> gridDataCounter;

    std::shared_ptr<eQ::gridFunction<double>> D11grid;
    std::shared_ptr<eQ::gridFunction<double>> D22grid;
    std::shared_ptr<eQ::gridFunction<double>> D12grid;

    double maxSpringCompression;


    bool aspectRatioInduction = false;


    //the list data structure that holds the cell objects:
//    std::forward_list<std::shared_ptr<eColi>>               cellList;
    std::vector<std::tuple<double,long,double,long,double>> divisionList;

private:
    using fli_t = std::forward_list<std::shared_ptr<eColi>>::iterator;

    std::shared_ptr<cpmHabitat>                 habitat;
    std::shared_ptr<cpmTrap>                    trap;
    long                                        cellID_factory=1;
    double                                      stabilityScaling;
    bool                                        hslSignaling=false;

    std::vector<fli_t> tiBegins, tiEnds;

    void        recordDivisionEvent(double, std::shared_ptr<eColi>, std::shared_ptr<eColi>);
    void        updateCells();
    double      rn();


    double      dt;
    size_t nodesHighData, nodesWideData;
    size_t nodesHigh, nodesWide;
    double trapWidthMicrons, trapHeightMicrons;
    double nodesPerMicron, nodesPerMicronData;

    size_t nodesToEdge;
};

#endif // EQABM_H
