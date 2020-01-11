#ifndef EQABM_H
#define EQABM_H

#include "../eQ.h"

#include "../Strain.h"

#include "../inputOutput.h"
#include "eColi.h"
#include "cpmEColi.h"
#include "cpmTrap.h"
#include "cpmHabitat.h"


class eQabm
{
public:

    //2D grid at signaling resolution (2D array of simulation lattice points)
    //This grid is used to lookup the Fenics DOF# from the lattice point location.
    //its size is determined at runtime and it is instantiated after data transfer of the DOF mappings via MPI
    typedef  std::shared_ptr<eQ::gridFunction<size_t>>  HSLgrid;
//    typedef  std::shared_ptr<std::vector<double>>       HSLsolution;
    typedef  std::vector<double>                        HSLsolution;
    typedef  std::shared_ptr<std::vector<double>>       channelHSLsolution;

	struct params
	{
		inputOutput		*fileIO;
        eQ::dataFiles_t *dataFiles;

        //declare as vector of pointers to these objects:
        std::vector<eQabm::HSLgrid>                 dofLookupTable;
        std::vector<eQabm::HSLsolution>             hslSolutionVector;
        std::vector<eQabm::channelHSLsolution>      channelSolutionVector;
        std::vector<double>                         membraneDiffusionRates;

        //these should be vectors:
//        eQabm::HSLgrid c4lookup;
//        eQabm::HSLgrid c14lookup;
//        std::vector<double> *c4grid;
//        std::vector<double> *c14grid;
        struct dso_parameters           *pA,*pR;
        size_t seedValue;
    };



    eQabm(const eQabm::params &);
	~eQabm();
    struct params Params;

    void initCells(int);
    void stepChipmunk();
    void updateCellPositions(double simTime);
    void updateCellModels(size_t threads);
    void finalizeDataRecording(std::string fpath);

//    std::vector<std::vector<std::pair<double, double>>> highResolutionDataVector;
    std::vector<std::vector<eColi::highResData_t>> highResolutionDataVector;

    bool timeSeriesDataTrigger = false;
    bool pressureInductionFlag = false;

    size_t cellCount;//must count these manually using a forward_list
    int eraseCounter;
    double strainRatio;

    double averagePointsPerCell;

    std::shared_ptr<eQ::gridFunction<std::shared_ptr<eColi>>> cellPointers;
    std::shared_ptr<eQ::gridFunction<size_t>> gridDataCounter;

    std::shared_ptr<eQ::gridFunction<double>> D11grid;
    std::shared_ptr<eQ::gridFunction<double>> D22grid;
    std::shared_ptr<eQ::gridFunction<double>> D12grid;

    std::vector<double> compressionTimeSeries;
    std::vector<std::vector<size_t>>angleBinTimeSeries;
    std::vector<std::vector<size_t>>angleBinTimeSeries2;
    const size_t numBins = 16;
    std::vector<size_t> binBuffer;
    std::vector<size_t> binBuffer2;
    double maxSpringCompression;

    bool mutantTriggerFlag = false;
    size_t mutantCellNumber;
    double mutant_xpos, mutant_ypos;

    bool aspectRatioInduction = false;

//    //array of unknown rates;  swept over range in Chen
//    double  basalRate[static_cast<unsigned int>(basalRates::numBasalRates)];
//    //strong, med, weak rates table;
//    struct rates  rates;
//    struct dso_parameters pA, pR;

    std::forward_list<std::shared_ptr<eColi>>   cellList;
    std::vector<std::tuple<double,long,double,long,double>> divisionList;

private:
    typedef
        std::forward_list<std::shared_ptr<eColi>>::iterator fli_t;

    std::shared_ptr<cpmHabitat>                 habitat;
    std::shared_ptr<cpmTrap>                    trap;
    long                                        cellID_factory=1;
    double      stabilityScaling;

    std::vector<fli_t> tiBegins, tiEnds;

    void        recordDivisionEvent(double, std::shared_ptr<eColi>, std::shared_ptr<eColi>);
    void        updateCells(fli_t begin, fli_t end);
    void        updateCellGrid(fli_t begin, fli_t end);
    double      rn();

    std::shared_ptr<eQ::uniformRandomNumber> zeroOne;

    double      dt;
    size_t nodesHighData, nodesWideData;
    size_t nodesHigh, nodesWide;
    double trapWidthMicrons, trapHeightMicrons;
    double nodesPerMicron, nodesPerMicronData;

    size_t nodesToEdge;
};

#endif // EQABM_H
