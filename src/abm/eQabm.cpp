#include "eQabm.h"

#include <cstdlib>
#include <unistd.h>

//NOTE:  ABM will have two spatial resolutions:
//one for full sub-micron position to map HSL
//and one for general data recording (1 micron or greater typically).
//the two resolutions are passed in simData structure as:
//nodesPerMicron (full resolution) and nodesPerMicronData (reduced for data recording)
//nodesHigh and nodesWide use full resolution for making the pointer array grid

double eQabm::rn()
{//interface for uniform random numbers in [0,1]
    return params.zeroOne->randomNumber();
}

eQabm::eQabm(const eQabm::Params &p)
    : eQBaseABM(), params(p)
{
    cellList.init();

    dt              = double(eQ::data::parameters["dt"]);

    trapHeightMicrons   = double(eQ::data::parameters["simulationTrapHeightMicrons"]);
    trapWidthMicrons    = double(eQ::data::parameters["simulationTrapWidthMicrons"]);
    nodesPerMicron      = double(eQ::data::parameters["nodesPerMicronSignaling"]);
    nodesPerMicronData  = double(eQ::data::parameters["nodesPerMicronData"]);

    nodesHigh       = size_t(trapHeightMicrons * nodesPerMicron) + 1;
    nodesWide       = size_t(trapWidthMicrons * nodesPerMicron) + 1;
    nodesHighData   = size_t(trapHeightMicrons * nodesPerMicronData) + 1;
    nodesWideData   = size_t(trapWidthMicrons * nodesPerMicronData) + 1;

    std::cout<<"Creating eQabm habitat and trap..."<<std::endl;

    //create the habitat:
    habitat = std::make_shared<cpmHabitat>(dt);

    //keep old params type for now:
    cpmTrap::Params trapParams;
        trapParams.space = habitat->get_cpSpace();
    //create the trap:
    trap = std::make_shared<cpmTrap>(trapParams);
    habitat->setTrap(trap);//don't omit this

    std::cout<<"Created eQabm habitat and trap..."<<std::endl;
    std::cout<<"eQ::w,h:"<<trapWidthMicrons<<", "<<trapHeightMicrons<<std::endl;
    stabilityScaling = dt/0.001;//needed for stability of the ABM growth, per model (see Phys.Bio paper appendix)
    std::cout<<"stabilityScaling: "<<stabilityScaling<<std::endl;
    std::cout<<"Number of data files: "<<(*params.dataFiles).size()<<std::endl;
    std::cout<<std::endl;

    gridDataCounter =
                std::make_shared<eQ::gridFunction<size_t>>(nodesHighData, nodesWideData);

    std::cout<<"creating cellPointers: "<<nodesHigh<<", "<<nodesWide<<std::endl;
    //Note: uses full resolution of grid (not the data recording resolution)
    cellPointers =
            std::make_shared<eQ::gridFunction<std::shared_ptr<Ecoli>>>(nodesHigh, nodesWide);

    auto numHSLGrids = size_t(eQ::data::parameters["D_HSL"].size());//size of diffusion vector is # of grids (known globally)
    hslSignaling = (numHSLGrids > 0);
    membraneDiffusionRates = std::vector<double>(eQ::data::parameters["membraneDiffusionRates"].get<std::vector<double>>());

    //grids for defining anisotropic diffusion tensor
    D11grid = std::make_shared<
            eQ::gridFunction<double>>(nodesHigh, nodesWide);
    D22grid = std::make_shared<
            eQ::gridFunction<double>>(nodesHigh, nodesWide);
    D12grid = std::make_shared<
            eQ::gridFunction<double>>(nodesHigh, nodesWide);

    //compute the max. # of nodes from a pole center to the boundary of the cell (for searching for interior points)
    nodesToEdge = size_t(round(nodesPerMicron * eQ::Cell::DEFAULT_CELL_WIDTH_MICRONS/2.0));

    maxSpringCompression = 0.0;
    averagePointsPerCell = 0.0;

}

void
eQabm::initCells(eQabm::initType howToInit, int numCellsToInit, std::vector<std::shared_ptr<Strain>> &strains)
{

    Ecoli::Params       cellParams;
    assignDefaultParameters(cellParams);

    auto numStrains = strains.size();

    double              rx,ry,ra;
    size_t              whichStrain = 0;

    cellParams.meanDivisionLength =
            double(eQ::data::parameters["defaultAspectRatioFactor"]) * eQ::Cell::DEFAULT_DIVISION_LENGTH_MICRONS;
    cellParams.divisionCorrelationAlpha =
            double(eQ::data::parameters["divisionCorrelationAlpha"]);


    overWrites.clear();

    auto getRandomCell = [&]()
    {
        rx = rn() * trapWidthMicrons;
        ry = rn() * trapHeightMicrons;
        ra = rn() * 2*M_PI;

        cellParams.x = rx;
        cellParams.y = ry;
        cellParams.angle = ra;

        //add some distribution over the intial cells' division length
//        double dn = double(eQ::data::parameters["divisionNoiseScale"]);
//        double noiseFactor = (rn()-0.5) * dn;
//        cellParams.divisionLength  = (1.0 + noiseFactor) * cellParams.meanDivisionLength;
        cellParams.divisionLength  = cellParams.meanDivisionLength;
        cellParams.length          = (1.0 + rn()) * 0.5 * cellParams.divisionLength;

        if(eQabm::initType::RANDOM == howToInit)
        {
            //randomly selects a strain
            whichStrain = size_t(floor(rn() * numStrains));
        }
        if(eQabm::initType::BANDED == howToInit)
        {
            //selects strain based on x location:
            whichStrain = size_t(floor(rx * double(numStrains)/trapWidthMicrons));
        }

        //clone the strain object from the objects passed in:
        cellParams.strain = strains[whichStrain]->clone();

        return std::make_shared<Ecoli>(cellParams);

    };

    for(int i(0); i<numCellsToInit; ++i)
    {
        //CREATE CELL OBJECT:
        auto thisCell = getRandomCell();
        cellList << thisCell;

        std::cout<<"cell: ("
                    <<thisCell->getCenter_x()<<","<<thisCell->getCenter_y()<<") "<<thisCell->getAngle()<<"=angle; strain = "
                    <<thisCell->strain->getStrainType()<<std::endl;

        thisCell->setCellID(cellID_factory++);
        thisCell->setParentID(-(thisCell->getCellID()));//initial cells set parent to negative cellID.
        recordDivisionEvent(0.0, nullptr, thisCell);

    }

    std::cout<<"Initial strain fractions: ";
    for(auto frac : cellList.strainFractions())
    {
        std::cout<<"  "<<frac;
    }
    std::cout<<std::endl;

}
void eQabm::stepSimulation()
{
    habitat->stepSimulation(dt);
}
void eQabm::updateCellData(double simTime)
{//this is done while waiting for HSL data transfer (single-threaded)

    //UPDATES PHYSICS MODEL RESULTS (cell length, position, angle, spring compression, velocity, ...)
    //REMOVES CELLS OUT-OF-BOUNDS
    //DIVIDES CELLS THAT REACH THEIR DIVISION LENGTH TARGET
    //UPDATES LIST DATA STRUCTURE, REMOVING (EXIT) OR ADDING (DIVISION) CELLS AS NEEDED

    std::vector<std::shared_ptr<Ecoli>> newDaughterCells;
    std::shared_ptr<Ecoli> cell;


    //TRAVERSE LIST OF CELLS: (uses C++ STL forward list data structure)
    cellList.beginIteration();
    while(cellList >> cell)//delays cellList increment until after possible removal of cell
    {
        //==============================================
        //  UPDATE PHYSICS MODEL AND CHECK FOR DIVISION:
        //==============================================
        //CELL MODEL:
        std::shared_ptr<Ecoli> daughter = cell->updatePhysicsModel(*params.zeroOne);//pass random number stack by reference
        if(nullptr != daughter)
        {//divided cell, check if outside trap (but record division to parent first):
            daughter->setCellID(cellID_factory++);
            daughter->setParentID(cell->getCellID());
            recordDivisionEvent(simTime, cell, daughter);
            habitat->outsideTrap(daughter->cpmCell)
                ? daughter.reset() : newDaughterCells.push_back(daughter);
        }
        //==============================================
        //  CHECK FOR OUTSIDE OF BOUNDARY:
        //==============================================
        (habitat->outsideTrap(cell->cpmCell))
                ? --cellList : ++cellList;
    }
    for(auto &newCell : newDaughterCells)
    {//add daughter cells to main cell list:
        cellList << newCell;
    }

    //new positions and daughter cells are set; update x,y grid for pointers:
    cellPointers->clear();
    eQ::nodePoint index;
    double springEnergy_L2 = 0.0;


    cellList.beginIteration();
    while( ++(cellList >> cell) )
    {
        index = eQ::data::ij_from_xy(cell->getCenter_x(),cell->getCenter_y(),nodesPerMicron);

        if(nullptr != (*cellPointers)[index])
        {
            auto where = eQ::data::ij_from_xy(cell->getCenter_x(),cell->getCenter_y(),nodesPerMicronData);
            ++overWrites[where];
        }

        (*cellPointers)[index] = cell;

        springEnergy_L2 += cell->cpmCell->compression * cell->cpmCell->compression;
        if(cell->cpmCell->compression > maxSpringCompression) maxSpringCompression = cell->cpmCell->compression;

    }
}


void eQabm::updateCellModels()
{//launches multi-threaded cell model updates (removed for now)
    //SINGLE-THREADED CELL UPDATE PATH:
    updateCells();
}
void eQabm::updateCells()
{//cell update function called with distributed partition of the cell list:
 //typically called multi-threaded with list range and no risk of race conditions for cells
    double x,y;

    eQ::nodePoint dataPoint, gridPoint;

    //ANISOTROPIC DIFFUSION CODE:
    double Dx = double(eQ::data::parameters["AnisotropicDiffusion_Axial"]);
    double Dy = double(eQ::data::parameters["AnisotropicDiffusion_Transverse"]);
    //init the grid to isotropic diffusion: (these grids copy directly to the .ufl file tensor)
    //note:  these are 1D vectors using standard translation for (i,j) entries using the eQ helper functions
    D11grid->assign(1.0);//diagonal scaling
    D22grid->assign(1.0);//diagonal scaling
    D12grid->assign(0.0);//symmetric, off-diagonal scaling

    averagePointsPerCell = 0.0;

    std::shared_ptr<Ecoli> thisCell;
    cellList.beginIteration();
    while( ++(cellList >> thisCell) )
    {
        double cellLength = thisCell->getLengthMicrons();
        x = thisCell->getCenter_x();
        y = thisCell->getCenter_y();

        gridPoint = eQ::data::ij_from_xy(x,y,nodesPerMicron);//uses full scaling
        dataPoint = eQ::data::ij_from_xy(x,y,nodesPerMicronData);//uses data scaling
        //keep tally of "hits" at this grid point for averaging:        
        ++(gridDataCounter->operator[](dataPoint));

        std::vector<eQ::nodePoint>  cellPoints;//vector of (i,j) pairs
        std::vector<double>         hslData;

        auto findInteriorPoints = [&](double npm, std::vector<eQ::nodePoint> &cellPoints)
        {
            //implicitly passed by reference: x, y, thisCell, nodesToEdge, nodesHigh, nodesWide

            //get location of pole centers:
            //NOTE: need to check if they are out of bounds! (only cell centers are checked for cell removal)
            auto poleA = eQ::data::ij_from_xy(thisCell->polePositionA, npm);//uses full grid scaling
            auto poleB = eQ::data::ij_from_xy(thisCell->polePositionB, npm);//uses full grid scaling

            //determine the search area for inside-cell points;  start with min, max i,j values:
            size_t i1, i2, j1, j2;
            if(poleA.first > poleB.first)   {i1 = poleB.first; i2 = poleA.first;}
            else                            {i1 = poleA.first;  i2 = poleB.first;}
            if(poleA.second > poleB.second) {j1 = poleB.second; j2 = poleA.second;}
            else                            {j1 = poleA.second;  j2 = poleB.second;}
            //now add a 1/2 cell width number of nodes to widen area to entire cell interior (range checking also):
            i1 = (i1 >= nodesToEdge) ? (i1 - nodesToEdge) : 0;
            j1 = (j1 >= nodesToEdge) ? (j1 - nodesToEdge) : 0;
            i2 = ((i2 + nodesToEdge) >= (nodesHigh-1)) ? nodesHigh-1 : (i2 + nodesToEdge);
            j2 = ((j2 + nodesToEdge) >= (nodesWide-1)) ? nodesWide-1 : (j2 + nodesToEdge);

            //now search the rectangular area computed above and call pointIsInCell() method on each point:
            for(size_t pointi(i1); pointi<=i2; pointi++)
            {
                for(size_t pointj(j1); pointj<=j2; pointj++)
                {
                    auto xypoint = eQ::data::xy_from_ij(pointi, pointj, npm);
                    if(thisCell->cpmCell->pointIsInCell(xypoint))
                        cellPoints.push_back(std::make_pair(pointi, pointj));
                }
            }
            if(cellPoints.empty())
            {
//                std::cout<<"empty cell points"<<std::endl;
                cellPoints.push_back(eQ::data::ij_from_xy(thisCell->getCenter_x(), thisCell->getCenter_y(), npm));
            }

        };
        auto setDiffusionTensor = [&](double theta, std::vector<std::pair<size_t, size_t>> &cellPoints)
        {
            //compute rotation matrix entries:
            //Dtensor = Rtheta * [D11 D12; D12 D22] * R-theta
            //  = Dx*cos^2theta + Dy*sin^2theta   (Dx-Dy)*sintheta*costheta;
            //    (Dx-Dy)*sincost                 Dx*sin^2theta + Dy*cos^2theta;
            double cos2t = cos(theta)*cos(theta);
            double sin2t = sin(theta)*sin(theta);
            double sincost = sin(theta)*cos(theta);

            for(auto &point : cellPoints)
            {
                if(D11grid->isValidIndex(point))
                {
                    D11grid->operator[](point) = Dx*cos2t + Dy*sin2t;
                    D22grid->operator[](point) = Dx*sin2t + Dy*cos2t;
                    D12grid->operator[](point) = (Dx-Dy)*sincost;
                }
            }
        };
        auto readHSL = [&](HSLsolution HSLvector, HSLgrid HSLlookup, std::vector<std::pair<size_t, size_t>> &cellPoints)
        {
            //read from each interior point of the cell (computed above) and equal-average values read:
            double HSL = 0.0;
            for(auto &point : cellPoints)
            {
                //fenics implementation:
                auto location = HSLlookup->operator[](point);//translates the i,j position to the fenics DOF entry in the solution vector
                HSL += HSLvector->at(location);
            }
            return HSL/double(cellPoints.size());
        };
        auto writeHSL = [&](double HSL, HSLsolution HSLvector, HSLgrid HSLlookup, std::vector<std::pair<size_t, size_t>> &cellPoints)
        {//NOTE:  HSL passed in is in units of concentration [nM]...must scale to extra-cellular volume to ensure convservation of molecule #
            //we need the fraction of volume outside of cell relative to total rectangular volume 1um wide x 1um tall x cell length (um)
            //[HSL] is the delta concentraion of HSL for the cell; need to determine the amount to put on one grid point:
            //[HSL]*cellVolume=#HSL; to update 1um^2 2D grid point with #HSL molecules and xum^2 extra-cellular volume per grid point:
            //H1, new concentration at grid point;  H0, previous.
            //H1 = [(H0*xum^2) + #HSL]/xum^2 = H0 + #HSL/xum^2, where #HSL = [HSL]*cellVolume

            //note:  use simulation units throughout since length scaling cancels:
            double numberHSL = eQ::Cell::nanoMolarToMoleculeNumber(HSL, cellLength);
            double updatePerSquareMicron = numberHSL/eQ::Cell::computeExtraCellularVolumeFraction(cellLength);
            double updateForOneGridPoint = updatePerSquareMicron * nodesPerMicron * nodesPerMicron;

            //distribute over interior points of the cell:
            double dHSL = updateForOneGridPoint/double(cellPoints.size());
            for(auto &point : cellPoints)
            {
                //fenics implementation:
                auto location = HSLlookup->operator[](point);//translates the i,j position to the fenics DOF entry in the solution vector
                HSLvector->at(location) += dHSL;
            }
        };
        //here, we populate the data structure (fenics recording node) with the cell data
        auto recordCelldata = [&]()
        {
            //std::shared_ptr<eQ::dataFiles_t> = = std::vector<gridData>
            for(auto &dataFile : *(params.dataFiles))
            {
                if(eQ::dataParameterType::HSL == dataFile.type)
                {
                    auto which = static_cast<Strain::hsl>(dataFile.index);
                    dataFile.data->grid[0]->operator[](dataPoint)
                            = eQ::Cell::moleculeNumberToNanoMolar(thisCell->strain->getHSL(which), cellLength);
                }
                if(eQ::dataParameterType::PROTEIN == dataFile.type)
                {
                    auto which = static_cast<Strain::concentrations>(dataFile.index);
                    dataFile.data->grid[0]->operator[](dataPoint)
                            = eQ::Cell::moleculeNumberToNanoMolar(thisCell->strain->getProteinNumber(which), cellLength);
                }
                if(eQ::dataParameterType::VECTOR == dataFile.type)
                {

                }
            }//end for
        };//end lambda recordCelldata
        auto updateCellModel = [&]()
        {
            hslData.clear();
            //HSL READ:
            for(size_t i(0); i< hslSolutionVector.size(); i++)
            {
//                        if(bool(eQ::data::parameters["PETSC_SIMULATION"]))
//                            hslData.push_back(readHSL(petscSolutionVector[i], petscLookupTable[i], cellPoints));
//                        else
                    hslData.push_back(readHSL(hslSolutionVector[i], dofLookupTable[i], cellPoints));
            }

            auto deltaHSL = thisCell->strain->computeProteins(hslData, membraneDiffusionRates, cellLength);

            //HSL WRITE:
            for(size_t i(0); i< hslSolutionVector.size(); i++)
            {
//                    if(bool(eQ::data::parameters["PETSC_SIMULATION"]))
//                        writeHSL(deltaHSL[0], petscSolutionVector[0], petscLookupTable[0], cellPoints);
//                    else
                writeHSL(deltaHSL[i], hslSolutionVector[i], dofLookupTable[i], cellPoints);
            }

            setDiffusionTensor(thisCell->getAngle(), cellPoints);

        };


        //now search the rectangular area computed above and call pointIsInCell() method on each point:
        findInteriorPoints(nodesPerMicron, cellPoints);
        //note: the vector contains points that are verified inside the cell
        averagePointsPerCell += cellPoints.size();

        updateCellModel();

        recordCelldata();

        //update the spring rest-length extension
        //note: possibly will be altered by genetic circuit, pressure, [protein], etc...
        //should pass growth rate, which should be computed in the ABM class
        thisCell->updateGrowth();
    }

    cellList.average(averagePointsPerCell);//divides RHS by the cell count (if count > 0)

/*
    if(bool(eQ::data::parameters["PETSC_SIMULATION"]))
    {
        //COPY FENICS I/O TO PETSC VECTOR FOR TESTING:
        for(size_t grid(0); grid< hslSolutionVector.size(); grid++)
        {
            for (size_t i(0); i < (nodesHigh); i++)
            {
                for (size_t j(0); j < (nodesWide); j++)
                {
                    auto f = dofLookupTable[grid]->grid[i][j];
                    auto p = petscLookupTable[grid]->grid[i][j];
//                    petscSolutionVector[grid].at(p) =  hslSolutionVector[grid].at(f);
                    hslSolutionVector[grid]->at(f) =  petscSolutionVector[grid]->at(p);
                }
            }
        }
    }
*/

}

eQabm::~eQabm()
{
    cellPointers.reset();
    cellList.clear();
    trap.reset();
    habitat.reset();
}
void eQabm::recordDivisionEvent(double timeStamp, std::shared_ptr<Ecoli> cell, std::shared_ptr<Ecoli> daughter)
{
    if(nullptr == cell)//initial cell at t=0: set parent to -ID (length 0)
    {
        divisionList.push_back(
                    std::make_tuple(
                            timeStamp,
                            daughter->getParentID(), 0.0,//parent length=0
                            daughter->getCellID(), daughter->getLengthMicrons()
                    )
        );
    }
    else
    {
        divisionList.push_back(
                    std::make_tuple(
                            timeStamp,
                            cell->getCellID(), cell->getLengthMicrons(),
                            daughter->getCellID(), daughter->getLengthMicrons()
                    )
        );
    }
}
namespace eQ {
    std::ostream &operator<<(std::ostream &os,  const eQ::Cell::strainType &type)
    {
        if(eQ::Cell::strainType::ACTIVATOR == type)
            os<<"ACTIVATOR";
        else if(eQ::Cell::strainType::REPRESSOR == type)
            os<<"REPRESSOR";
        else if(eQ::Cell::strainType::X == type)
            os<<"X";
        else if(eQ::Cell::strainType::Y == type)
            os<<"Y";
        else if(eQ::Cell::strainType::Z == type)
            os<<"Z";
        else os<< "";
        return os;
    }
}
