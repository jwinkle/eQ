#include "eQabm.h"


//NOTE:  ABM will have two spatial resolutions:
//one for full sub-micron position to map HSL
//and one for general data recording (1 micron or greater typically).
//the two resolutions are passed in simData structure as:
//nodesPerMicron (full resolution) and nodesPerMicronData (reduced for data recording)
//nodesHigh and nodesWide use full resolution for making the pointer array grid

double eQabm::rn()
{//interface for uniform random numbers in [0,1]
    return zeroOne->randomNumber();
}

eQabm::eQabm(const struct eQabm::params &p)
        :Params(p)
{
    dt = double(eQ::parameters["dt"]);
    nodesHigh = size_t(eQ::parameters["simulationTrapHeightMicrons"])*size_t(eQ::parameters["nodesPerMicronSignaling"])+ 1;
    nodesWide = size_t(eQ::parameters["simulationTrapWidthMicrons"])*size_t(eQ::parameters["nodesPerMicronSignaling"])+ 1;
    nodesHighData = size_t(eQ::parameters["simulationTrapHeightMicrons"])*size_t(eQ::parameters["nodesPerMicronData"])+ 1;
    nodesWideData = size_t(eQ::parameters["simulationTrapWidthMicrons"])*size_t(eQ::parameters["nodesPerMicronData"])+ 1;

    trapHeightMicrons   = double(eQ::parameters["simulationTrapHeightMicrons"]);
    trapWidthMicrons    = double(eQ::parameters["simulationTrapWidthMicrons"]);
    nodesPerMicron      = double(eQ::parameters["nodesPerMicronSignaling"]);
    nodesPerMicronData  = double(eQ::parameters["nodesPerMicronData"]);

    std::cout<<"Creating eQabm habitat and trap..."<<std::endl;

    //create the habitat:
    habitat = std::make_shared<cpmHabitat>(dt);

    //keep old params type for now:
    cpmTrap::params trapParams;
        trapParams.space = habitat->_space;
    //create the trap:
    trap = std::make_shared<cpmTrap>(trapParams);
    habitat->setTrap(trap);//don't omit this

    zeroOne = std::make_shared<eQ::uniformRandomNumber>(Params.seedValue);//[0,1]

    std::cout<<"Created eQabm habitat and trap..."<<std::endl;
    std::cout<<"eQ::w,h:"<<trapWidthMicrons<<", "<<trapHeightMicrons<<std::endl;
    std::cout<<"habitat->_trap->w,h:"<<habitat->_trap->w<<", "<<habitat->_trap->h<<std::endl;
    stabilityScaling = dt/0.001;//needed for stability of the ABM growth, per model (see Phys.Bio paper appendix)
//    stabilityScaling = Params.simData.dt/0.001;//needed for stability of the ABM growth, per model (see Phys.Bio paper appendix)
    std::cout<<"stabilityScaling: "<<stabilityScaling<<std::endl;
    std::cout<<"Number of data files: "<<(*Params.dataFiles).size()<<std::endl;
    std::cout<<std::endl;

        gridDataCounter =
                    std::make_shared<eQ::gridFunction<size_t>>(nodesHighData, nodesWideData);

    cellCount = 0;
    eraseCounter=0;
    std::cout<<"creating cellPointers: "<<nodesHigh<<", "<<nodesWide<<std::endl;
        //Note: uses full resolution of grid (not the data recording resolution)
        cellPointers =
                std::make_shared<eQ::gridFunction<std::shared_ptr<eColi>>>(nodesHigh, nodesWide);

    //grids for defining anisotropic diffusion tensor
    D11grid = std::make_shared<
            eQ::gridFunction<double>>(nodesHigh, nodesWide);
    D22grid = std::make_shared<
            eQ::gridFunction<double>>(nodesHigh, nodesWide);
    D12grid = std::make_shared<
            eQ::gridFunction<double>>(nodesHigh, nodesWide);

    //compute the max. # of nodes from a pole center to the boundary of the cell (for searching for interior points)
    nodesToEdge = size_t(floor(nodesPerMicron * DEFAULT_CELL_WIDTH_MICRONS/2.0));

    maxSpringCompression = 0.0;

}

void eQabm::initCells(int numCellsToInit)
{
    if(!cellList.empty())
        std::cout<<"cellList not empty!"<<std::endl;
    cellList.clear();
    cellCount = 0;
    binBuffer.assign(numBins, 0);
    binBuffer2.assign(numBins, 0);

    eColi::params cellParams;
    cellParams.space = habitat->_space;

    double rx, ry, ra;

    //set default initial cell parameters
    cellParams.width = DEFAULT_CELL_WIDTH_MICRONS;
    cellParams.vx = 0.0;
    cellParams.vy = 0.0;
    cellParams.dRL = 0.1*dt;//nominal growth rate is 0.1 um/min.
    cellParams.stabilityScaling = stabilityScaling;
    cellParams.meanDivisionLength       = DEFAULT_DIVISION_LENGTH_MICRONS;//default, possibly reset below:
    cellParams.divisionLength           = DEFAULT_DIVISION_LENGTH_MICRONS;//default, possibly reset below:
    cellParams.doublingPeriodMinutes    = DEFAULT_CELL_DOUBLING_PERIOD_MINUTES;

    int countA=0, countB=0;

    //to speed up cell init:
    int initColumns=0, initRows=0, columnCounter=0, rowCounter=0;
    if( ("MODULUS_1" == eQ::parameters["simType"])
            || ("MODULUS_2" == eQ::parameters["simType"])
            || ("ASPECTRATIO_INVASION" == eQ::parameters["simType"])
            )
    {
        initRows = (int(trapHeightMicrons) - 0)/4;
        initColumns = (int(trapWidthMicrons) - 0)/2;
        if(eQ::parameters["trapType"] == "H_TRAP")
        {
            initRows = 1;
            initColumns = (int(trapWidthMicrons) - 0)/4;
        }
        std::cout<<"initRows initColumns: "<<initRows<<", "<<initColumns<<std::endl;
    }    

    std::cout<<"void eQabm::initCells(int numCellsToInit)"<<std::endl;

    for(int i(0); i< numCellsToInit; i++)
    {
        rx = rn() * trapWidthMicrons;
        ry = rn() * trapHeightMicrons;
        ra = rn() * 2*M_PI;
            cellParams.x = rx;
            cellParams.y = ry;
            cellParams.angle = ra;
        //set default division length (set in main) for all cells (may be reset below)
        cellParams.meanDivisionLength =
                double(eQ::parameters["defaultAspectRatioFactor"]) * DEFAULT_DIVISION_LENGTH_MICRONS;

        if( ("MODULUS_1" == eQ::parameters["simType"])
            || ("MODULUS_2" == eQ::parameters["simType"])
            || ("ASPECTRATIO_INVASION" == eQ::parameters["simType"])
            )
        {
            //code to populate the cells on a grid to speed-up init:
            if(initRows == rowCounter) continue;//skips to next for loop iteration, so runs that loop out
            cellParams.x =  5 + (columnCounter++)*2;
            cellParams.y =  5 + (rowCounter)*4;
            if(initColumns == columnCounter){rowCounter++; columnCounter = 0;}
            cellParams.angle = M_PI/2.0;
            if(eQ::parameters["trapType"] == "H_TRAP")
            {
                cellParams.angle = 0.0;
                cellParams.x =  (columnCounter++)*3;
                cellParams.y =  trapHeightMicrons/2.0;
            }

            if("ASPECTRATIO_INVASION" == eQ::parameters["simType"])
            {
                double halfWidth = trapWidthMicrons/2.0;
                cellParams.strainType = (cellParams.x > halfWidth)
//                        ? eQ::strainType::REPRESSOR : eQ::strainType::ACTIVATOR;
                        ? eQ::strainType::ACTIVATOR : eQ::strainType::ACTIVATOR;
            }
            else
            {
                cellParams.strainType = eQ::strainType::ACTIVATOR;
                cellParams.dsoParams = Params.pA;
            }
        }
        //todo:  separate these:
        else if(
           ("DUALSTRAIN_OSCILLATOR" == eQ::parameters["simType"])
            || ("SENDER_RECEIVER" == eQ::parameters["simType"])
            || ("INDUCED_SENDER_RECEIVER" == eQ::parameters["simType"])
            || ("DUAL_SENDER_RECEIVER" == eQ::parameters["simType"])
           || ("STATIC_ASPECTRATIO" == eQ::parameters["simType"])
                || ("INDUCED_DYNAMIC_ASPECTRATIO" == eQ::parameters["simType"])
           )
        {
            if("ACTIVATOR_ONLY" == eQ::parameters["cellInitType"])
            {
                cellParams.strainType = eQ::strainType::ACTIVATOR;
            }
            if("RANDOM" == eQ::parameters["cellInitType"])
            {
                if(rn() > 0.5)
                {
                    cellParams.strainType = eQ::strainType::ACTIVATOR;
                    cellParams.dsoParams = Params.pA;
                }
                else
                {
                    cellParams.strainType = eQ::strainType::REPRESSOR;
                    cellParams.dsoParams = Params.pR;
                }
                std::cout<<"Strain: "<<cellParams.strainType<<" ";
            }
            else if("AB_HALF" == eQ::parameters["cellInitType"])
            {
                double halfWidth = trapWidthMicrons/2.0;
                if(rn() > 0.5)
                {
                    cellParams.strainType = eQ::strainType::ACTIVATOR;
                    cellParams.dsoParams = Params.pA;
                    //note:this is ignored if send/recv since factor=1.0:
                    if("STATIC_ASPECTRATIO" == eQ::parameters["simType"])
                        cellParams.meanDivisionLength = double( eQ::parameters["mutantAspectRatioScale"])*DEFAULT_DIVISION_LENGTH_MICRONS;
                    //set the two strains on opposite sides of the trap:
                    if(cellParams.x > halfWidth)
                        cellParams.x -= halfWidth;
                }
                else
                {
                    cellParams.strainType = eQ::strainType::REPRESSOR;
                    cellParams.dsoParams = Params.pR;
                    if(cellParams.x < halfWidth)
                        cellParams.x += halfWidth;
                }
            }
            else if("ABA_THIRDS" == eQ::parameters["cellInitType"])
            {
                double thirdWidth = trapWidthMicrons/3.0;
                if(rn() < 0.67)
                {//2/3 are flanking activator/mutant
                    cellParams.strainType = eQ::strainType::ACTIVATOR;
                    cellParams.dsoParams = Params.pA;
                    //note:this is ignored if send/recv since factor=1.0:
                    if("STATIC_ASPECTRATIO" == eQ::parameters["simType"])
                        cellParams.meanDivisionLength = double( eQ::parameters["mutantAspectRatioScale"])*DEFAULT_DIVISION_LENGTH_MICRONS;
                    //set the two strains on opposite sides of the trap:
                    if( (cellParams.x > thirdWidth) && ((cellParams.x < 2.0*thirdWidth)) )
                        cellParams.x += (rn() > 0.5) ? thirdWidth:-thirdWidth;
                }
                else
                {
                    cellParams.strainType = eQ::strainType::REPRESSOR;
                    cellParams.dsoParams = Params.pR;
                    if(cellParams.x < thirdWidth)
                        cellParams.x += thirdWidth;
                    else if(cellParams.x > 2.0*thirdWidth)
                        cellParams.x -= thirdWidth;
                }
            }
        }
        else if("SINGLEMUTANT_ASPECTRATIO" == eQ::parameters["simType"])
        {
            cellParams.strainType = eQ::strainType::REPRESSOR;
        }
        else
        {
            cellParams.strainType = eQ::strainType::ACTIVATOR;
        }


        //add some distribution over the intial cells' division length
        double dn = double(eQ::parameters["divisionNoiseScale"]);
        double noiseFactor = (rn()-0.5) * dn;

        cellParams.divisionLength = (1.0 + noiseFactor) * cellParams.meanDivisionLength;
        cellParams.length = 0.5*cellParams.divisionLength;

        std::cout<<"cell: "<<i<<" ("<<cellParams.x<<","<<cellParams.y<<") "<<cellParams.angle<<"=angle; strain = "
                <<cellParams.strainType<<std::endl;


        if("ASPECTRATIO_INVASION" == eQ::parameters["simType"])
        {
            cellParams.strain = std::make_shared<aspectRatioInvasionStrain>(
                        cellParams.strainType, double(eQ::parameters["dt"]), nodesPerMicron, Params.hslSolutionVector.size());
        }
        //create lineages, proteins, strains...:
        else if("DUALSTRAIN_OSCILLATOR" == eQ::parameters["simType"])//    if(eQ::simType::DUALSTRAIN_OSCILLATOR ==  Params.simData.simtype)
        {
            cellParams.strain = std::make_shared<Strain>(
                        cellParams.strainType, cellParams.dsoParams, double(eQ::parameters["dt"]), nodesPerMicron);
        }
        else if( ("SENDER_RECEIVER" == eQ::parameters["simType"])
            || ("DUAL_SENDER_RECEIVER" == eQ::parameters["simType"])
            || ("INDUCED_SENDER_RECEIVER" == eQ::parameters["simType"]) )
        {
//            cellParams.strain = std::make_shared<sendRecvStrain>(
//                        cellParams.strainType, cellParams.dsoParams, double(eQ::parameters["dt"]), nodesPerMicron);
            cellParams.strain = std::make_shared<sendRecvStrain>(
                        cellParams.strainType, double(eQ::parameters["dt"]), nodesPerMicron, Params.hslSolutionVector.size());
        }
        else if( ("MODULUS_1" == eQ::parameters["simType"])
            || ("MODULUS_2" == eQ::parameters["simType"]) )
        {
//            cellParams.strain = std::make_shared<MODULUSmodule>(
//                        cellParams.dsoParams, double(eQ::parameters["dt"]), nodesPerMicron);
            cellParams.strain = std::make_shared<MODULUSmodule>(
                        double(eQ::parameters["dt"]), nodesPerMicron, Params.hslSolutionVector.size());
        }


        //CREATE CELL OBJECT:
        auto thisCell = std::make_shared<eColi>(cellParams);

        (cellParams.strainType == eQ::strainType::ACTIVATOR)
                ? ++countA : ++countB;

        cellList.push_front(thisCell);
        ++cellCount;
        thisCell->setCellID(cellID_factory++);
        thisCell->parentID = -(thisCell->getCellID());//initial cells set parent to negative cellID.
        recordDivisionEvent(0.0, nullptr, thisCell);
    }
    auto totalCells = countA+countB;
    if(totalCells > 0)
    {
        strainRatio = double(countA)/double(totalCells);
        std::cout<<"Initial strain ratio: "<<strainRatio<<std::endl;
    }
}
void eQabm::stepChipmunk()
{
    habitat->stepSimulation(dt);
}
void eQabm::updateCellPositions(double simTime)
{//this is done while waiting for HSL data transfer (single-threaded)

    //UPDATES PHYSICS MODEL RESULTS (cell length, position, angle, spring compression, velocity, ...)
    //REMOVES CELLS OUT-OF-BOUNDS
    //DIVIDES CELLS THAT REACH THEIR DIVISION LENGTH TARGET
    //UPDATES LIST DATA STRUCTURE, REMOVING (EXIT) OR ADDING (DIVISION) CELLS AS NEEDED

    std::vector<std::shared_ptr<eColi>> newDaughterCells;
//    refreshRandomNumberStack();

    auto icells     = cellList.begin();
    auto icellsPrev = cellList.before_begin();//needed for forward_list C++ data structure
    while(icells != cellList.end())
    {
        //TRAVERSE LIST OF CELLS: (uses C++ STL forward list data structure)

        std::shared_ptr<eColi> cell = (*icells);
        if(nullptr == cell){std::cout<<"(nullptr == cell)"<<std::endl; continue;}

        //==============================================
        //  UPDATE PHYSICS MODEL AND CHECK FOR DIVISION:
        //==============================================
        //CELL MODEL:
//        std::shared_ptr<eColi> daughter = cell->updatePhysicsModel(rns);//pass random number stack by reference
        std::shared_ptr<eColi> daughter = cell->updatePhysicsModel(*zeroOne);//pass random number stack by reference
        if(nullptr != daughter)
        {//divided cell, check if outside trap (but record division to parent first):
            daughter->setCellID(cellID_factory++);
            daughter->parentID = cell->getCellID();
            recordDivisionEvent(simTime, cell, daughter);
            if(habitat->outsideTrap(daughter->cpmCell))
            {
                daughter.reset();
            }
            else
            {
                newDaughterCells.push_back(daughter);
            }
        }
        //==============================================
        //  CHECK FOR OUTSIDE OF BOUNDARY:
        //==============================================
        if(habitat->outsideTrap(cell->cpmCell))
        {
            //push on test data before cell is deleted:
//            cell->finalizeCell(highResolutionDataVector);
            icells = cellList.erase_after(icellsPrev);
            --cellCount;
            ++eraseCounter;
        }
        else
        {//necessary for forward_list to have 2 pointers:
            icellsPrev = icells;
            ++icells;
        }
    }
    for(auto &cell : newDaughterCells)
    {//add daughter cells to main cell list:
        cellList.push_front(cell);
        ++cellCount;
    }
    if("SINGLEMUTANT_ASPECTRATIO" == eQ::parameters["simType"])
    {
        if(mutantTriggerFlag)
        {
            //randomly select a cell from the count (flag reset below)
            mutantCellNumber = size_t(std::floor(rn()*cellCount));
        }
    }

    //new positions and daughter cells are set; update x,y grid for pointers:
//    std::cout<<"updateCellGrid()"<<std::endl;
    updateCellGrid(cellList.begin(), cellList.end());
//    std::cout<<"updateCellGrid() END"<<std::endl;
}
void eQabm::updateCellGrid(fli_t begin, fli_t end)//passes forward list iterator type (fli_t)
{//populates the 2D array of cell pointers based on cell position (full grid resolution)
    cellPointers->clear();
    double x,y;    size_t i,j;    std::pair<size_t, size_t> index;
    double springEnergy_L2 = 0.0;
    int countA=0, countB=0;
    size_t counter=0;


    //iterate cell list structure
    for(auto &cell(begin); cell!=end; cell++)
    {
        x = (*cell)->getCenter_x();
        y = (*cell)->getCenter_y();
        index = eQ::ij_from_xy(x,y,nodesPerMicron);
        i = index.first;  j = index.second;

        if(nullptr != cellPointers->grid[i][j])
            std::cout<<"overwrite of cell grid pointer: (x,y) = "<<x<<", "<<y<<std::endl;
        cellPointers->grid[i][j] = (*cell);


        springEnergy_L2 += (*cell)->cpmCell->compression * (*cell)->cpmCell->compression;
        (eQ::strainType::ACTIVATOR == (*cell)->Params.strainType) ? countA++:countB++;
        if((*cell)->cpmCell->compression > maxSpringCompression) maxSpringCompression = (*cell)->cpmCell->compression;
    }
    if(timeSeriesDataTrigger)
    {
        timeSeriesDataTrigger = false;
        compressionTimeSeries.push_back(springEnergy_L2);
        angleBinTimeSeries.push_back(binBuffer);
        angleBinTimeSeries2.push_back(binBuffer2);
    }
    auto totalCells = countA+countB;
    if(totalCells > 0)
    {
        strainRatio = double(countA)/double(totalCells);
    }
}
void eQabm::updateCellModels(size_t numThreads)
{//launches multi-threaded cell model updates
    for(auto &dataFile : *(Params.dataFiles))
    {
        if( (eQ::dataParameter::SPRING_COMPRESSION != dataFile.second)
                && (eQ::dataParameter::CELL_VELOCITY != dataFile.second) )
        {
            dataFile.first->dataGrid->clear();
        }
    }
    //SINGLE-THREADED CELL UPDATE PATH:
    if(true)
//        if(1 == numThreads)
    {
        updateCells(cellList.begin(), cellList.end());
        return;
    }

    //MULTI-THREADED TRAVERSAL OF UPDATED CELL LIST, TO UPDATE ABM CELL MODELS (this works and is tested):
    size_t chunk = cellCount/numThreads;
    size_t counter=0;
    size_t chunks=0;//counter of how many chunks

    auto cell = cellList.begin();
    tiBegins.clear();  tiEnds.clear();
    tiBegins.push_back(cell);//first begin is front of list
    while(cell != cellList.end())
    {
        ++cell;
        if(++counter == chunk)
        {
            if(++chunks == numThreads)
                cell = cellList.end();
            else
            {
                tiEnds.push_back(cell);
                tiBegins.push_back(cell);
            }
            counter = 0;
        }
    }
    tiEnds.push_back(cell);

    std::vector<std::thread> threads;
    for(size_t i(0); i<numThreads; i++)
    {
    //            auto start = i*chunk;
    //            auto end = (i+1)*chunk;
        //creates and launches threads, calling function passed w/ parameters
        threads.push_back(std::thread(&eQabm::updateCells, this, tiBegins[i], tiEnds[i]));
    }

    //the calling thread runs the final chunk:
    //        stepCells(numThreads*chunk, rows, simulation.get());
    //wait for all threads to complete:
    for (auto& th : threads) th.join();

}
void eQabm::updateCells(fli_t begin, fli_t end)
{//cell update function called with distributed partition of the cell list:
 //typically called multi-threaded with list range and no risk of race conditions for cells
    double x,y;
    size_t i,j;
    std::pair<size_t, size_t> index;

    //ANISOTROPIC DIFFUSION CODE:
    double Dx = double(eQ::parameters["AnisotropicDiffusion_Axial"]);
    double Dy = double(eQ::parameters["AnisotropicDiffusion_Transverse"]);
    //init the grid to isotropic diffusion: (these grids copy directly to the .ufl file tensor)
    //note:  these are 1D vectors using standard translation for (i,j) entries using the eQ helper functions
    D11grid->assign(1.0);//diagonal scaling
    D22grid->assign(1.0);//diagonal scaling
    D12grid->assign(0.0);//symmetric, off-diagonal scaling

    averagePointsPerCell = 0.0;
    int cellCounter = 0;

    for(auto &cell(begin); cell!=end; cell++)
    {
        auto thisCell = (*cell);
        double cellLength = thisCell->getLengthMicrons();
        x = thisCell->getCenter_x();
        y = thisCell->getCenter_y();
        cellCounter++;
            index = eQ::ij_from_xy(x,y,nodesPerMicronData);//uses data scaling
            i = index.first;  j = index.second;
            //keep tally of "hits" at this grid point for averaging:
            gridDataCounter->grid[i][j]++;

        auto findInteriorPoints = [&](double npm, std::vector<std::pair<size_t, size_t>> &cellPoints)
        {
            //implicitly passed by reference: x, y, thisCell, nodesToEdge, nodesHigh, nodesWide
            size_t i,j;
            //CONVERT THE GRID (X,Y) VALUE TO A (J,I) GRID POSITION FOR HSL SIGNALING INTERFACE
            index = eQ::ij_from_xy(x,y,npm);//uses  grid scaling as passed by "nodes per micron" npm
            i = index.first;  j = index.second;//location of cell center, stored for use in all that follows...

            //get location of pole centers:
            //NOTE: need to check if they are out of bounds! (only cell centers are checked for cell removal)
            auto poleA = eQ::ij_from_xy(thisCell->polePositionA, npm);//uses full grid scaling
            auto poleB = eQ::ij_from_xy(thisCell->polePositionB, npm);//uses full grid scaling

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
//            std::vector<std::pair<size_t, size_t>> cellPoints;//vector of (i,j) pairs
            for(size_t pointi(i1); pointi<=i2; pointi++)
            {
                for(size_t pointj(j1); pointj<=j2; pointj++)
                {
                    auto xypoint = eQ::xy_from_ij(pointi, pointj, npm);
                    if(thisCell->cpmCell->pointIsInCell(xypoint))
                        cellPoints.push_back(std::make_pair(pointi, pointj));
                }
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
                    D11grid->grid[point.first][point.second] = Dx*cos2t + Dy*sin2t;
                    D22grid->grid[point.first][point.second] = Dx*sin2t + Dy*cos2t;
                    D12grid->grid[point.first][point.second] = (Dx-Dy)*sincost;
                }
            }
        };
        auto readHSL = [&](std::vector<double> &HSLvector, eQabm::HSLgrid &HSLlookup, std::vector<std::pair<size_t, size_t>> &cellPoints)
        {
            //read from each interior point of the cell (computed above) and equal-average values read:
            double HSL = 0.0;
            for(auto &point : cellPoints)
            {
                //fenics implementation:
                auto location = HSLlookup->grid[point.first][point.second];//translates the i,j position to the fenics DOF entry in the solution vector
                HSL += HSLvector.at(location);
            }
            return HSL/double(cellPoints.size());
        };
        auto writeHSL = [&](double HSL, std::vector<double> &HSLgrid, eQabm::HSLgrid &HSLlookup, std::vector<std::pair<size_t, size_t>> &cellPoints)
        {//NOTE:  HSL passed in is in units of concentration [nM]...must scale to extra-cellular volume to ensure convservation of molecule #
            //we need the fraction of volume outside of cell relative to total rectangular volume 1um wide x 1um tall x cell length (um)
            //[HSL] is the delta concentraion of HSL for the cell; need to determine the amount to put on one grid point:
            //[HSL]*cellVolume=#HSL; to update 1um^2 2D grid point with #HSL molecules and xum^2 extra-cellular volume per grid point:
            //H1, new concentration at grid point;  H0, previous.
            //H1 = [(H0*xum^2) + #HSL]/xum^2 = H0 + #HSL/xum^2, where #HSL = [HSL]*cellVolume

            double oneNanoMolarToMoleculeNumber = 1.0/eQ::proteinNumberToNanoMolar(1.0, cellLength);
            double numberHSL = HSL * oneNanoMolarToMoleculeNumber;
            double updatePerSquareMicron = numberHSL/eQ::computeExtraCellularVolumeFraction(cellLength);
            double updateForOneGridPoint = updatePerSquareMicron * nodesPerMicron * nodesPerMicron;

            //distribute over interior points of the cell:
            double dHSL = updateForOneGridPoint/double(cellPoints.size());
            for(auto &point : cellPoints)
            {
                //fenics implementation:
                auto location = HSLlookup->grid[point.first][point.second];//translates the i,j position to the fenics DOF entry in the solution vector
                HSLgrid.at(location) += dHSL;
            }
        };
        //Note: define new parameters and set filenames in eQ.h
        //here, we populate the data structure (fenics recording node) with the cell data
        //some data is averaged between recording to file, others are one-shot (check +=)
        auto recordCelldata = [&]()
        {
            index = eQ::ij_from_xy(x,y,nodesPerMicronData);//uses data scaling
            i = index.first;  j = index.second;

            for(auto &dataFile : *(Params.dataFiles))
            {
                switch (dataFile.second){
                //spring compression and velocity are averaged for now:
                case eQ::dataParameter::SPRING_COMPRESSION:
                    dataFile.first->dataGrid->grid[i][j] += (*cell)->getSpringCompression();
                break;
                case eQ::dataParameter::CELL_VELOCITY:
                    dataFile.first->xdataGrid->grid[i][j] += (*cell)->cpmCell->velocity.x;
                    dataFile.first->ydataGrid->grid[i][j] += (*cell)->cpmCell->velocity.y;
                break;
                case eQ::dataParameter::C4:
                    dataFile.first->dataGrid->grid[i][j]
//                            = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(H), cellLength);
//                            = eQ::proteinNumberToNanoMolar(100.0, cellLength);
                            = eQ::proteinNumberToNanoMolar((*cell)->strain->iHSL[0], cellLength);
//                            = (*cell)->strain->tHSL[0];
                    break;
                case eQ::dataParameter::C14:
                    dataFile.first->dataGrid->grid[i][j]
//                            = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(I), cellLength);
//                            = eQ::proteinNumberToNanoMolar(200.0, cellLength);
                            = eQ::proteinNumberToNanoMolar((*cell)->strain->iHSL[1], cellLength);
//                            = (*cell)->strain->tHSL[2];
                break;
                case eQ::dataParameter::FP:
                    if( ("MODULUS_1" == eQ::parameters["simType"]) ||
                        ("MODULUS_2" == eQ::parameters["simType"]))
                    {
//                        auto pmodulus = std::dynamic_pointer_cast<MODULUSmodule>((*cell)->strain);
                        dataFile.first->dataGrid->grid[i][j]
//                                = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(FP), pmodulus->filteredCellLength);
                                = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(FP), cellLength);
                    }
                    else
                        dataFile.first->dataGrid->grid[i][j]
                            = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(FP), cellLength);
                break;
                case eQ::dataParameter::CFP:
                    dataFile.first->dataGrid->grid[i][j] = (eQ::strainType::ACTIVATOR == (*cell)->strain->whichType)
                            ? eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(FP), cellLength) : 0.0;
                break;
                case eQ::dataParameter::YFP:
                    dataFile.first->dataGrid->grid[i][j] = (eQ::strainType::REPRESSOR == (*cell)->strain->whichType)
                            ? eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(FP), cellLength) : 0.0;
                break;
                case eQ::dataParameter::SYN:
                    dataFile.first->dataGrid->grid[i][j]
                            = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(S), cellLength);
                break;
                case eQ::dataParameter::LACI:
                    dataFile.first->dataGrid->grid[i][j]
                            = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(L), cellLength);
                break;

                case eQ::dataParameter::AIIA:
                    dataFile.first->dataGrid->grid[i][j]
                            = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(A), cellLength);
                break;
                case eQ::dataParameter::MFP:
                    dataFile.first->dataGrid->grid[i][j]
                            = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(M), cellLength);
                break;
                case eQ::dataParameter::C4RHL:
                    dataFile.first->dataGrid->grid[i][j]
                            = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(HR), cellLength);
                break;
                case eQ::dataParameter::RHL_T:
                    dataFile.first->dataGrid->grid[i][j]
                            = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(R_T), cellLength);
                break;
                case eQ::dataParameter::MODULUS_S:
                    dataFile.first->dataGrid->grid[i][j]
                            = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(S), cellLength);
                break;
                case eQ::dataParameter::MODULUS_H:
                    dataFile.first->dataGrid->grid[i][j]
                            = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(LEGI_A), cellLength);
                break;
                case eQ::dataParameter::MODULUS_R:
                    dataFile.first->dataGrid->grid[i][j]
                            = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(FP), cellLength);
                    break;
                case eQ::dataParameter::DTENSOR_11:
                    dataFile.first->dataGrid->grid[i][j]
                            = D11grid->grid[i][j];
//                    = (*cell)->strain->tHSL[1];
                    break;
                case eQ::dataParameter::DTENSOR_22:
                    dataFile.first->dataGrid->grid[i][j]
                            = D22grid->grid[i][j];
//                    = (*cell)->strain->tHSL[3];
                    break;
                case eQ::dataParameter::DTENSOR_12:
                    dataFile.first->dataGrid->grid[i][j]
                            = D12grid->grid[i][j];
                    break;
                default: break;
                }//end switch
            }//end for
        };//end lambda recordCelldata

        if("NOTRAP" == eQ::parameters["trapType"])                  {}
        if("NO_SIGNALING" == eQ::parameters["simType"])             {}
        else if("STATIC_ASPECTRATIO" == eQ::parameters["simType"])  {}
        else
        {
            //CONVERT THE GRID (X,Y) VALUE TO A (J,I) GRID POSITION FOR HSL SIGNALING INTERFACE
            index = eQ::ij_from_xy(x,y,nodesPerMicron);//uses full grid scaling
            i = index.first;  j = index.second;//location of cell center, stored for use in all that follows...

            //JW: compute petsc index from x,y here:
            //double petscLocationC4 = ...

            //now search the rectangular area computed above and call pointIsInCell() method on each point:
            std::vector<std::pair<size_t, size_t>> cellPoints;//vector of (i,j) pairs
            findInteriorPoints(nodesPerMicron, cellPoints);
            //note: the vector contains points that are verified inside the cell
            //for checking this is working:
            averagePointsPerCell += cellPoints.size();


            if("INDUCED_DYNAMIC_ASPECTRATIO" == eQ::parameters["simType"])
            {
                if( (aspectRatioInduction) && (eQ::strainType::ACTIVATOR == (*cell)->Params.strainType) )
                {
                    (*cell)->Params.meanDivisionLength = double( eQ::parameters["mutantAspectRatioScale"])
                                                       * double( eQ::parameters["defaultAspectRatioFactor"])
                                                       * DEFAULT_DIVISION_LENGTH_MICRONS;
                    (*cell)->Params.divisionLength = (*cell)->Params.meanDivisionLength;//set to divde at this length immediately
                }
            }
            if("DUALSTRAIN_OSCILLATOR" == eQ::parameters["simType"])
            {
                //fenics implementation:
//                auto c4 = readHSL(Params.c4grid, Params.c4lookup, cellPoints);
//                auto c14 = readHSL(Params.c14grid, Params.c14lookup, cellPoints);
                auto c4 = readHSL(Params.hslSolutionVector[0], Params.dofLookupTable[0], cellPoints);
                auto c14 = readHSL(Params.hslSolutionVector[1], Params.dofLookupTable[1], cellPoints);

                auto deltaHSL = (*cell)->strain->computeProteins(c4, c14,(*cell)->getLengthMicrons());

                writeHSL(deltaHSL[0], Params.hslSolutionVector[0], Params.dofLookupTable[0], cellPoints);
                writeHSL(deltaHSL[1], Params.hslSolutionVector[1], Params.dofLookupTable[1], cellPoints);

                //13May.2019:  cell size modulation connected to dso:
//                if(eQ::strainType::REPRESSOR == (*cell)->Params.strainType)
                if(true)
                {
                    double threshold = double(eQ::parameters["dso_hslThresh"]);
                    double newDivisionLength = (*cell)->Params.meanDivisionLength;

                    double ftsZ = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(H), cellLength);
//                    double ftsZ = (*cell)->strain->conc[H]/(*cell)->getLengthMicrons();

                    if (ftsZ > threshold)
                     newDivisionLength *= (2.0/3.0);//scale

                    (*cell)->Params.divisionLength = newDivisionLength;
                }
            }
            else if("ASPECTRATIO_INVASION" == eQ::parameters["simType"])
            {
                double aspectRatioThresh = double( eQ::parameters["aspectRatioThresholdHSL"]);
                double aspectRatioScaling = double( eQ::parameters["defaultAspectRatioFactor"]);

                if(aspectRatioInduction)
                {
                    double hslValue = (eQ::strainType::REPRESSOR == thisCell->Params.strainType)
                            ? thisCell->strain->iHSL[1] : thisCell->strain->iHSL[3];

                    if(hslValue > aspectRatioThresh)
                        aspectRatioScaling *= double(eQ::parameters["mutantAspectRatioScale"]);
                }                

                thisCell->Params.meanDivisionLength = aspectRatioScaling * DEFAULT_DIVISION_LENGTH_MICRONS;
//                    (*cell)->Params.divisionLength = (*cell)->Params.meanDivisionLength;//set to divde at this length immediately


                std::vector<double> hslData;
                std::vector<double> membraneDiffusion;
                //HSL READ:
                for(size_t i(0); i< Params.hslSolutionVector.size(); i++)
                {
//                    hslData.push_back(readHSL(Params.hslSolutionVector[i].get(), Params.dofLookupTable[i], cellPoints));
                    hslData.push_back(readHSL(Params.hslSolutionVector[i], Params.dofLookupTable[i], cellPoints));
                }
                //TODO: move these to main where they are defined with the diffusion coeff. of each HSL
                membraneDiffusion.push_back(3.0);
                membraneDiffusion.push_back(2.1);
                membraneDiffusion.push_back(3.0);
                membraneDiffusion.push_back(2.1);

                auto deltaHSL = thisCell->strain->computeProteins(hslData, membraneDiffusion, cellLength);

                //HSL WRITE:
                for(size_t i(0); i< Params.hslSolutionVector.size(); i++)
                {
//                    writeHSL(deltaHSL[i], Params.hslSolutionVector[i].get(), Params.dofLookupTable[i], cellPoints);
                    writeHSL(deltaHSL[i], Params.hslSolutionVector[i], Params.dofLookupTable[i], cellPoints);
                }

                setDiffusionTensor(thisCell->getAngle(), cellPoints);

            }
            else if (
                     ("SENDER_RECEIVER" == eQ::parameters["simType"])
                    || ("INDUCED_SENDER_RECEIVER" == eQ::parameters["simType"])
                    || ("DUAL_SENDER_RECEIVER" == eQ::parameters["simType"]) )

            {
//                (*cell)->strain->setPressureK50(10.0);
//                if(pressureInductionFlag)
//                    (*cell)->strain->setPressureValue((*cell)->getSpringCompression());
//                else
//                    (*cell)->strain->setPressureValue(0.0);

                if ("DUAL_SENDER_RECEIVER" == eQ::parameters["simType"])
                {
//                    auto c4 = readHSL(Params.c4grid, Params.c4lookup, cellPoints);
//                    auto c14 = readHSL(Params.c14grid, Params.c14lookup, cellPoints);
                    auto c4 = readHSL(Params.hslSolutionVector[0], Params.dofLookupTable[0], cellPoints);
                    auto c14 = readHSL(Params.hslSolutionVector[1], Params.dofLookupTable[1], cellPoints);

                    auto deltaHSL = (*cell)->strain->computeProteins(c4, c14,(*cell)->getLengthMicrons());

//                    writeHSL(deltaHSL[0], Params.c4grid, Params.c4lookup, cellPoints);
//                    writeHSL(deltaHSL[1], Params.c14grid, Params.c14lookup, cellPoints);
                    writeHSL(deltaHSL[0], Params.hslSolutionVector[0], Params.dofLookupTable[0], cellPoints);
                    writeHSL(deltaHSL[1], Params.hslSolutionVector[1], Params.dofLookupTable[1], cellPoints);

                    enum concentrations  whichHSL = (eQ::strainType::REPRESSOR == (*cell)->Params.strainType) ? H : I;
                    double threshold = double(eQ::parameters["hslThresh"]);
                    double newDivisionLength = (*cell)->Params.meanDivisionLength;

                    double ftsZ = eQ::proteinNumberToNanoMolar((*cell)->strain->getProteinNumber(whichHSL), cellLength);

                    if (ftsZ > threshold)
                     newDivisionLength *= (2.0/3.0);//scale

                    (*cell)->Params.divisionLength = newDivisionLength;
                }
                else
                {
                    double iptg = (x/trapWidthMicrons);
                    eQ::parameters["MODULUS_IPTG"] = iptg;

                    std::vector<double> hslData;
                    std::vector<double> membraneDiffusionRates;
                    //HSL READ:
                    for(size_t i(0); i< Params.hslSolutionVector.size(); i++)
                    {
                        if(bool(eQ::parameters["PETSC_SIMULATION"]))
                            hslData.push_back(readHSL(Params.petscSolutionVector[i], Params.petscLookupTable[i], cellPoints));
                        else
                            hslData.push_back(readHSL(Params.hslSolutionVector[i], Params.dofLookupTable[i], cellPoints));
                    }

                    //TODO: move these to main where they are defined with the diffusion coeff. of each HSL
                    membraneDiffusionRates = std::vector<double>(eQ::parameters["membraneDiffusionRates"].get<std::vector<double>>());

                    auto deltaHSL = thisCell->strain->computeProteins(hslData, membraneDiffusionRates, cellLength);

                    //HSL WRITE:
                    if(bool(eQ::parameters["PETSC_SIMULATION"]))
                        writeHSL(deltaHSL[0], Params.petscSolutionVector[0], Params.petscLookupTable[0], cellPoints);
                    else
                        writeHSL(deltaHSL[0], Params.hslSolutionVector[0], Params.dofLookupTable[0], cellPoints);


                    setDiffusionTensor(thisCell->getAngle(), cellPoints);
                }
            }
            else if( ("MODULUS_1" == eQ::parameters["simType"]) || ("MODULUS_2" == eQ::parameters["simType"]) )
            {
                //    a = v(4)/D;
                //    cxR = c0 / (1 - exp(a*L))  ...
                //        * (exp(a*x) - exp(a*L)) + offset;
                double cl = double(eQ::parameters["channelLengthMicronsLeft"]);//already scaled
                double cr = double(eQ::parameters["channelLengthMicronsRight"]);
                double tw = trapWidthMicrons;
                double xt = cl+cr+tw;

                double vx = double(eQ::parameters["trapChannelLinearFlowRate"])/double(eQ::parameters["lengthScaling"]);
                //um^2/min for IPTG (relative to C4HSL a la Pai and You, using molecular weight vs. C4HSL)
                double D = 3.0e4 * sqrt(159.0/238.0) * double(eQ::parameters["diffusionScaling"]);
                double a = vx/D;
                double expaL = exp(a*xt);//exp(ax) evaluated at x=L (right side of device)
                double expax = exp(a*(x+cl));//exp(ax) evaluated at cell position x (add left-side offset)

                double iptgOffset   = 0.0;//low side of gradient at media channel right (Modulus c0 offset of input gradient)
                double iptgLeft     = +100.0;//high side of gradient at media channel left (+ difference from right side)

                double thisC = (iptgLeft/(1.0 - expaL)) * (expax - expaL) + iptgOffset;

                    if ("MODULUS_2" == eQ::parameters["simType"])
                    {
//                        eQ::parameters["MODULUS_IPTG"] = thisC;
//                        const double iptgOffset = 0.0;
//                        const double iptgOffset = 10.0;
//                        const double iptgOffset = 20.0;
//                        const double iptgOffset = 30.0;
//                        const double iptgOffset = 40.0;
//                        const double iptgOffset = 50.0;
                        const double iptgOffset = 60.0;
                        double iptgMin = 0.0;
                        double iptgMax = 30.0;
//                        double iptgMax = 00.0;//control input
//                        const double iptgMin = 20.0;
//                        const double iptgMax = 80.0;
//                        const double iptgMin = 10.0;
//                        const double iptgMax = 100.0;
//                        const double iptgMin = 100.0;
//                        const double iptgMax = 10.0;
//                        const double iptgMin = 30.0;
//                        const double iptgMax = 30.0;
                        iptgMin += iptgOffset;
                        iptgMax += iptgOffset;
                        double iptg = (x/trapWidthMicrons)*(iptgMax-iptgMin) + iptgMin;
                        eQ::parameters["MODULUS_IPTG"] = iptg;
                    }
                    std::vector<double> hslData;
                    std::vector<double> membraneDiffusion;
                    //HSL READ:
                    for(size_t i(0); i< Params.hslSolutionVector.size(); i++)
                    {
                        if(bool(eQ::parameters["PETSC_SIMULATION"]))
                            hslData.push_back(readHSL(Params.petscSolutionVector[i], Params.petscLookupTable[i], cellPoints));
                        else
                            hslData.push_back(readHSL(Params.hslSolutionVector[i], Params.dofLookupTable[i], cellPoints));
                    }

                    //TODO: move these to main where they are defined with the diffusion coeff. of each HSL
                    membraneDiffusion.push_back(3.0);

                    auto deltaHSL = thisCell->strain->computeProteins(hslData, membraneDiffusion, cellLength);

                    //HSL WRITE:
                    if(bool(eQ::parameters["PETSC_SIMULATION"]))
                        writeHSL(deltaHSL[0], Params.petscSolutionVector[0], Params.petscLookupTable[0], cellPoints);
                    else
                        writeHSL(deltaHSL[0], Params.hslSolutionVector[0], Params.dofLookupTable[0], cellPoints);

                    setDiffusionTensor(thisCell->getAngle(), cellPoints);
                if ("MODULUS_2" == eQ::parameters["simType"])
                {//HSL WRITE:
//                    writeHSL(deltaHSL[1], Params.hslSolutionVector[1], Params.dofLookupTable[1], cellPoints);
                }
            }
        }

        recordCelldata();

        //update the spring rest-length extension
        //note: possibly will be altered by genetic circuit, pressure, [protein], etc...
        //should pass growth rate, which should be computed in the ABM class
        (*cell)->updateGrowth();
    }

    if(bool(eQ::parameters["PETSC_SIMULATION"]))
    {
        //COPY FENICS I/O TO PETSC VECTOR FOR TESTING:
        for(size_t grid(0); grid< Params.hslSolutionVector.size(); grid++)
        {
            for (size_t i(0); i < (nodesHigh); i++)
            {
                for (size_t j(0); j < (nodesWide); j++)
                {
                    auto f = Params.dofLookupTable[grid]->grid[i][j];
                    auto p = Params.petscLookupTable[grid]->grid[i][j];
//                    Params.petscSolutionVector[grid].at(p) =  Params.hslSolutionVector[grid].at(f);
                    Params.hslSolutionVector[grid].at(f) =  Params.petscSolutionVector[grid].at(p);
                }
            }
        }
    }

    if(cellCounter != 0)  averagePointsPerCell /= double(cellCounter);
}
void eQabm::finalizeDataRecording(std::string fpath)
{//done at the end of the simulation for exporting ABM acquired data

    //data recording for debugging;  disable here by direct return:
    return;

    std::stringstream sstream;
    std::ofstream logFile;


    size_t maxIndex = 0;
    size_t maxSize = 0;
    for(size_t i(0); i<highResolutionDataVector.size(); i++)
    {
        if(highResolutionDataVector[i].size() > maxSize)
        {
            maxIndex = i;
            maxSize = highResolutionDataVector[i].size();
        }
    }
    std::cout<<"highResolutionDataVector maxIndex: "<< maxIndex << " size: "<<maxSize<<std::endl;
    sstream.str("");//to clear the string, not "clear()" which is only error flags
    int count=0;

    for(auto &cell : highResolutionDataVector)
    {
        sstream.str("");//to clear the string, not "clear()" which is only error flags
        for(auto &data : cell)
        {
                sstream <<
                           std::get<0>(data) << "    " <<
                           std::get<1>(data) << "    " <<
                           std::get<2>(data) << "    " <<
                           std::get<3>(data);
//            sstream << data.first << "    " << data.second;
            sstream<<std::endl;
        }
        logFile.open(fpath + "highRes_" + std::to_string(count++) + ".txt", std::ios::trunc);
        logFile << sstream.str();
        logFile.flush();
        logFile.close();
    }

//    return;


    for(auto comp : compressionTimeSeries)
        sstream << comp << std::endl;

    logFile.open(fpath + "_compEnergy.txt", std::ios::trunc);
    logFile << sstream.str();
    logFile.flush();
    logFile.close();
    std::cout<<"Compression Energy data written to file..."<<std::endl;

    sstream.str("");//to clear the string, not "clear()" which is only error flags
    for(auto angleData : angleBinTimeSeries)
    {
        for(auto bin : angleData)
            sstream << bin << ",";
        sstream<<std::endl;
    }
    logFile.open(fpath + "_binAngleData.txt", std::ios::trunc);
    logFile << sstream.str();
    logFile.flush();
    logFile.close();
        sstream.str("");//to clear the string, not "clear()" which is only error flags
        for(auto angleData : angleBinTimeSeries2)
        {
            for(auto bin : angleData)
                sstream << bin << ",";
            sstream<<std::endl;
        }
        logFile.open(fpath + "_binAngleData2.txt", std::ios::trunc);
        logFile << sstream.str();
        logFile.flush();
        logFile.close();
    std::cout<<"binAngleData written to file..."<<std::endl;

}
eQabm::~eQabm()
{
//    for(auto &cell : cellList)
//    {
////        std::cout<<"cell length: "<<cell->cpmCell->springRestLength<<std::endl;
//        cell.reset();
//    }
    cellPointers.reset();
    cellList.clear();
    trap.reset();
    habitat.reset();
}
void eQabm::recordDivisionEvent(double timeStamp, std::shared_ptr<eColi> cell, std::shared_ptr<eColi> daughter)
{
    if(nullptr == cell)//initial cell at t=0: set parent to -ID (length 0)
    {
        divisionList.push_back(
                    std::make_tuple(
                            timeStamp,
                            daughter->parentID, 0.0,//parent length=0
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
std::ostream &operator<<(std::ostream &os,  const eQ::strainType &type)
{
    if(eQ::strainType::ACTIVATOR == type)
        os<<"ACTIVATOR";
    else if(eQ::strainType::REPRESSOR == type)
        os<<"REPRESSOR";
    else if(eQ::strainType::X == type)
        os<<"X";
    else if(eQ::strainType::Y == type)
        os<<"Y";
    else if(eQ::strainType::Z == type)
        os<<"Z";
    else os<< "";
    return os;
}

/*
void eQabm::stepCells(fli_t begin, fli_t end,
                      eQ::gridFunction<size_t> &c4lookup, eQ::gridFunction<size_t> &c14lookup,
                      std::vector<double> &c4, std::vector<double> &c14)
{
//    std::shared_ptr<eQ::gridFunction<std::shared_ptr<Strain>>>
//            AGrid = sim->dsoGrid[0];
//    std::shared_ptr<eQ::gridFunction<std::shared_ptr<Strain>>>
//            RGrid = sim->dsoGrid[1];
    size_t dofC4, dofC14;
    double x,y;
    size_t i,j;
    std::pair<size_t, size_t> index;

    for(auto cell(begin); cell!=end; cell++)
    {
        x = (*cell)->getCenter_x();
        y = (*cell)->getCenter_y();
        index = eQ::scalarDataSource::ij_from_xy(x,y);
        i = index.first;  j = index.second;
        dofC4 = c4lookup.grid[i][j];
        dofC14 = c14lookup.grid[i][j];
            double &c4val = c4[dofC4];
            double &c14val = c14[dofC14];
            std::vector<double &> proteins = {c4val, c14val};
            (*cell)->updateProteins(proteins);
//            AGrid->grid[i][j]->computeProteins(c4, c14, 1.0);
//            RGrid->grid[i][j]->computeProteins(c4, c14, 1.0);
    }
}
*/
