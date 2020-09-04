#include "Ecoli.h"


Ecoli::Ecoli(const Ecoli::Params &p)
        : Cell(), params(p)
{
    //build the chipmunk cell model parameters:
    cpmEcoli::Params cell;

    cell.space                  = params.space;
    cell.baseData               = params;//will implicitly cast to eQ::Cell::Params to only copy that part
    cell.kspring                = Ecoli::DEFAULT_KSPRING_MINS;

    //create the chipmunk class
    cpmCell = std::make_shared<cpmEcoli>(cell);

    //26Feb.2019: scale damping by cell length to remove discontinuity across cell divisions
    setDampingGamma();
    updatePoleCenters();

    //note:  strain object has already been created
    strain = params.strain;
    strain->params.baseData = &params;//allows the Strain class to access the eQ::Cell::Params by reference
    strain->init();
}
void Ecoli::setDampingGamma()
{
//    double gammaScale = 1.0;
//    double gammaScale = 0.5;
    double gammaScale = 0.25;
//    double gammaScale = 0.1;
//    double gammaScale = 0.05;
    cpmCell->parameterData.gammaFluid
            = gammaScale * Ecoli::DEFAULT_GAMMA_FLUID_MINS * params.stabilityScaling * getLengthMicrons();
}
void Ecoli::updatePoleCenters(void)
{
    double trapWidth = double(eQ::data::parameters["simulationTrapWidthMicrons"]);
    double trapHeight = double(eQ::data::parameters["simulationTrapHeightMicrons"]);

    cpVect polePosition;
    int r = 1;
    double a = cpmCell->angle;
    polePosition = cpmCell->center + cpvmult ( cpv ( cos ( a ), sin ( a ) ),
                                               (r)*0.5*getLengthMicrons() -  getWidthMicrons()/2.0);

    if(polePosition.x < 0.0) polePosition.x = 0.0;
    if(polePosition.y < 0.0) polePosition.y = 0.0;
    if(polePosition.x > trapWidth) polePosition.x = trapWidth;
    if(polePosition.y > trapHeight) polePosition.y = trapHeight;

    polePositionA =    std::make_pair(double(polePosition.x), double(polePosition.y));

    //change to  '-r':
    polePosition = cpmCell->center + cpvmult ( cpv ( cos ( a ), sin ( a ) ),
                                               (-r)*0.5*getLengthMicrons() -  getWidthMicrons()/2.0);
    if(polePosition.x < 0.0) polePosition.x = 0.0;
    if(polePosition.y < 0.0) polePosition.y = 0.0;
    if(polePosition.x > trapWidth) polePosition.x = trapWidth;
    if(polePosition.y > trapHeight) polePosition.y = trapHeight;

    polePositionB =    std::make_pair(double(polePosition.x), double(polePosition.y));
}
void Ecoli::updateGrowth()
{
    //allow the strain to scale the growth rate ( \in [0,1] )
    double growthScale = strain->growthRateScaling();

    //call the chipmunk model
    int status;
    //set the growth force, with (generally) updated rest-length expansion rate
    //Exponential growth: l(t) = l0 * 2^(t/Td) ==> dl = dt * l * 2^(dt/20) * ln(2)/20
    //NOTE: 2^(dt/20) =~ 1 so skip that in the calculation
    const double expScale = (log(2.0)/params.doublingPeriodMinutes) * double(eQ::data::parameters["dt"]);
    status = cpmCell->updateExpansionForce(growthScale * expScale * getLengthMicrons());
}

std::shared_ptr<Ecoli> Ecoli::updatePhysicsModel(eQ::uniformRandomNumber  &zeroOne)
{
      //call the chipmunk model
      int status;
      status = cpmCell->updateModel();

      //update the parameters structure (manually for now)
      params.x         = getCenter_x();
      params.y         = getCenter_y();
      params.length    = getLengthMicrons();
      params.angle     = getAngle();
      //velocity if needed

      setDampingGamma();
      updatePoleCenters();

//=================================================================================================
      //CHECK DIVISION:
//=================================================================================================
    double parentLength = getLengthMicrons();
    if(parentLength < params.divisionLength)
        return nullptr;


    int r = -1;

    double dn = double(eQ::data::parameters["divisionNoiseScale"]);
    double frac = 0.5 + dn * ( zeroOne.randomNumber() - 0.5 );
    double fracToDaughter = 1.0-frac;

    //============================================================
    //      CREATE NEW PARAMETERS FOR PARENT CELL
    //============================================================

    double newLength = frac * parentLength;
    double a = cpmCell->angle;
    double da = 0.0;//angle deviation;  set to zero to avoid spurious oscillations!

    cpVect oldpos = cpmCell->center;
    cpVect newpos = cpmCell->center + cpvmult ( cpv ( cos ( a - r*da ), sin ( a - r*da ) ),
                                               (-r)*0.5*parentLength*(1-frac) );
    //assign the old cell-half velocities to new center-of-mass velocities
    //NOTE:  cell will not be expanding at first subsequent time step
    cpVect vel = (1 == r) ?  cpmCell->velB:cpmCell->velA;

    cpmEcoli::Params cell;
        cell.space      = params.space;
        cell.baseData   = params;

        cell.baseData.x      = newpos.x;
        cell.baseData.y      = newpos.y;
        cell.baseData.angle  = a - r*da;
        cell.baseData.length = newLength;
        cell.baseData.vx     = vel.x;
        cell.baseData.vy     = vel.y;

        cell.kspring         = cpmCell->k_spring;

    auto newCell = std::make_shared<cpmEcoli>(cell);

    //need to use averaged values to avoid noise in daughter cells:
    double newRestLength = (cpmCell->springRestLength - cpmCell->length) + cell.baseData.length;
    newCell->setSpringRestLength(newRestLength);

    //determine whether it is the "right" or "left" cell;  r=1 is "right"=bodyB
    cpVect newV1, newV2;
    if(1 == r)
    {
        newV1 = cpmCell->velB;
        newV2 = cpmCell->velA;
    }
    else
    {
        newV1 = cpmCell->velA;
        newV2 = cpmCell->velB;
    }

    newCell->setVelocity(newV1);

    //replace the old cell with the new one:
    cpmCell.reset();
    cpmCell = newCell;

    computeDivisionLength(newLength);
    setDampingGamma();
    updatePoleCenters();
  //============================================================
  //  new daughter cell
  //============================================================
    double dLength = fracToDaughter*parentLength;

    Ecoli::Params dparams = params;//copy parameters to new cell (sets all defaults)
        //update parameters different for daughter cell
        dparams.x = oldpos.x + r*0.5*parentLength*frac*cos ( a + r*da );
        dparams.y = oldpos.y + r*0.5*parentLength*frac*sin ( a + r*da );
        dparams.angle = a + r*da;
        dparams.length = dLength;
        dparams.vx = newV2.x;
        dparams.vy = newV2.y;

    dparams.strain = strain->clone();

    auto daughter = std::make_shared<Ecoli>(dparams);

    daughter->computeDivisionLength(dLength);
    daughter->cpmCell->setVelocity(newV2);//set explicitly to cpm model
    daughter->cpmCell->setSpringRestLength(newRestLength - cell.baseData.length + dLength);
    daughter->setDampingGamma();
    daughter->updatePoleCenters();

    //divide proteins to daughter cell:
    //will call over-ride function if defined (note: passes base class pointer as daughter strain)
    strain->divideProteins(daughter->strain, fracToDaughter);

    return daughter;

}

Ecoli::~Ecoli()
{
    cpmCell.reset();
}
