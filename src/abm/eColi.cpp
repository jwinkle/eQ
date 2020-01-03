#include "eColi.h"


eColi::eColi(const eColi::params &p)
        : Params(p)
{
    //build the chipmunk cell model parameters:
    //TODO: define a cell geometry structure to copy directly to the physics model layer (m,moment,x,y,v,theta,w,l)
    cpmEColi::params cell;
    cell.space = Params.space;
        cell.mass = 1.0;  cell.moment = 100.0;
        cell.x = Params.x;  cell.y = Params.y;
        cell.angle = Params.angle;
        cell.length = Params.length;//in microns
        cell.width = DEFAULT_CELL_WIDTH_MICRONS;//over-ride width for now...
        cell.vx = Params.vx;  cell.vy = Params.vy;
        cell.av = 0.0;
        cell.dRL = Params.dRL;
        cell.kspring                = DEFAULT_KSPRING_MINS;

    //create the chipmunk class
    cpmCell = std::make_shared<cpmEColi>(cell);

    //26Feb.2019: scale damping by cell length to remove discontinuity across cell divisions
    setDampingGamma();
    updatePoleCenters();

    strain = Params.strain;


}
void eColi::computeDivisionLength(double birthLength)
{
    //15Feb.2019:  cell size regulation using ~"incremental model" by Ariel Amir (PRL 112, 2014)
    double alpha = 0.5;
    Params.divisionLength
            = computeDivisionLength(birthLength, 0.5*Params.meanDivisionLength, alpha);
}
double eColi::computeDivisionLength(double lengthAtBirth, double meanLengthAtBirth, double alpha)
{
        //15Feb.2019:  cell size regulation using ~"incremental model" by Ariel Amir (PRL 112, 2014)
        return 2.0 * pow(lengthAtBirth, 1.0-alpha) * pow(meanLengthAtBirth, alpha);
    //"sizer" regulation, alpha=1
//    return 2.0*meanLengthAtBirth;
}
void eColi::setDampingGamma()
{
//    double gammaScale = 1.0;
//    double gammaScale = 0.5;
    double gammaScale = 0.25;
//    double gammaScale = 0.1;
//    double gammaScale = 0.05;
    cpmCell->parameterData.gammaFluid
            = gammaScale * DEFAULT_GAMMA_FLUID_MINS * Params.stabilityScaling * getLengthMicrons();
//            = DEFAULT_GAMMA_FLUID_MINS * Params.stabilityScaling;
}
void eColi::updateProteins(std::vector<double &> &protein)
{

}
//void eColi::finalizeCell(std::vector<std::vector<std::pair<double, double>>> &data)
void eColi::finalizeCell(std::vector<std::vector<std::tuple<double, double, double, double>>> &data)
{
//    const double lifeTimeThreshold = 100.0 * (1.0/double(eQ::parameters["dt"]));
//    if(highResolutionData.size() > lifeTimeThreshold)
//        data.push_back(highResolutionData);

}
void eColi::updatePoleCenters(void)
{
    double trapWidth = double(eQ::parameters["physicalTrapWidth_X_Microns"]);
    double trapHeight = double(eQ::parameters["simulationTrapHeightMicrons"]);

    cpVect polePosition;
    int r = 1;
    double a = cpmCell->angle;
    polePosition = cpmCell->center + cpvmult ( cpv ( cos ( a ), sin ( a ) ),
                                               (r)*0.5*getLengthMicrons() -  DEFAULT_CELL_WIDTH_MICRONS/2.0);

    if(polePosition.x < 0.0) polePosition.x = 0.0;
    if(polePosition.y < 0.0) polePosition.y = 0.0;
    if(polePosition.x > trapWidth) polePosition.x = trapWidth;
    if(polePosition.y > trapHeight) polePosition.y = trapHeight;

    polePositionA =    std::make_pair(double(polePosition.x), double(polePosition.y));

    //change to  '-r':
    polePosition = cpmCell->center + cpvmult ( cpv ( cos ( a ), sin ( a ) ),
                                               (-r)*0.5*getLengthMicrons() -  DEFAULT_CELL_WIDTH_MICRONS/2.0);
    if(polePosition.x < 0.0) polePosition.x = 0.0;
    if(polePosition.y < 0.0) polePosition.y = 0.0;
    if(polePosition.x > trapWidth) polePosition.x = trapWidth;
    if(polePosition.y > trapHeight) polePosition.y = trapHeight;

    polePositionB =    std::make_pair(double(polePosition.x), double(polePosition.y));
}
void eColi::updateGrowth()
{
    //call the chipmunk model
    int status;
    //set the growth force, with (generally) updated dRL rest-length expansion
    //Exponential growth: l(t) = l0 * 2^(t/Td) ==> dl = dt * l * 2^(dt/20) * ln(2)/20
    //NOTE: 2^(dt/20) =~ 1 so skip that in the calculation
    const double expScale = (log(2.0)/Params.doublingPeriodMinutes) * double(eQ::parameters["dt"]);
//    status = cpmCell->updateExpansionForce(Params.dRL);
    status = cpmCell->updateExpansionForce(expScale * getLengthMicrons());
}

//std::shared_ptr<eColi> eColi::updatePhysicsModel(std::vector<double> &randomNumbers)
std::shared_ptr<eColi> eColi::updatePhysicsModel(eQ::uniformRandomNumber  &zeroOne)
{
      //call the chipmunk model
      int status;
      status = cpmCell->updateModel();
      setDampingGamma();
      updatePoleCenters();

      if( ("MODULUS_1" == eQ::parameters["simType"]) ||
          ("MODULUS_2" == eQ::parameters["simType"]))
      {
//          auto pmodulus = std::dynamic_pointer_cast<MODULUSmodule>(strain);
//          auto thisValue = pmodulus->computeAverageCellLength();

//          auto thisData =  std::make_tuple(
//                      getLengthMicrons(),
//                      strain->getProteinNumber(FP),
//                      getExpansionSpeed(),
//                      getSpringCompression());
//          highResolutionData.push_back(thisData);
      }


//=================================================================================================
      //CHECK DIVISION:
//=================================================================================================
    double parentLength = getLengthMicrons();
    if(parentLength < Params.divisionLength)
        return nullptr;

    //random assignment of cell ID to cell pole:
    // int r = frand() > 0.5 ? 1 : -1;
    //fixed assignment:
    int r = -1;

    double dn = double(eQ::parameters["divisionNoiseScale"]);
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

    //eventually read these from the parent cell:
    cpmEColi::params cell;//need copy constructor
        cell.space  = Params.space;
        cell.mass   = cpmCell->initialMass;
        cell.moment = cpmCell->initialMoment;
        cell.x      = newpos.x;
        cell.y      = newpos.y;
        cell.angle  = a - r*da;
        cell.length = newLength;
        cell.width  = DEFAULT_CELL_WIDTH_MICRONS;
        cell.vx     = vel.x;
        cell.vy     = vel.y;
        cell.dRL    = cpmCell->dRL;
        cell.av     = 0.0;
//        cell.gammaFluidParameter    = DEFAULT_GAMMA_FLUID_MINS * Params.stabilityScaling;
        cell.gammaFluidParameter    = DEFAULT_GAMMA_FLUID_MINS * Params.stabilityScaling * cell.length;
        cell.kspring                = DEFAULT_KSPRING_MINS;


    auto newCell = std::make_shared<cpmEColi>(cell);

    //need to use averaged values to avoid noise in daughter cells:
    double newRestLength = (cpmCell->springRestLength - cpmCell->length) + cell.length;
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

    eColi::params dparams = Params;//copy parameters to new cell (sets all defaults)
            //update parameters different for daughter cell
            dparams.x = oldpos.x + r*0.5*parentLength*frac*cos ( a + r*da );
            dparams.y = oldpos.y + r*0.5*parentLength*frac*sin ( a + r*da );
            dparams.angle = a + r*da;
            dparams.length = dLength;
            dparams.vx = newV2.x;
            dparams.vy = newV2.y;
            dparams.dRL = Params.dRL;//nominal growth rate is 0.1 um/min.
            //placeholders here, in case should change for daughter on division:
            dparams.stabilityScaling = Params.stabilityScaling;
            dparams.width = DEFAULT_CELL_WIDTH_MICRONS;
//            dparams.divisionLength = DEFAULT_DIVISION_LENGTH_MICRONS;//"sizer" regulation if set constant
            dparams.meanDivisionLength = Params.meanDivisionLength;//propagate the div. length
            dparams.doublingPeriodMinutes = Params.doublingPeriodMinutes;//propagate the growth rate.


    dparams.strain = strain->clone();

    auto daughter = std::make_shared<eColi>(dparams);

    daughter->computeDivisionLength(dLength);
    daughter->cpmCell->setVelocity(newV2);//set explicitly to cpm model
    daughter->cpmCell->setSpringRestLength(newRestLength - cell.length + dLength);
    daughter->setDampingGamma();
    daughter->updatePoleCenters();

    //divide proteins to daughter cell:
    //will call over-ride function if defined (note: passes base class pointer as daughter strain)
    strain->divideProteins(daughter->strain, fracToDaughter);

    return daughter;

}

eColi::~eColi()
{
    cpmCell.reset();
}
//  protein = 0.0;//will be reset later
//  // plasmidCopyNumber = 40;
//  plasmidCopyNumber = 20;//expected # plasmids for daughter cells (just after division)
//   // plasmidCopyNumber = 10;
//  plasmidCount = plasmidCopyNumber;//will be possibly reset later

//  pAlpha = PROTEIN_ALPHA;
//  pBeta = PROTEIN_BETA;
//  pDil = PROTEIN_DILUTION;
//  xFPmax = (16.0*pAlpha)/(pDil+pBeta);//max protein for binning protein data

//  lineage = new lineageTracker(this);

//  // Strain::type thisCircuitType = (ECOLITYPE_01 == ) ? Strain::type::ACTIVATOR : Strain::type::REPRESSOR;

//  circuit = new Strain(cellType, NULL, dt);

//void jEColi::setPlasmidCopyNumber(double rate)
//{
//  //only called on init of space (not on divide/daughter)
//  plasmidCopyNumber = rate;
//  plasmidCount = rate;
//}

//void jEColi::jupdate(double &C4ext, double &C14ext)
//{
//  //exclude w/ non HSL sims:
//  // circuit->computeProteins(C4ext, C14ext, getLengthMicrons());

//    // if(Strain::type::ACTIVATOR == cellType)
//    // {
//    //   double ftsZ = circuit->conc[L]/getLengthMicrons();
//    //   if (ftsZ > 10.0)
//    //     div_vol = DIV_VOL*(2.0/3.0);//scale
//    //   else
//    //     div_vol = DIV_VOL;
//    // }

//  if(circuit->inductionActive())
//    div_vol = DIV_VOL*(2.0/3.0);//scale
//  else
//    div_vol = DIV_VOL;
//}
//void jEColi::jupdate(double &C4ext)
//{
//  // int status = jupdate();
//  if(circuit->inductionActive())
//  {
//    if(Strain::type::REPRESSOR == cellType)
//    {
//      double ftsZ = circuit->conc[L]/getLengthMicrons();
//      if (ftsZ > 10.0)
//        div_vol = DIV_VOL*(2.0/3.0);//scale
//      else
//        div_vol = DIV_VOL;
//    }
//  }
//  double devnull;
//  circuit->computeProteins(C4ext, devnull, getLengthMicrons());
//}


//int jEColi::jupdate( size_t index )
//{

//  #define SCALE_NUMS 13
//  extern bool g_inductionFlag;

//  const double scale0A[SCALE_NUMS] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9, 0.7, 0.5, 0.5, 0.5, 0.5};
//  const double scale0B[SCALE_NUMS] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9, 0.7, 0.5};

//  const double scale1A[SCALE_NUMS] = {1.0, 0.9, 0.7, 0.5, 0.5, 0.5, 0.5, 0.9, 0.7, 0.5, 0.5, 0.5, 0.5};
//  const double scale1B[SCALE_NUMS] = {1.0, 1.0, 1.0, 1.0, 0.9, 0.7, 0.5, 1.0, 1.0, 1.0, 0.9, 0.7, 0.5};

//  div_vol = DIV_VOL;
//  if (index >= SCALE_NUMS)//error over-range index
//    return jupdate();

//  if(g_inductionFlag)
//  {
//    div_vol  *= (Strain::type::ACTIVATOR == cellType) ? scale1A[index] : scale1B[index];
//  }
//  else
//  {
//    div_vol  *= (Strain::type::ACTIVATOR == cellType) ? scale0A[index] : scale0B[index];
//  }
