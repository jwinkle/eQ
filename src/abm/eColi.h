#ifndef ECOLI_H
#define ECOLI_H

//jMODS:
#include <iostream>
#include <fstream>
#include <thread>
#include <vector>

#include "../eQ.h"
#include "cpmEColi.h"
#include "../Strain.h"


//#define PIXELS_PER_MICRON 10.0
#define DEFAULT_CELL_WIDTH_MICRONS              1.0
#define DEFAULT_DIVISION_LENGTH_MICRONS         4.2
//#define DEFAULT_DIVISION_LENGTH_MICRONS         4.
#define DEFAULT_CELL_DOUBLING_PERIOD_MINUTES    20.0

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#define KSPRING_SECS        1.0                     //normalize to 1.0 in SI units ([kg] = 1)
#define GAMMA_FLUID_SECS    1.0                     //normalize to 1.0 ""
#define KSPRING        (KSPRING_SECS*3600.0)    //convert to min^-2
#define GAMMA_FLUID    (GAMMA_FLUID_SECS*60.0) //convert to min^-1

#define DEFAULT_KSPRING_MINS        (KSPRING_SECS*3600.0)    //convert to min^-2
#define DEFAULT_GAMMA_FLUID_MINS    (KSPRING_SECS*60.0)

// #define GAMMA_FLUID         (SCALE_DT * 0.25e0 * KSPRING_SECS*60.0) //dt=0.001
//#define GAMMA_FLUID         (SCALE_DT  * KSPRING_SECS*60.0) //dt=0.001
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

class eColi  {
public:

    struct params
    {
//        struct eQ::simData  simData;
        cpSpace     *space;
        double mass, moment, length, width;
        double x, y, angle, vx, vy;
        double dRL;
        double divisionLength;
        double meanDivisionLength;
        double doublingPeriodMinutes;
        double stabilityScaling;
        //these should be coupled to a class:
        eQ::strainType                  strainType;
        struct dso_parameters           *dsoParams;
        std::shared_ptr<Strain>         strain;

    };


    eColi(const eColi::params &p);
    ~eColi();
    eColi::params Params;

//std::shared_ptr<eColi>      updatePhysicsModel(std::vector<double> &randomNumbers);
std::shared_ptr<eColi>      updatePhysicsModel(eQ::uniformRandomNumber  &zeroOne);
    void                        updateProteins(std::vector<double&> &protein);
    void                        updateGrowth();
//    void                    finalizeCell(std::vector<std::vector<std::pair<double, double>>> &data);
    void                    finalizeCell(std::vector<std::vector<std::tuple<double, double, double, double>>> &data);

    typedef  std::tuple<double, double, double, double> highResData_t;
//    std::vector<std::pair<double, double>> highResolutionData;
    std::vector<highResData_t> highResolutionData;

//    std::shared_ptr<eColi>      checkDivision(std::vector<double> &randomNumbers);
//    std::shared_ptr<eColi>      checkDivision(eQ::uniformRandomNumber  &zeroOne);
    void       setCellID ( long i ) { id = i; }
    long       getCellID ( void ) { return long(id); }
    double     getSpringCompression(){ return double(cpmCell->compression);}//length in microns
    double     getLengthMicrons(){ return double(cpmCell->length);}//length in microns
    double     getCenter_x ( void ) { return double(cpmCell->center.x);}
    double     getCenter_y ( void ) { return double(cpmCell->center.y);}
    double     getAngle(void) { return double(cpmCell->angle);}
    double     getExpansionSpeed(void) { return double(cpmCell->separationSpeed);}
    std::vector<double> getCellData(void)
    {
        std::vector<double> data = {getCenter_x(), getCenter_y(), getAngle(), getLengthMicrons(), getSpringCompression()};
        return data;
    }



    void updatePoleCenters(void);
    std::pair<double, double> polePositionA, polePositionB;


    std::shared_ptr<cpmEColi>       cpmCell;

    std::shared_ptr<Strain>         strain;

    long                            parentID;
//    std::vector<long>               daughterIDs;
private:
    void setDampingGamma();
    void computeDivisionLength(double);
    double computeDivisionLength(double lengthAtBirth, double meanLengthAtBirth, double alpha);
    double _dt;
    size_t id;

};

#endif
