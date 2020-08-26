#ifndef Ecoli_H
#define Ecoli_H

#include <iostream>
#include <fstream>
#include <thread>
#include <vector>

#include "../eQ.h"
#include "../eQrandom.h"
#include "../Strain.h"
#include "cpmEcoli.h"


class Ecoli : public eQ::Cell
{//Ecoli cell that uses Chipmunk2D physics engine for cell growth
public:
        static constexpr double  KSPRING_SECS               =       1.0;//normalize to 1.0 in SI units ([kg] = 1)
        static constexpr double  GAMMA_FLUID_SECS           =       1.0;//normalize to 1.0 ""
        static constexpr double  KSPRING                    =       (Ecoli::KSPRING_SECS*3600.0);    //convert to min^-2
        static constexpr double  GAMMA_FLUID                =       (Ecoli::GAMMA_FLUID_SECS*60.0); //convert to min^-1
        static constexpr double  DEFAULT_KSPRING_MINS       =       (Ecoli::KSPRING_SECS*3600.0);    //convert to min^-2
        static constexpr double  DEFAULT_GAMMA_FLUID_MINS   =       (Ecoli::KSPRING_SECS*60.0);

    struct Params : eQ::Cell::Params //public by default for structs
    {
        cpSpace                         *space;
        double                          stabilityScaling;
        std::shared_ptr<Strain>         strain;
    };

    Ecoli(const Ecoli::Params &p);
    ~Ecoli() override;
    Ecoli::Params params;

    std::shared_ptr<Ecoli>      updatePhysicsModel(eQ::uniformRandomNumber  &zeroOne);
    void                        updateGrowth();

    double     getLengthMicrons()   override {return double(cpmCell->length);}//length in microns
    double     getCenter_x()        override {return double(cpmCell->center.x);}
    double     getCenter_y()        override {return double(cpmCell->center.y);}
    double     getAngle()           override {return double(cpmCell->angle);}
    double     getWidthMicrons()    override {return params.width;}//length in microns


    using eQ::Cell::computeDivisionLength;
    void computeDivisionLength(double birthLength)
    {
        params.divisionLength =  eQ::Cell::computeDivisionLength(birthLength, 0.5*params.meanDivisionLength, params.divisionCorrelationAlpha);
    };

    double     getSpringCompression()   {return double(cpmCell->compression);}//length in microns
    double     getExpansionSpeed()      {return double(cpmCell->separationSpeed);}
    void       updatePoleCenters();

    std::pair<double, double>       polePositionA, polePositionB;

    std::shared_ptr<cpmEcoli>       cpmCell;
    std::shared_ptr<Strain>         strain;

private:
    void    setDampingGamma();

};

#endif
