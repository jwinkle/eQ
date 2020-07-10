#ifndef ECOLI_H
#define ECOLI_H

#include <iostream>
#include <fstream>
#include <thread>
#include <vector>

#include "../eQ.h"
#include "../eQrandom.h"
#include "../Strain.h"
#include "cpmEColi.h"


class eColi : public eQ::Cell
{//Ecoli cell that uses Chipmunk2D physics engine for cell growth
public:
        static constexpr double  KSPRING_SECS               =       1.0;//normalize to 1.0 in SI units ([kg] = 1)
        static constexpr double  GAMMA_FLUID_SECS           =       1.0;//normalize to 1.0 ""
        static constexpr double  KSPRING                    =       (eColi::KSPRING_SECS*3600.0);    //convert to min^-2
        static constexpr double  GAMMA_FLUID                =       (eColi::GAMMA_FLUID_SECS*60.0); //convert to min^-1
        static constexpr double  DEFAULT_KSPRING_MINS       =       (eColi::KSPRING_SECS*3600.0);    //convert to min^-2
        static constexpr double  DEFAULT_GAMMA_FLUID_MINS   =       (eColi::KSPRING_SECS*60.0);

    struct Params
    {
        Cell::Params                    baseData;//updates each time step (check if all variables are updated)
        cpSpace                         *space;
        double                          stabilityScaling;
        std::shared_ptr<Strain>         strain;
//        dso_parameters                  *dsoParams;
    };

    eColi(const eColi::Params &p);
    ~eColi() override;
    eColi::Params params;

    std::shared_ptr<eColi> clone(std::shared_ptr<Strain> &newStrain);

    std::shared_ptr<eColi>      updatePhysicsModel(eQ::uniformRandomNumber  &zeroOne);
    void                        updateGrowth();

    double     getLengthMicrons()   override {return double(cpmCell->length);}//length in microns
    double     getCenter_x()        override {return double(cpmCell->center.x);}
    double     getCenter_y()        override {return double(cpmCell->center.y);}
    double     getAngle()           override {return double(cpmCell->angle);}
    double     getWidthMicrons()    override {return params.baseData.width;}//length in microns

    double     getSpringCompression()   {return double(cpmCell->compression);}//length in microns
    double     getExpansionSpeed()      {return double(cpmCell->separationSpeed);}
    void       updatePoleCenters();

    std::pair<double, double>       polePositionA, polePositionB;

    std::shared_ptr<cpmEColi>       cpmCell;
    std::shared_ptr<Strain>         strain;

private:
    void    setDampingGamma();

};

#endif
