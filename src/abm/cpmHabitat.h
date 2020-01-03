#ifndef CPM_HABITAT_H
#define CPM_HABITAT_H

#include "cpm.h"
#include "cpmTrap.h"
#include "cpmEColi.h"


class cpmHabitat
{
public:
	cpmHabitat(double sim_dt);
	~cpmHabitat();
//    jChipmunk_Habitat(double sim_dt,
//					  unsigned int numThreads,
//					  unsigned int numIterations);
//	unsigned long getNumThreads(void);

    void setTrap(std::shared_ptr<cpmTrap> trap){_trap = trap;}
    void stepSimulation(double dt);
    bool outsideTrap(std::shared_ptr<cpmEColi>);
    bool updateCell(std::shared_ptr<cpmEColi>);

    double getSimTime(){return simTime;}

private:
    std::shared_ptr<cpmTrap>			_trap;
    cpSpace         *_space;
    cpFloat         dt;
    double          simTime;

////allow cells to access the space pointer:
    friend class eQabm;
//    friend class jEColi;
//    friend class jChipmunk_Trap;
//    friend bool setBinData(jEColi * c, int which);

};

#endif // CPM_HABITAT_H
