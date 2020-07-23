#ifndef CPM_HABITAT_H
#define CPM_HABITAT_H

#include "cpm.h"
#include "cpmTrap.h"
#include "cpmEcoli.h"


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
    bool outsideTrap(std::shared_ptr<cpmEcoli>);
    bool updateCell(std::shared_ptr<cpmEcoli>);

    double getSimTime(){return simTime;}

    cpSpace *get_cpSpace(){return _space;}
private:
    std::shared_ptr<cpmTrap>			_trap;
    cpSpace         *_space;
    cpFloat         dt;
    double          simTime;

////allow cells to access the space pointer:
//    friend class eQabm;
//    friend class jEColi;
//    friend class jChipmunk_Trap;
//    friend bool setBinData(jEColi * c, int which);

};

#endif // CPM_HABITAT_H
