#include "cpm.h"
#include "cpmHabitat.h"
#include "cpmTrap.h"

#define jITERATIONS 10
//#define jITERATIONS 20
// #define jITERATIONS 50

cpmHabitat::cpmHabitat(double sim_dt)
    : dt(sim_dt), _trap(nullptr)
{
    //default settings:
    _space = cpSpaceNew();
    cpSpaceSetCollisionSlop(_space, 0.1);
    cpSpaceSetCollisionBias(_space, 0.1);
    cpSpaceUseSpatialHash(_space, 10.0, 20000);
    cpSpaceSetIterations(_space, 10);
    cpSpaceSetDamping(_space, 1.0);//1.0 => no damping
}
//jChipmunk_Habitat::jChipmunk_Habitat(double sim_dt, unsigned int nthreads, unsigned int niterations)
//    : dt(sim_dt), simTime(0.0),  trap(nullptr), numThreads(nthreads)
//    {
//    //parent = w;
//    space = cpSpaceNew();
//    // space = cpHastySpaceNew();
//        // cpHastySpaceSetThreads(space, (unsigned long)nthreads);
//        // numThreads = cpHastySpaceGetThreads(space);

////    cpSpaceSetCollisionSlop(space, 0.01f);
//    cpSpaceSetCollisionSlop(space, 0.1);
//    //    cpSpaceSetCollisionSlop(space, 0.2f);
////        cpSpaceSetCollisionSlop(space, 0.5f);

//    cpSpaceSetCollisionBias(space, 0.1);
////    cpSpaceSetCollisionBias(space, 0.5);
////    cpSpaceSetCollisionPersistence(space, 6);//cpTimestamp value)

//    cpSpaceUseSpatialHash(space, 10.0, 20000);
////    cpSpaceUseSpatialHash(space, 15.0f, 10000);
//    // cpSpaceSetIterations(space, jITERATIONS);
//    cpSpaceSetIterations(space, int(niterations));

////a lower number for damping is more damping:
//    //  1.0 ==> no damping.
//    // epsilon->0.0 ==> full damping (must be > 0).
//    // cpSpaceSetDamping(space, 0.8);
////    cpSpaceSetDamping(space, 1.0e-6f);
//   cpSpaceSetDamping(space, 1.0);

//}

cpmHabitat::~cpmHabitat()
{
    cpSpaceFree(_space);
    // cpHastySpaceFree(space);
}
bool cpmHabitat::outsideTrap(std::shared_ptr<cpmEColi> cell)
{
    return (nullptr != _trap) ? _trap->outsideTrap(cell) : false;
}

bool cpmHabitat::updateCell(std::shared_ptr<cpmEColi> cell)
{
    return (nullptr != _trap) ? _trap->updateModel(cell) : false;
}

void cpmHabitat::stepSimulation(double dt)
{
    cpSpaceStep(_space, cpFloat(dt));
    simTime += dt;
    // cpHastySpaceStep(space, cpFloat(dt));
}
