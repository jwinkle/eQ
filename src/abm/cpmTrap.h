#ifndef CPM_TRAP_H
#define CPM_TRAP_H

#include "../eQ.h"
#include "cpm.h"
#include "cpmEColi.h"


class cpmTrap
{
public:
    struct Params
    {
        cpSpace     *space;
//        bool        removeOutsideTrap;
    };

    cpmTrap(const cpmTrap::Params &);
    ~cpmTrap();
    const cpmTrap::Params &myParams;

    bool outsideTrap(std::shared_ptr<cpmEColi> cell);
    bool updateModel(std::shared_ptr<cpmEColi> cell);
    

    double w;
    double h;

    cpBody      *staticBody;
    std::vector<cpShape *> trap;

};

#endif // JCHIPMUNK_TRAP_H
