#ifndef CPM_ECOLI_H
#define CPM_ECOLI_H

#include "cpm.h"
#include <utility>
#include "eQ.h"
#include "eQcell.h"


class eColi;
typedef struct{
	int which;
	eColi *parent;
} body_t;

//too low elast. results in overlap?:
#define SHAPE_ELASTICITY     0.1
//#define SHAPE_ELASTICITY     1.0f
#define SHAPE_FRICTION     0.0
//#define SHAPE_FRICTION     0.8f

#define GROOVE_JOINT_MAXFORCE   1.0e6
//#define GROOVE_JOINT_MAXFORCE   5000.0
//#define GROOVE_JOINT_MAXFORCE   2500.0
//#define GROOVE_JOINT_MAXFORCE   1024.0
#define GROOVE_JOINT_ERRORBIAS   0.5
#define GROOVE_JOINT_MAXBIAS    40.0

class cpmEColi
{
public:
    struct Params
	{
        eQ::Cell::Params        baseData;
        cpSpace                 *space;
        double                  kspring;
        double                  gammaFluidParameter;
	};

    cpmEColi(const cpmEColi::Params &p);
	~cpmEColi(void);

    typedef struct bodyData_t
    {
//        double gammaFluidInverse;
        double gammaFluid;
    }bodyData_t;

    bodyData_t parameterData;

	int     updateModel(void);
    void    setVelocity(cpVect vel);
    void    setBodyVelocities(cpVect vA, cpVect vB);
    void    applyForce(cpVect);
	double   getScaledExpansionForce(double scale);
    int     updateExpansionForce(double);
	void    setSpringRestLength(double);
	double getSpringRestLength(void);
	void    incSpringRestLength(double deltaRL);

// private:
	void    applyGrowthForce(double growthRate, bool expanding);
    void    setStaticCell();
    void    resetDynamicCell();
    bool    cellIsStatic;
    bool    pointIsInCell(std::pair<double, double> point);
    void    computeCellBoudaryEdges();

    double              k_spring;

    cpVect              edges[4];

    cpVect              vertsA[4];
    cpVect              vertsB[4];
    cpFloat             springRestLength;
    cpFloat             dRL;
    cpFloat             compression;

    cpVect               posA, posB;
    cpVect               velA, velB;
    cpFloat              separationDistance, separationAngle, separationVelocityAngle;
    cpFloat              separationSpeed, projectionSpeed;
    bool                 cellExpanding;
    cpVect               center;
    cpFloat              angleA,angleB,angle;
    cpFloat              length;
    cpVect               velocity, separationVelocity;
    cpFloat              angularVelocity;
    cpVect              forceVE;
    body_t              pA,pB;

    cpFloat             radius, offset;
    cpFloat             nextRatchet;
    cpFloat             shapeLength;
    cpFloat             shapeHeight;
    cpFloat             initialLength;
    cpFloat             initialMass, initialMoment;

//    jEColi              *parent;
//    jChipmunk_Habitat   *habitat;
    cpSpace             *space;
    cpBody              *bodyA, *bodyB;

        cpShape           *shapes[32];
        unsigned int      shapeCount;
    cpShape             *capA, *capB;
    cpShape             *capA2, *capB2;
    cpShape             *boxA, *boxB;
        cpConstraint        *constraints[32];
        unsigned int      constraintCount;
    cpConstraint        *grooveJointA, *grooveJointB;
    cpConstraint        *grooveJointA2, *grooveJointB2;
    cpConstraint        *gearJoint;
    cpConstraint        *springJoint1, *springJoint2;

    cpFloat             cpSpringRestLength;


//    bool        firstIterationFlag;
//    cpFloat     vc, v0;

};


#endif // CPM_ECOLI_H
