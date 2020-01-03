#include "cpmTrap.h"
#include "cpmEColi.h"

//#define TRAP_ELASTICITY     1.0f
#define TRAP_ELASTICITY     0.0
#define TRAP_FRICTION       0.0

// #define TRAP_BORDER_WIDTH  2.0f
//#define TRAP_BORDER_WIDTH  10.0
#define TRAP_BORDER_WIDTH  1.0//converted to microns

//static  inline cpFloat frand(void) { return (cpFloat)rand()/(cpFloat)RAND_MAX; }
//static  inline cpFloat frand(void) { return cpFloat(rand()/RAND_MAX); }


cpmTrap::cpmTrap(const cpmTrap::params &p)
//cpmTrap::cpmTrap(const cpmTrap::params &p)
       : myParams(p)
{
    staticBody = cpSpaceGetStaticBody(myParams.space);

    if("NOWALLED" == eQ::parameters["trapType"])
        return;

    w = double(eQ::parameters["simulationTrapWidthMicrons"]);
    h = double(eQ::parameters["simulationTrapHeightMicrons"]);

    if("NOWALLED" != eQ::parameters["trapType"])
    {//1,2,or 3
        cpShape *shape;

        if("LEFTCORNER_HALF" == eQ::parameters["trapType"])
        {
            //LEFT VERTICAL WALL:
            shape = cpSpaceAddShape(
                        myParams.space, cpSegmentShapeNew(
                            staticBody, cpv(0.0, 0.0), cpv(0.0, h), TRAP_BORDER_WIDTH));
                trap.push_back(shape);
            //TOP 1/2 WALL:
            shape = cpSpaceAddShape(
                        myParams.space, cpSegmentShapeNew(
                            staticBody, cpv(0.0, h), cpv(w/2.0, h), TRAP_BORDER_WIDTH));
                trap.push_back(shape);
        }
        else if("OPPOSITE_CORNERS" == eQ::parameters["trapType"])
        {
            //LEFT VERTICAL WALL: TOP HALF
            shape = cpSpaceAddShape(
                        myParams.space, cpSegmentShapeNew(
                            staticBody, cpv(0.0, h/2.0), cpv(0.0, h), TRAP_BORDER_WIDTH));
                trap.push_back(shape);
            //TOP 1/4 WALL:
            shape = cpSpaceAddShape(
                        myParams.space, cpSegmentShapeNew(
                            staticBody, cpv(0.0, h), cpv(w/4.0, h), TRAP_BORDER_WIDTH));
                trap.push_back(shape);
            //RIGHT VERTICAL WALL:
            shape = cpSpaceAddShape(
                        myParams.space, cpSegmentShapeNew(
                            staticBody, cpv(w, 0.0), cpv(w, h), TRAP_BORDER_WIDTH));
                trap.push_back(shape);
//            //RIGHT VERTICAL WALL: BOTTOM HALF
//            shape = cpSpaceAddShape(
//                        myParams.space, cpSegmentShapeNew(
//                            staticBody, cpv(w, h/2.0), cpv(w,0.0), TRAP_BORDER_WIDTH));
//                trap.push_back(shape);
//            //BOTTOM 1/4 WALL:
//            shape = cpSpaceAddShape(
//                        myParams.space, cpSegmentShapeNew(
//                            staticBody, cpv(w,0.0), cpv(3.0*w/4.0, 0.0), TRAP_BORDER_WIDTH));
//                trap.push_back(shape);
        }
        else if("ONEWALLED_LEFT" == eQ::parameters["trapType"])
        {
            //LEFT VERTICAL WALL:
            shape = cpSpaceAddShape(
                        myParams.space, cpSegmentShapeNew(
                            staticBody, cpv(0.0, 0.0), cpv(0.0, h), TRAP_BORDER_WIDTH));
                trap.push_back(shape);
        }
        else if("H_TRAP" == eQ::parameters["trapType"])
        {
            //BOTTOM HORIZONTAL WALL:
            shape = cpSpaceAddShape(
                        myParams.space, cpSegmentShapeNew(
                            staticBody, cpv(0.0, 0.0), cpv(w, 0.0), TRAP_BORDER_WIDTH));
                trap.push_back(shape);
            //TOP HORIZONTAL WALL:
            shape = cpSpaceAddShape(
                        myParams.space, cpSegmentShapeNew(
                            staticBody, cpv(0.0, h), cpv(w, h), TRAP_BORDER_WIDTH));
                trap.push_back(shape);
        }
        else
        {
            if("ONEWALLED" != eQ::parameters["trapType"])
            {//2,3
                //LEFT VERTICAL WALL:
                shape = cpSpaceAddShape(
                            myParams.space, cpSegmentShapeNew(
                                staticBody, cpv(0.0, 0.0), cpv(0.0, h), TRAP_BORDER_WIDTH));
                    trap.push_back(shape);
                //RIGHT VERTICAL WALL:
                shape = cpSpaceAddShape(
                            myParams.space, cpSegmentShapeNew(
                                staticBody, cpv(w, 0.0), cpv(w, h), TRAP_BORDER_WIDTH));
                    trap.push_back(shape);
            }
            //close the top of the trap:
            if("TWOWALLED" != eQ::parameters["trapType"])
            {//1,3
                shape = cpSpaceAddShape(
                            myParams.space, cpSegmentShapeNew(
                                staticBody, cpv(0.0, h), cpv(w, h), TRAP_BORDER_WIDTH));
                    trap.push_back(shape);
            }
           //  //CODE TO ADD "BUMPS" TO THE TOP OF THE TRAP:
           //  // #define TRAPBUMPS
           //  #define NUMBUMPS 30
           // int tw = 2*width;
           // for(int i=0;i<tw/NUMBUMPS;i++)
           // {
           //     shape = cpSpaceAddShape(
           //                 space, cpCircleShapeNew(
           //                     staticBody, TRAP_BORDER_WIDTH/2.0, cpv(-w + NUMBUMPS*i, -h+TRAP_BORDER_WIDTH)));
           //         trap[trapSegments++] = shape;
           // }
         }
    }
    else{}//NO-WALLS CASE

    for(auto &shape : trap)
    {
        cpShapeSetElasticity(shape, TRAP_ELASTICITY);
        cpShapeSetFriction(shape, TRAP_FRICTION);
    }
}

bool cpmTrap::outsideTrap(std::shared_ptr<cpmEColi> cell)
{
    if("NOTRAP" == eQ::parameters["trapType"])
        return false;

    double y = double(cell->center.y);
    double x = double(cell->center.x);

    //should compute this in eQ class:  todo
    return (
            (x > w) ||
                (x < 0.0) ||
            (y > h) ||
                (y < 0.0)
            );
}
cpmTrap::~cpmTrap()
{
    for(auto &shape : trap)
        cpShapeFree(shape);
}
bool cpmTrap::updateModel(std::shared_ptr<cpmEColi> cell)
{
    return outsideTrap(cell);
}

//bool cpmTrap::updateModel(jChipmunk_EColi *cell)
//{
//    if(noTrap) return false;

//    // if(outsideTrap2(cell))
//    if(outsideTrap(cell))
//    {
//        if(removeOutsideTrap)
//        {
//            ;
//        }
//        //else:  TODO:
//        cpVect flow;
//        double y = double(cell->center.y);

//        if (y>highChannel)
//        {
//            flow = flowForceHigh;
//        }
//        else if (y<lowChannel)
//        {
//            flow = flowForceLow;
//        }
//    //        return true;
//        cell->applyForce(flow);
//        return true;
//    }
//    else
//    {
//        if(2==whichTrap)
//        {
//            if (cell->center.y < 0)
//            {
//                compAccumulator += cell->compression;
//            }
//        }

//    //space-filling parabolas: y = 2a(1 - x^2) - 1
//    //
//            //normalize to [-1,1] X [-1,1]
//            // float x = cell->center.x/width;
//            double y = cell->center.y/height;
//            // float yr = (2.0*frand()-1.0);

//            // float alpha = (-y + 1.0)/(2.0*(1.0 - x*x));

//            //y = 2a(1-x^2/a^2) - 1     y' = -4x/a
//            // float alpha2 = 0.25*((-y + 1) + sqrt((-y + 1)*(-y+1) + 16.0*x*x));


//            //only add force if point inside or on a=1
//    //        if(alpha2 <= 1.0)
//    //        if( (alpha2 <= 0.75) && (alpha2 > 0.05))
//        if( (frand() > 0.5) && (-y < 0.5) )
//        {
//    #define FORCE_FACTOR 50.0
//    //        cpVect force = cpvnormalize(cpv(xr, yr*4.0*alpha*x)) * cell->getScaledExpansionForce(10.5);
//    //            cpVect force = cpvnormalize(cpv(1.0, 4.0*x/alpha2)) * cell->getScaledExpansionForce(FORCE_FACTOR);
//            // cpVect force = cpv(xr, y+1.0+frand()) * cell->getScaledExpansionForce(50.0);

//    //FORCE OPTION (LATEST)
//        // float xr = (2.0*frand()-1.0);
//            // cpVect force = cpv(xr, y+1.0+frand()) * cell->getScaledExpansionForce(50.0);
//            // cell->applyForce(force);



//    //            cell->applyForce(force*frand());
//            //turn ON/OFF flow:
//        }
//        return false;
//    }
//}

////JAMMING TRAP CONSTRUCTOR, WHICH=1,2,... (OR NO TRAP, WHICH=0)
//cpmTrap::cpmTrap(cpmHabitat *habitat, int which)
//{
//    noTrap=(0==which) ? true:false;
//    space = habitat->space;
//    habitat->trap = this;

//    trapSegments = 0;
//    cpShape *shape;


//    whichTrap = which;
//    //no trap condition
//    if(0==which) return;
//    //else, add jamming trap skeleton:

//    compAccumulator = 0.0;
//    //JAMMING TRAP:
//    staticBody = cpSpaceGetStaticBody(space);
//    #define JAMHEIGHT 125
//    #define JAMWIDTH  150
//    #define EDGEOUTLETX 35//50
//    #define EDGE1X  55//75//150
//    #define EDGE1Y  (JAMHEIGHT+50) //350
//    #define EDGE2Y  (EDGE1Y+JAMHEIGHT+50)
//    //THREE WALLS TO TRAP:
//        shape = cpSpaceAddShape(
//                    space, cpSegmentShapeNew(
//                        staticBody, cpv(-JAMWIDTH/2,0), cpv(-JAMWIDTH/2,-JAMHEIGHT), TRAP_BORDER_WIDTH));
//            trap[trapSegments++] = shape;
//        shape = cpSpaceAddShape(
//                    space, cpSegmentShapeNew(
//                        staticBody, cpv(-JAMWIDTH/2,-JAMHEIGHT), cpv(JAMWIDTH/2,-JAMHEIGHT), TRAP_BORDER_WIDTH));
//            trap[trapSegments++] = shape;
//        shape = cpSpaceAddShape(
//                    space, cpSegmentShapeNew(
//                        staticBody, cpv(JAMWIDTH/2,-JAMHEIGHT), cpv(JAMWIDTH/2,0), TRAP_BORDER_WIDTH));
//            trap[trapSegments++] = shape;
//    //lower edges to outlet:
//        shape = cpSpaceAddShape(
//                    space, cpSegmentShapeNew(
//                        staticBody, cpv(-JAMWIDTH/2,0), cpv(-EDGEOUTLETX,0), TRAP_BORDER_WIDTH));
//            trap[trapSegments++] = shape;
//        shape = cpSpaceAddShape(
//                    space, cpSegmentShapeNew(
//                        staticBody, cpv(JAMWIDTH/2,0), cpv(EDGEOUTLETX,0), TRAP_BORDER_WIDTH));
//            trap[trapSegments++] = shape;
//    //initial path from outlet:  0->250
//        shape = cpSpaceAddShape(
//                    space, cpSegmentShapeNew(
//                        staticBody, cpv(-EDGEOUTLETX,0), cpv(-EDGEOUTLETX,JAMHEIGHT), TRAP_BORDER_WIDTH));
//            trap[trapSegments++] = shape;
//        shape = cpSpaceAddShape(
//                    space, cpSegmentShapeNew(
//                        staticBody, cpv(EDGEOUTLETX,0), cpv(EDGEOUTLETX,JAMHEIGHT), TRAP_BORDER_WIDTH));
//            trap[trapSegments++] = shape;
//    //final outlet: 10 um  500->600
//        shape = cpSpaceAddShape(
//                    space, cpSegmentShapeNew(
//                        staticBody, cpv(-EDGEOUTLETX,EDGE2Y), cpv(-EDGEOUTLETX,EDGE2Y+100), TRAP_BORDER_WIDTH));
//            trap[trapSegments++] = shape;
//        shape = cpSpaceAddShape(
//                    space, cpSegmentShapeNew(
//                        staticBody, cpv(EDGEOUTLETX,EDGE2Y), cpv(EDGEOUTLETX,EDGE2Y+100), TRAP_BORDER_WIDTH));
//            trap[trapSegments++] = shape;

//    if(1==which)
//    {
//    //edges from outlet:  no valve:  250->500
//        shape = cpSpaceAddShape(
//                    space, cpSegmentShapeNew(
//                        staticBody, cpv(-EDGEOUTLETX,250), cpv(-EDGEOUTLETX,500), TRAP_BORDER_WIDTH));
//            trap[trapSegments++] = shape;
//        shape = cpSpaceAddShape(
//                    space, cpSegmentShapeNew(
//                        staticBody, cpv(EDGEOUTLETX,250), cpv(EDGEOUTLETX,500), TRAP_BORDER_WIDTH));
//            trap[trapSegments++] = shape;
//    }
//    else if(2==which)
//    {
//            //left side
//                shape = cpSpaceAddShape(
//                            space, cpSegmentShapeNew(
//                                staticBody, cpv(-EDGEOUTLETX,JAMHEIGHT), cpv(-EDGE1X,EDGE1Y), TRAP_BORDER_WIDTH));
//                    trap[trapSegments++] = shape;
//                shape = cpSpaceAddShape(
//                            space, cpSegmentShapeNew(
//                                staticBody, cpv(-EDGE1X,EDGE1Y), cpv(-EDGE1X,EDGE2Y), TRAP_BORDER_WIDTH));
//                    trap[trapSegments++] = shape;
//                shape = cpSpaceAddShape(
//                            space, cpSegmentShapeNew(
//                                staticBody, cpv(-EDGE1X,EDGE2Y), cpv(-EDGEOUTLETX,EDGE2Y), TRAP_BORDER_WIDTH));
//                    trap[trapSegments++] = shape;
//            //right side:
//                shape = cpSpaceAddShape(
//                            space, cpSegmentShapeNew(
//                                staticBody, cpv(EDGEOUTLETX,JAMHEIGHT), cpv(EDGE1X,EDGE1Y), TRAP_BORDER_WIDTH));
//                    trap[trapSegments++] = shape;
//                shape = cpSpaceAddShape(
//                            space, cpSegmentShapeNew(
//                                staticBody, cpv(EDGE1X,EDGE1Y), cpv(EDGE1X,EDGE2Y), TRAP_BORDER_WIDTH));
//                    trap[trapSegments++] = shape;
//                shape = cpSpaceAddShape(
//                            space, cpSegmentShapeNew(
//                                staticBody, cpv(EDGE1X,EDGE2Y), cpv(EDGEOUTLETX,EDGE2Y), TRAP_BORDER_WIDTH));
//                    trap[trapSegments++] = shape;
//    }
//    for(int i=0; i<trapSegments; i++)
//    {
//        cpShapeSetElasticity(trap[i], TRAP_ELASTICITY);
//        cpShapeSetFriction(trap[i], TRAP_FRICTION);
//    }
//}


