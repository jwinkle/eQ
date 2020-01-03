
//needed to set body shape polygon points and update them each step:
#include "../../Chipmunk-7.0.1/include/chipmunk/chipmunk_private.h"
#include "../../Chipmunk-7.0.1/include/chipmunk/chipmunk_unsafe.h"

//this must be included *after* the chipmunk *private/unsafe.h files, above
#include "cpmEColi.h"



//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//updated to microns:
#define RATCHET_QUANTUM    0.05
#define COMPRESSION_GAP    0.1
//old: pixels
//#define RATCHET_QUANTUM    0.5
//#define COMPRESSION_GAP    1.0
////#define RATCHET_QUANTUM    0.1
////#define COMPRESSION_GAP    0.1

//TODO: set these to private members and get/set
extern double       g_lambdaFitness;
extern double       g_compressionLimit;
extern unsigned int g_indexValue;
extern unsigned int g_overComp;


//==============================================================================
//==============================================================================
void jcpBodyUpdateVelocity(cpBody *body, cpVect gravity, cpFloat damping, cpFloat dt)
{
    // Skip kinematic bodies.
    if(cpBodyGetType(body) == CP_BODY_TYPE_KINEMATIC) return;

    cpAssertSoft(body->m > 0.0f && body->i > 0.0f, 
        "Body's mass and moment must be positive to simulate. (Mass: %f Moment: %f)", 
        body->m, body->i);

    //get the user data for this cell
    cpmEColi::bodyData_t *pBodyData
            =  (cpmEColi::bodyData_t *)cpBodyGetUserData(body);

    // cpVect Fbym = body->f * body->m_inv;
    cpVect Fbym     = body->f;
//    body->v         = Fbym * (1.0/GAMMA_FLUID);
//    body->v         = Fbym * pBodyData->gammaFluidInverse;
    body->v         = Fbym * (1.0/pBodyData->gammaFluid);
    body->w         = body->w*damping + body->t*body->i_inv*dt;

    // Reset forces.
    body->f = cpvzero;
    body->t = cpFloat(0.0);
//    cpAssertSaneBody(body);
}
//==============================================================================
//==============================================================================
void cpmEColi::setStaticCell()
{//sets the cell body halves to "static" -- they should not be updated, so set flag
    cpBodySetType(bodyA, CP_BODY_TYPE_STATIC);
    cpBodySetType(bodyB, CP_BODY_TYPE_STATIC);
    cellIsStatic = true;
}
void cpmEColi::resetDynamicCell()
{//sets the cell body halves back to "dynamic"
    cpBodySetType(bodyA, CP_BODY_TYPE_DYNAMIC);
    cpBodySetType(bodyB, CP_BODY_TYPE_DYNAMIC);
    cellIsStatic = false;
}
//==============================================================================
void cpmEColi::computeCellBoudaryEdges()
{
    edges[0] = vertsA[1] - vertsA[0];
    edges[1] = vertsA[2] - vertsA[1];
    edges[2] = vertsA[3] - vertsA[2];
    edges[3] = vertsA[0] - vertsA[3];
//    edges[0] = vertsB[1] - vertsA[0] + cpv(separationDistance, 0.0);
//    edges[1] = vertsB[2] - vertsB[1];
//    edges[2] = vertsA[3] - vertsB[2] + cpv(separationDistance, 0.0);
//    edges[3] = vertsA[0] - vertsA[3];
}
//==============================================================================
cpmEColi::cpmEColi(const cpmEColi::params &ecoli)
            : space(ecoli.space)
{
    /*
         * Constructor for the chipmunk model.
         * Input data is through Init structure pointer
         * Takes 'length' parameter as actual cell length and forms the
         * model from this and a 'width' -- fixed at 1um.
         *
    */

//    parameterData.gammaFluidInverse = 1.0/ecoli.gammaFluidParameter;
    parameterData.gammaFluid = ecoli.gammaFluidParameter;
    cellIsStatic = false;

    k_spring            = ecoli.kspring;

    initialMass         = ecoli.mass;
    initialMoment       = ecoli.moment;
    shapeLength         = ecoli.length - ecoli.width;
    shapeHeight         = ecoli.width;
    length              = ecoli.length;
    initialLength       = ecoli.length;
    center              = cpv(ecoli.x, ecoli.y);
    angle               = ecoli.angle;
    velocity            = cpv(ecoli.vx, ecoli.vy);
    angularVelocity     = ecoli.av;

    posA                = center;
    posB                = center;
    separationDistance  = 0.0;
    separationVelocity  = cpvzero;

    compression         = 0.0;

    //rest-length increment per time step for spring model
    dRL = ecoli.dRL;

    radius = shapeHeight * 0.5;//radius of a cell pole
    offset = shapeLength * 0.5;//offset to center of pole from cell center
    //    nextRatchet = 0.0;
    nextRatchet = RATCHET_QUANTUM;

    //GROOVE JOINT:
    cpVect upperLeft = cpv(-offset, radius);
    cpVect upperRight = cpv(offset, radius);
    cpVect lowerRight = cpv(offset, -radius);
    cpVect lowerLeft = cpv(-offset, -radius);

    // //DAMPED SPRING:
    // cpFloat springx = offset;
    // cpFloat springy = radius;

    //    cpVect springAnchor1A = cpv(-springx, springy);
    //    cpVect springAnchor1B = cpv(springx, springy);
    //    cpVect springAnchor2A = cpv(-springx, -springy);
    //    cpVect springAnchor2B = cpv(springx, -springy);

    //chipmunk spring RL is set relative to box edges; set so RL=shapeLength here:
    cpSpringRestLength = shapeLength;
   

    //==============================================================================
    //      CELL BODY HALVES:  mass,moment,center,angle,vel.,angular_vel
    //==============================================================================
    bodyA = cpSpaceAddBody(
                space, cpBodyNew(
                    initialMass, initialMoment));
    bodyB = cpSpaceAddBody(
                space, cpBodyNew(
                    initialMass, initialMoment));
    cpBodySetPosition(bodyA, center);
        cpBodySetAngle(bodyA, angle);
            cpBodySetVelocity(bodyA, velocity);
    cpBodySetPosition(bodyB, center);
        cpBodySetAngle(bodyB, angle);
            cpBodySetVelocity(bodyB, velocity);

    cpBodySetVelocityUpdateFunc(bodyA,cpBodyVelocityFunc(jcpBodyUpdateVelocity));
    cpBodySetVelocityUpdateFunc(bodyB,cpBodyVelocityFunc(jcpBodyUpdateVelocity));

    cpBodySetUserData(bodyA,  (cpDataPointer *) &parameterData);
    cpBodySetUserData(bodyB,  (cpDataPointer *) &parameterData);

    //lab-frame velocity
    velA = velocity;
    velB = velocity;
    //==============================================================================
    //      CELL SHAPES:  box(length,width), pole(radius,offset)
    //==============================================================================
    //NOTE: SHAPE "A" IS THE LHS SHAPE, WHOSE PRIMARY POLE POINTS LEFT;  SHAPE B IS LIKEWISE RHS
    //EACH SHAPE HAS A "BACK FILLED" TRAILING RECTANGLE THAT IS UPDATED VIA THE RATCHET ALGORITHM
    //THERE IS ALSO A BACKSIDE POLE THAT IS ATTACHED THAT PREVENTS COMPRESSION OF THE TWO CELL HALVES A,B
    //vertices of initial box, for updates to growth:
    vertsA[0]=upperLeft;
    vertsA[1]=upperRight;
    vertsA[2]=lowerRight;
    vertsA[3]=lowerLeft;
        vertsB[0]=upperLeft;
        vertsB[1]=upperRight;
        vertsB[2]=lowerRight;
        vertsB[3]=lowerLeft;

    computeCellBoudaryEdges();

    //BOX CENTER FOR CELL:
        //from: https://chipmunk-physics.net/release/ChipmunkLatest-Docs/
//        Because boxes are so common in physics games, Chipmunk provides shortcuts to create box shaped polygons.
//                The boxes will always be centered at the center of gravity of the body you are attaching them to.
//                Adding a small radius will bevel the corners and can significantly reduce problems where the box gets stuck on seams in your geometry.
//                If you want to create an off-center box, you will need to use cpPolyShapeNew() or cpPolyShapeInit().
    //create and add the box and circle shapes on the body:
    shapeCount = 0;
    boxA = cpSpaceAddShape(
                space, cpBoxShapeNew(
                  bodyA, shapeLength, shapeHeight, cpFloat(0.0)));
    shapes[shapeCount++] = boxA;
    boxB = cpSpaceAddShape(
                space, cpBoxShapeNew(
                  bodyB, shapeLength, shapeHeight, cpFloat(0.0)));
    shapes[shapeCount++] = boxB;
    //CIRC ENDS FOR CELL:
        //from: https://chipmunk-physics.net/release/ChipmunkLatest-Docs/
//            body is the body to attach the circle to, offset is the offset from the body’s center of gravity in body local coordinates.
    //"primary" poles for the cell halves (expansion-side):
    capA = cpSpaceAddShape(
                space, cpCircleShapeNew(
                    bodyA, radius, cpv(-offset, cpFloat(0.0))));
    shapes[shapeCount++] = capA;
    capB = cpSpaceAddShape(
                space, cpCircleShapeNew(
                    bodyB, radius, cpv(offset, cpFloat(0.0))));
    shapes[shapeCount++] = capB;
    //"secondary" poles (back-filled side):
    capA2 = cpSpaceAddShape(
                space, cpCircleShapeNew(
                    bodyA, radius, cpv(offset, cpFloat(0.0))));
    shapes[shapeCount++] = capA2;
    capB2 = cpSpaceAddShape(
                space, cpCircleShapeNew(
                    bodyB, radius, cpv(-offset, cpFloat(0.0))));
    shapes[shapeCount++] = capB2;
    for(unsigned int i=0;i<shapeCount;i++)
    {
        cpShapeSetElasticity(shapes[i], SHAPE_ELASTICITY);
        cpShapeSetFriction(shapes[i], SHAPE_FRICTION);
    }
    //============================================================
    //      CELL GROOVE JOINTS:
    //============================================================
    constraintCount = 0;
    //GROOVE JOINTS
    //lock the two bodies together using symmetric groove joints (and set to not collide the two bodies):
    grooveJointA = cpSpaceAddConstraint(
                space, cpGrooveJointNew(
                    bodyA,bodyB, upperLeft, upperRight, upperLeft));
    constraints[constraintCount++] = grooveJointA;
    grooveJointB = cpSpaceAddConstraint(
                space, cpGrooveJointNew(
                    bodyB,bodyA, lowerLeft, lowerRight, lowerRight));
    constraints[constraintCount++] = grooveJointB;

    grooveJointA2 = cpSpaceAddConstraint(
                space, cpGrooveJointNew(
                    bodyA,bodyB, lowerLeft, lowerRight, lowerLeft));
    constraints[constraintCount++] = grooveJointA2;
    grooveJointB2 = cpSpaceAddConstraint(
                space, cpGrooveJointNew(
                    bodyB,bodyA, upperLeft, upperRight, upperRight));
    constraints[constraintCount++] = grooveJointB2;


        // //============================================================
        // //      CELL SPRING EXPANSION:
        // //============================================================
        // //SPRING JOINTS:
        // // #ifndef NOSPRING
        // #define cpSPRING_DAMPING    1.0f
        // #define cpSPRING_STIFFNESS  (100.0)
        // // #endif
        // springJoint1 = cpSpaceAddConstraint(
        //             space, cpDampedSpringNew(
        //                 bodyA, bodyB, springAnchor1A, springAnchor1B,
        //                 cpSpringRestLength, cpFloat(cpSPRING_STIFFNESS), cpFloat(cpSPRING_DAMPING)));
        // constraints[constraintCount++] = springJoint1;
        // // cpConstraintSetMaxBias(springJoint1, 0.1);

        // springJoint2 = cpSpaceAddConstraint(
        //             space, cpDampedSpringNew(
        //                 bodyA, bodyB, springAnchor2A, springAnchor2B,
        //                 cpSpringRestLength, cpFloat(cpSPRING_STIFFNESS), cpFloat(cpSPRING_DAMPING)));
        // constraints[constraintCount++] = springJoint2;
        // // cpConstraintSetMaxBias(springJoint2, 0.1);



    for(unsigned int i=0;i<constraintCount;i++)
    {
        cpConstraintSetCollideBodies(constraints[i], cpFalse);
        cpConstraintSetMaxForce(constraints[i], GROOVE_JOINT_MAXFORCE);
        cpConstraintSetErrorBias(constraints[i], GROOVE_JOINT_ERRORBIAS);
        cpConstraintSetMaxBias(constraints[i], GROOVE_JOINT_MAXBIAS);
    }

    setSpringRestLength(length);
}//end jChipmunk_EColi::jChipmunk_EColi()
//==============================================================================
//==============================================================================

//==============================================================================
//==============================================================================
cpmEColi::~cpmEColi(void)
{
    for(unsigned int i=0;i<constraintCount;i++)
    {
        cpSpaceRemoveConstraint(space, constraints[i]);
        cpConstraintFree(constraints[i]);
    }

    for(unsigned int i=0;i<shapeCount;i++)
    {
        cpSpaceRemoveShape(space, shapes[i]);
        cpShapeFree(shapes[i]);
    }
    cpSpaceRemoveBody(space, bodyA);
    cpSpaceRemoveBody(space, bodyB);
    cpBodyFree(bodyA);
    cpBodyFree(bodyB);
}

bool cpmEColi::pointIsInCell(std::pair<double, double> point)
{
    cpVect lpoint = cpBodyWorldToLocal(bodyA, cpv(point.first, point.second));

//    cpVect p0 = lpoint - vertsB[1] + cpv(separationDistance, 0.0);
//    cpVect p1 = lpoint - vertsB[2] + cpv(separationDistance, 0.0);
    cpVect p0 = lpoint - vertsA[1];
    cpVect p1 = lpoint - vertsA[2];
    cpVect p2 = lpoint - vertsA[3];
    cpVect p3 = lpoint - vertsA[0];

    return ((cpvcross(edges[0], p0) < 0.0)
            && (cpvcross(edges[1], p1) < 0.0)
            && (cpvcross(edges[2], p2) < 0.0)
            && (cpvcross(edges[3], p3) < 0.0));

    //Find the distance from point to shape. If the point is inside of the shape, the distance will be negative and equal to the depth of the point.
//    struct cpPointQueryInfo info;
//    cpShapeCacheBB(boxA);
//    cpFloat responseA = cpShapePointQuery(boxA, p, &info);
//    if (responseA < 0.0) return true;
//    cpShapeCacheBB(boxB);
//    cpFloat responseB = cpShapePointQuery(boxB, p, &info);
//    if (responseB < 0.0) return true;
//    return false;
}
void cpmEColi::setVelocity(cpVect vel)
{
    cpBodySetVelocity(bodyA, vel); velA = vel;
    cpBodySetVelocity(bodyB, vel); velB = vel;
}
void cpmEColi::setBodyVelocities(cpVect vA, cpVect vB)
{
    cpBodySetVelocity(bodyA, vA); velA =vA;
    cpBodySetVelocity(bodyB, vB); velB=vB;
}

void cpmEColi::applyForce(cpVect force)
{//adds to the current force on the object
    cpBodyApplyForceAtWorldPoint(bodyA,  force, cpBodyLocalToWorld(bodyA, cpvzero));
    cpBodyApplyForceAtWorldPoint(bodyB,  force, cpBodyLocalToWorld(bodyB, cpvzero));
}

void cpmEColi::setSpringRestLength(cpFloat RL)
{
    //viruaal spring RL is set rel. to the full length of the cell
    springRestLength = RL;

    // //chipmunk (damped) spring implementation:
    // cpSpringRestLength = RL - shapeHeight;
    // cpDampedSpringSetRestLength(springJoint1, cpSpringRestLength);
    // cpDampedSpringSetRestLength(springJoint2, cpSpringRestLength);
}
void cpmEColi::incSpringRestLength(cpFloat deltaRL)
{
    setSpringRestLength(springRestLength + deltaRL);
}
cpFloat cpmEColi::getSpringRestLength(void)
{
    // return cpDampedSpringGetRestLength(springJoint1);
    return springRestLength;
}
//==============================================================================
//==============================================================================
int cpmEColi::updateModel(void)
{
    //read the post-step data from the chipm. model:
    //position, angle, velocity, angularVelocity

    posA = cpBodyGetPosition(bodyA);
    posB = cpBodyGetPosition(bodyB);
        center = (posA + posB) * 0.5;                   //center of the cell wrt. the body positions
        separationDistance = cpvdist(posA, posB);       //distance between centers of mass
        length = separationDistance + initialLength;    //length of the whole cell, incl. poles
    angleA = (cpBodyGetAngle(bodyA));
    angleB = (cpBodyGetAngle(bodyB));
        angle = (angleA + angleB) * 0.5;        //compass heading wrt. body angles (average)
    velA = cpBodyGetVelocity(bodyA);
    velB = cpBodyGetVelocity(bodyB);
        velocity = (velA + velB) * 0.5;                 //vel. of center of mass of cell
    angularVelocity = (cpBodyGetAngularVelocity(bodyA)
                       + cpBodyGetAngularVelocity(bodyB)) * 0.5;

    separationAngle         = cpvtoangle(posB - posA);          //compass heading wrt. centers of mass
    separationVelocity      = (velB - velA); //vel. wrt. cell frame (expansion velocity); uses bodyB as + dir.
    separationVelocityAngle = cpvtoangle(velB - velA);//compass heading wrt sepraration vel.
    separationSpeed         = cpvlength(separationVelocity);
    projectionSpeed         = cpvdot(separationVelocity, cpvnormalize(velB));//signed speed of separation
    cellExpanding           = (projectionSpeed >= 0.0) ? true:false;
    //    if(false == cellExpanding) separationSpeed = -separationSpeed;

    //============================================================
    //  RATCHET ALGORITHM:
    //============================================================
    if(separationDistance > (nextRatchet + COMPRESSION_GAP))
    //        if(separationDistance > (nextRatchet + gaps[icell]))
    {
        cpFloat newOffset = offset + nextRatchet;
        nextRatchet = nextRatchet + RATCHET_QUANTUM;
    //update the rectangle for each cell half (only back-filled half):
    //        vertsA[0]=upperLeft;
            vertsA[1]=cpv(newOffset, radius);
            vertsA[2]=cpv(newOffset, -radius);
    //        vertsA[3]=lowerLeft;
            vertsB[0]=cpv(-newOffset, radius);;
    //        vertsB[1]=upperRight;
    //        vertsB[2]=lowerRight;
            vertsB[3]=cpv(-newOffset, -radius);;

        computeCellBoudaryEdges();

    // #include "chipmunk_unsafe.h"  <--included in header
    //note:  see  https://chipmunk-physics.net/release/ChipmunkLatest-Docs/#cpShape-Modifying
        //quoting on "Modifying Shapes":
        //"The short answer is that you can’t because the changes would be only picked up as a change
            //to the position of the shape’s surface, but not its velocity. The long answer is that you
            //can using the “unsafe” API as long as you realize that doing so will result in unrealistic
            //physical behavior. These extra functions are defined in a separate header chipmunk_unsafe.h."
        //redraw the shapes. NB: these should not generate new collisions (on non-expanding side)
        cpPolyShapeSetVertsRaw(boxA, 4, vertsA);
        cpPolyShapeSetVertsRaw(boxB, 4, vertsB);
        cpCircleShapeSetOffset(capA2, cpv(newOffset, 0.0));
        cpCircleShapeSetOffset(capB2, cpv(-newOffset, 0.0));
        //expand the grooves and anchors to match new lengths for shapes:
        cpGrooveJointSetGrooveB(grooveJointA, vertsA[1]);
            cpGrooveJointSetAnchorB(grooveJointA, vertsB[0]);
        cpGrooveJointSetGrooveB(grooveJointA2, vertsA[2]);
            cpGrooveJointSetAnchorB(grooveJointA2, vertsB[3]);
        cpGrooveJointSetGrooveA(grooveJointB, vertsB[3]);
            cpGrooveJointSetAnchorB(grooveJointB, vertsA[2]);
        cpGrooveJointSetGrooveA(grooveJointB2, vertsB[0]);
            cpGrooveJointSetAnchorB(grooveJointB2, vertsA[1]);
    }

    compression = springRestLength - length;

    return 0;
}
//============================================================
//  SPRING EXPANSION FORCE ALGORITHM:
//============================================================
int cpmEColi::updateExpansionForce(double restLengthIncrement)
{//cell update will have been called before this (updates compression)
//at dt=0.001 min, dRL increment is 0.0001 um => 20k x dRL = 2um (doubling) in 20 minutes

    if(cellIsStatic) return 0;//don't update cell if it has been frozen "static"

//to scale by cell length, we want 2um/20 minutes = 20pixels/20min = 1pixel/min = dt pixels/dt mins, say at 3um length
    cpFloat climit = g_compressionLimit;

    //Exponential growth: l(t) = l0 * 2^(t/Td) ==> dl = dt * l * 2^(dt/20) * ln(2)/20
    const double expScale = log(2.0)/20.0 * pow(2.0,dRL/20.0);//recall dRL is dt, so a little of a hack here for the source of dt
    //need to fix this to read length locally:
//    double dRL_exp = dRL * parent->jget_length() * expScale;

        if(true)//bypass growth rate regulation
        // if (compression < climit)
        {
             incSpringRestLength(restLengthIncrement);
//            incSpringRestLength(dRL_exp);
            // springForce = cpv(KSPRING * (compression+dRL), 0.0f);
            // springRestLength += dRL;
        }
        else
        {
            g_overComp++;
            cpFloat c = (compression-climit)/climit;
            if (c > 1.0) c=1.0;
             cpFloat scaled_dRL = (1.0-c)*dRL;
//            cpFloat scaled_dRL = (1.0-c)*dRL_exp;
            incSpringRestLength(scaled_dRL);
            // springForce = cpv(KSPRING * (compression + scaled_dRL), 0.0f);
            // springRestLength += scaled_dRL;
        }

    //VIRTUAL SPRING IMPLEMENTATION:
//        cpVect springForce = cpv(KSPRING * (springRestLength - length), cpFloat(0.0));
        cpVect springForce = cpv(k_spring * (springRestLength - length), cpFloat(0.0));
    cpBodyApplyForceAtLocalPoint(bodyA,  cpvneg(springForce), cpvzero);
    cpBodyApplyForceAtLocalPoint(bodyB,  springForce, cpvzero);

    return 0;
}
