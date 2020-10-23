#include "Strain.h"


//static data member for each strain type (can be initialized via enum within scope of class):
std::vector<bool> aspectRatioInvasionStrain::inductionFlags = {aspectRatioInvasionStrain::NUM_INDUCTIONFLAGS, false};
std::vector<bool> sendRecvStrain::inductionFlags            = {sendRecvStrain::NUM_INDUCTIONFLAGS, false};
std::vector<bool> MODULUSmodule::inductionFlags             = {MODULUSmodule::NUM_INDUCTIONFLAGS, false};
std::vector<bool> synchronousOscillator::inductionFlags     = {synchronousOscillator::NUM_INDUCTIONFLAGS, false};
synchronousOscillator::Data synchronousOscillator::environmentData;

std::vector<double>
synchronousOscillator::computeProteins
    (const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
{
    static bool setInitialValues = false;

    if( inductionFlags[SET_INITIAL_SYNTHASE_CONC] )
    {//set all cells:
        double xpos = params.baseData->x;
        double leftBoundary = environmentData.trapWidthMicrons/2 - environmentData.centerSliceWidth;
        double rightBoundary = environmentData.trapWidthMicrons/2 + environmentData.centerSliceWidth;
//        std::cout<<"checking cell at: "<<xpos<<"  "<<leftBoundary<<", "<<rightBoundary<<std::endl;

        if( (xpos > leftBoundary) && (xpos < rightBoundary) )
        {
            std::queue<double> q;
            for(size_t j=0; j<queueDepth; j++) {q.push(1.0e3);}
            HSL_tau.assign(HSL_tau.size(), q);
            //set the synthase conc. to 1.0 (max)
            iPROTEIN[RHLI] = eQ::Cell::nanoMolarToMoleculeNumber(1.0, lengthMicrons);
            std::cout<<"set synthase to 1 at: "<<xpos<<std::endl;
        }
        setInitialValues = true;
    }
    else if(!setInitialValues)
    {
        dHSL.assign(dHSL.size(), 0);
        return dHSL;
    }

    if(eHSL.empty()) return {};//error with no HSL, return empty vector

    //===========================================================================
    //===========================================================================

    computeConcentrations(eHSL, membraneRate, lengthMicrons);

    double alpha = double(eQ::data::parameters["hslProductionRate_C4"]);
    double delta = double(eQ::data::parameters["hslLeakProduction"]);
    double HK = 400.;
    double hn = 4.;
//    double AK = 10;
    double AKd = 400.;
    double gamma_dil = (log(2.0)/20.0);
    double gamma_deg = gamma_dil * double(eQ::data::parameters["gammaDegradationScale"]);

    ratio_H_tau = pow(tHSL[C4HSL]/HK, hn);
//    ratio_A     = iPROTEIN[AIIA]/AK;
    ratio_AKd   = iPROTEIN[AIIA]/AKd;

    deltaHSL[C4HSL]  = params.dt * (
                        alpha * iPROTEIN[RHLI]
                    -  20.0 * ratio_AKd/(1.0 + ratio_AKd) * iHSL[C4HSL])
                    - dHSL[C4HSL]                                  //MEMBRANE DIFFUSION
    ;
    deltaPROTEIN[RHLI]  = params.dt * (
                        delta
                    + gamma_dil * ratio_H_tau/(1.0 + ratio_H_tau)
    );
    deltaPROTEIN[AIIA]  = params.dt * (
                        0.25*(alpha * ratio_H_tau/(1.0 + ratio_H_tau))
                    - gamma_deg * iPROTEIN[AIIA]                  //set above
    );
    deltaPROTEIN[FP] = params.dt * (
                        gamma_dil * ratio_H_tau/(1.0 + ratio_H_tau)
    );

    conc[S]     = iPROTEIN[RHLI];//copy for data recording
    conc[A]     = iPROTEIN[AIIA];//copy for data recording
    conc[GFP]   = iPROTEIN[FP];//copy for data recording

//            ratio_I_tau = pow(I_tau/params->K.H, params->nH);//change to match that of C4
//            //C14 synthase:
//            delta[S] = dt * (
//                      (params->eta.S0 + params->eta.S1 * ratio_I_tau)/(1.0 + ratio_I_tau)
//                    - degradationProtein * conc[S] );
//            delta[I] = dt * (//C14HSL
//                    params->phi * conc[S]									//PRODUCTION
//                    - degradationHSL * conc[I] )							//DEGRADATION+DILUTION
//                    - dC14HSL;                                               //MEMBRANE DIFFUSION
            //hard code the receiver Prhl K50 for now:
            //C14 production from sender strain "priming signal" (want high sensitivity and output)

//            delta[FP] = params.dt * (
//                                gamma_d * ratio_H_tau/(1.0 + ratio_H_tau)
//                        );

    pushConcentrations();
    return dHSL;
}

std::vector<double>
sendRecvStrain::computeProteins
    (const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
{
    //===========================================================================
    //SENDER_RECEIVER
    //===========================================================================
    if( inductionFlags[INDUCTION] )
    {//set all cells:
        params.baseData->meanDivisionLength
                = double( eQ::data::parameters["mutantAspectRatioScale"])
                                           * double( eQ::data::parameters["defaultAspectRatioFactor"])
                                           * eQ::Cell::DEFAULT_DIVISION_LENGTH_MICRONS;
    }

    double trapWidthMicrons = double(eQ::data::parameters["simulationTrapWidthMicrons"]);
    double iptg = params.baseData->x / trapWidthMicrons;

    if(eHSL.empty()) return {};//error with no HSL, return empty vector

    //===========================================================================
    //===========================================================================

    computeConcentrations(eHSL, membraneRate, lengthMicrons);

    if(eQ::Cell::strainType::ACTIVATOR == params.whichType)
    {
//        double offset = 0.0;
//        double offset = 0.25;
//        double offset = 0.5;

        double scale = 1.0;
//        scale = iptg + offset;

        double alpha = scale * double(eQ::data::parameters["hslProductionRate_C4"]);

        //Sender cells produce C4HSL at maximum rate:
//        deltaHSL[0]  = dt * (double(eQ::data::parameters["hslProductionRate_C4"]) * 1.0)//PRODUCTION maximum=1.0
        deltaHSL[0]  = params.dt * (alpha)//PRODUCTION maximum=1.0
                - dHSL[0];                                               //MEMBRANE DIFFUSION
    }
    else if(eQ::Cell::strainType::REPRESSOR == params.whichType)
    {
        if("DUAL_SENDER_RECEIVER" == eQ::data::parameters["simType"])
        {
//            ratio_I_tau = pow(I_tau/params->K.H, params->nH);//change to match that of C4
//            //C14 synthase:
//            delta[S] = dt * (
//                      (params->eta.S0 + params->eta.S1 * ratio_I_tau)/(1.0 + ratio_I_tau)
//                    - degradationProtein * conc[S] );
//            delta[I] = dt * (//C14HSL
//                    params->phi * conc[S]									//PRODUCTION
//                    - degradationHSL * conc[I] )							//DEGRADATION+DILUTION
//                    - dC14HSL;                                               //MEMBRANE DIFFUSION
        }
        else
        {
            //hard code the receiver Prhl K50 for now:
            //C14 production from sender strain "priming signal" (want high sensitivity and output)
            double HK = 900.0;
            double hn = 15.0;
            ratio_H_tau = pow(tHSL[0]/HK, hn);

            double gamma_d = (log(2)/20.0);

            deltaHSL[0]  =      - dHSL[0];                                               //MEMBRANE DIFFUSION
            delta[FP] = params.dt * (
                                gamma_d * ratio_H_tau/(1.0 + ratio_H_tau)
                        );
        }
    }

    pushConcentrations();
    return dHSL;
}

std::vector<double>
aspectRatioInvasionStrain::computeProteins
    (const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
{
    double dummy = membraneRate.size() * eHSL.size() * lengthMicrons;
    dummy *= dummy;//removes "unused variable" wornings
    //===========================================================================
    //STATIC_ASPECTRATIO
    //===========================================================================

    if("STATIC_ASPECTRATIO" == eQ::data::parameters["simType"])
    {
        if( inductionFlags[ASPECTRATIO_INDUCTION] )
        {
            params.baseData->meanDivisionLength
                    = double( eQ::data::parameters["mutantAspectRatioScale"])
                                               * double( eQ::data::parameters["defaultAspectRatioFactor"])
                                               * eQ::Cell::DEFAULT_DIVISION_LENGTH_MICRONS;
    //        params.baseData->divisionLength
    //                = params.baseData->meanDivisionLength;//set to divde at this length immediately
        }

    }
    //===========================================================================
    //INDUCED_DYNAMIC_ASPECTRATIO
    //===========================================================================
    else if("INDUCED_DYNAMIC_ASPECTRATIO" == eQ::data::parameters["simType"])
    {
        if( inductionFlags[ASPECTRATIO_INDUCTION] && (eQ::Cell::strainType::ACTIVATOR == getStrainType()) )
        {
            params.baseData->meanDivisionLength
                    = double( eQ::data::parameters["mutantAspectRatioScale"])
                                               * double( eQ::data::parameters["defaultAspectRatioFactor"])
                                               * eQ::Cell::DEFAULT_DIVISION_LENGTH_MICRONS;
    //        params.baseData->divisionLength
    //                = params.baseData->meanDivisionLength;//set to divde at this length immediately
        }
    }

    return std::vector<double>(eHSL.size(), 0.0);//no gene circuit to update, return empty vector
}
std::vector<double>
aspectRatioOscillator::computeProteins
    (const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
{
    double aspectRatioThresh    = double( eQ::data::parameters["aspectRatioThresholdHSL"]);
    double aspectRatioScaling   = double( eQ::data::parameters["defaultAspectRatioFactor"]);
    double mutantScaling        = double(eQ::data::parameters["mutantAspectRatioScale"]);

    static bool flagCheck=false;

//    if(inductionFlags[ASPECTRATIO_INDUCTION])
//    {
        double hslValue = (eQ::Cell::strainType::REPRESSOR == getStrainType())
                ? getDelayedHSL(Strain::hsl::C4) : getDelayedHSL(Strain::hsl::C14);

        if(hslValue > aspectRatioThresh)
        {
            aspectRatioScaling *= mutantScaling;
            if(!flagCheck)
            {
                flagCheck = true;
                std::cout<<"\t hit aspect ratio flag...changing to: "<<aspectRatioScaling<<std::endl;
            }
        }


    params.baseData->meanDivisionLength = aspectRatioScaling * eQ::Cell::DEFAULT_DIVISION_LENGTH_MICRONS;


    //===========================================================================
    //===========================================================================

    computeConcentrations(eHSL, membraneRate, lengthMicrons);


    double gamma_d = (log(2)/20.0);

//    double pHinThresh           = 200.0;
//    double pHinThresh           = double(eQ::data::parameters["aspectRatioThresholdHSL"]);
    double pHinThresh           = aspectRatioThresh;

//    double pHControlThresh      = 1200.0;
//    double pLacControlThresh    = 1000.0;
//    double placThresh           = 0.25;
    double placThresh           = 0.5;

//    double ratio_L = pow(L_tau/placThresh, 10.0);

//    double responseScale = 30.0;
    //for using just 2 HSLs:
    const size_t    Ain = hslType::C14HSL;
    const size_t    Rin = hslType::C4HSL;


    if(eQ::Cell::strainType::ACTIVATOR == params.whichType)
    {
//        //HSL production from other strain "priming signal" (want high sensitivity and output)
//        double ratio_Hcontrol   = pow(tHSL[2]/pHControlThresh, 10.0);
//        //intra-strain signal
//        deltaHSL[3] = dt * responseScale * (double(eQ::parameters["hslProductionRate_C14"]) * ratio_Hcontrol/(1.0 + ratio_Hcontrol))
//                - dHSL[3];
//        //LacI
        double ratio_Hin        = pow(iHSL[Ain]/pHinThresh, 10.0);
        delta[L] = params.dt * (gamma_d * ratio_Hin/(1.0 + ratio_Hin));
//        double ratio_Hlac = pow(tHSL[3]/pLacControlThresh, 10.0);  //want just sensitive enough to switch on to ensure switching off
//        delta[L] = dt * (gamma_d * ratio_Hlac/(1.0 + ratio_Hlac));

        //HSL signal to other strain (inhibited by other strain via toggle switch topology)
        double ratio_L = pow(conc[L]/placThresh, 10.0);
        deltaHSL[Rin] = params.dt * (double(eQ::data::parameters["hslProductionRate_C4"]))
                *  1.0/(1.0 + ratio_L)
                - dHSL[Rin];

        //report HSL from receiver cells:
        deltaHSL[Ain] = - dHSL[Ain];                                         //MEMBRANE DIFFUSION
    }
    else if(eQ::Cell::strainType::REPRESSOR == params.whichType)
    {
//        double ratio_Hcontrol = pow(tHSL[0]/pHControlThresh, 10.0);
//        //intra-strain signal
//        deltaHSL[1] = dt * responseScale * (double(eQ::parameters["hslProductionRate_C14"]) * ratio_Hcontrol/(1.0 + ratio_Hcontrol))
//                - dHSL[1];
        //LacI
        double ratio_Hin        = pow(iHSL[Rin]/pHinThresh, 10.0);
        delta[L] = params.dt * (gamma_d * ratio_Hin/(1.0 + ratio_Hin));
//        double ratio_Hlac = pow(tHSL[1]/pLacControlThresh, 10.0);  //want just sensitive enough to switch on to ensure switching off
//        delta[L] = dt * (gamma_d * ratio_Hlac/(1.0 + ratio_Hlac));


        //priming signal output to other strain:
        double ratio_L = pow(conc[L]/placThresh, 10.0);
        deltaHSL[Ain] = params.dt * (double(eQ::data::parameters["hslProductionRate_C14"]))
                *  1.0/(1.0 + ratio_L)
                - dHSL[Ain];

        //report HSL from receiver cells:
        deltaHSL[Rin] = - dHSL[Rin];                                         //MEMBRANE DIFFUSION
    }

    pushConcentrations();
    return dHSL;
}



double modulusUpdate(double x);
std::vector<double>
MODULUSmodule::computeProteins
    (const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
{

    double iptg = modulusUpdate(params.baseData->x);

    computeConcentrations(eHSL, membraneRate, lengthMicrons);

    const double hn = 2.0;//hill exponent
        double gamma_d = (log(2)/20.0);


        if("MODULUS_1" == eQ::data::parameters["simType"])
        {
            //convert signal "c" to fixed activating signal for Plac
            double Klac = 30.0; //nanoMolar Plac K50 via Lutz et al.

            //FIXED K50:
//            double ratioPlac = pow(iptg/Klac, hn);
            //use lognormal Plac K50:
//            double ratioPlac = pow(iptg/lnKL, hn);
            //use lognormal INPUT: lnKL is now the input IPTG value (set on cell birth/division):
            double ratioPlac = pow(lnKL/Klac, hn);

            //FIXED K50 for output promoter:
            ratio_H_tau = pow(H_tau/double(eQ::data::parameters["lnmean"]), hn);
//            //FIXED K50 = 0.5 for feedback promoter:
//            double ratio_H_tau2 = pow(H_tau/0.5, hn);

            //use generated KH values from lognormal distributions (set on cell birth):
//            ratio_H_tau = pow(H_tau/lnKH, hn);
            //for feedback circuit:
//            double ratio_H_tau2 = pow(H_tau/lnKH2, hn);

            double fbalpha = double(eQ::data::parameters["MODULUS_FEEDBACK_FRACTION"]);//convex combination version
            double fbscale = double(eQ::data::parameters["MODULUS_FEEDBACK_STRENGTH"]);//multiple for feedback production

            //WITH DIFFUSION:
            if( (eQ::data::parameters["MODULUS_option"] == "+D+F") || (eQ::data::parameters["MODULUS_option"] == "+D-F") )
            {
                if(eQ::data::parameters["MODULUS_option"] == "+D-F")
                {//NO FEEDBACK:
                    //C4 synthase: [S] \in [0, 1]
                    delta[S] = params.dt * ( //rhlI under PL-lac promoter (direct activation by "c"=iptg, ratio set explicitly above)
                             (gamma_d * ratioPlac)/(1.0  + ratioPlac));
                }
                else
                {//WITH FEEDBACK:
                    //rhlI under PL-lac promoter (direct activation by "c"=iptg, ratio set explicitly above)
                    delta[S] = params.dt * (
                           gamma_d * ratioPlac/(1.0  + ratioPlac)
                        +  fbscale * gamma_d * ratio_H_tau/(1.0 + ratio_H_tau));
//                    (1.0 - fbalpha)*gamma_d * ratioPlac/(1.0  + ratioPlac)
//                 +  fbalpha * gamma_d * ratio_H_tau/(1.0 + ratio_H_tau));
                }

                //READOUT:
                delta[FP] = params.dt * ( //Prhl promoter
                           gamma_d * (ratio_H_tau)/(1.0 + ratio_H_tau));

//                dC4HSL = params.dt * (//positive direction is outwards from cell
//                        params->DH * (conc[H] - eHSL[0]));    //MEMBRANE DIFFUSION

                delta[H] = params.dt * (double(eQ::data::parameters["hslProductionRate_C4"]) * conc[S])//PRODUCTION
                        - dHSL[0];                                               //MEMBRANE DIFFUSION
//                        - dC4HSL;                                               //MEMBRANE DIFFUSION
            }
            //NO DIFFUSION: use AraC
            else if( (eQ::data::parameters["MODULUS_option"] == "-D+F") || (eQ::data::parameters["MODULUS_option"] == "-D-F") )
            {
                //for recording only:
                delta[L] = params.dt * ( gamma_d * iptg );

                //AraC production:
                if(eQ::data::parameters["MODULUS_option"] == "-D-F")
                {//NO FEEDBACK:
                    delta[S] = params.dt * ( //S=AraC under PL-lac promoter (H is identically S here)
                             (gamma_d * ratioPlac)/(1.0  + ratioPlac));
                }
                //with positive feedback, add AraC under Pbad promoter
                else if(eQ::data::parameters["MODULUS_option"] == "-D+F")
                {//WITH FEEDBACK:
                    delta[S] = params.dt * (
                            gamma_d * ratioPlac/(1.0  + ratioPlac)
                         +  fbscale * gamma_d * ratio_H_tau/(1.0 + ratio_H_tau));
//                           (1.0 - fbalpha)*gamma_d * ratioPlac/(1.0  + ratioPlac)
//                        +  fbalpha * gamma_d * ratio_H_tau/(1.0 + ratio_H_tau));
                }
                //READOUT
                delta[FP] = params.dt * ( //Pbad promoter for AraC
                    //normalize by the expected division time (to cancel dilution due to growth)
                    gamma_d * (ratio_H_tau)/(1.0 + ratio_H_tau));

                //MEMBRANE DIFFUSION OF HSL (degradation via spatial diffusion is done on grid)
                dHSL[0] = 0.0;    //MEMBRANE DIFFUSION
//                dC4HSL = 0.0;    //MEMBRANE DIFFUSION
                //H=AraC here:
                double gammaT = double(eQ::data::parameters["gammaT_C4"]);//need to scale back out the HSL decay rate:
                delta[H] = params.dt * (gamma_d * double(eQ::data::parameters["hslProductionRate_C4"])/gammaT * conc[S]);//PRODUCTION
            }
            else ;//bad option
        }
        else if("MODULUS_2" == eQ::data::parameters["simType"])
        {
            //convert signal "c" to fixed activating signal for Plac
            double Klac = 30.0; //nanoMolar Plac K50 via Lutz et al.
            Klac = 60.0;
            double ratioPlac = pow(iptg/Klac, hn);
            //FIXED K50 for output promoter:
//            double HK = double(eQ::data::parameters["lnmean"]) * 8.0;
//            double HK = 1.8e4;
            double HK = 3.0e4;
            ratio_H_tau = pow(tHSL[0]/HK, hn);

            //C4 synthase:
            delta[S] = params.dt * ( //rhlI under PL-lac promoter (direct activation by "c"=iptg, ratio set explicitly above)
                     (gamma_d * ratioPlac)/(1.0  + ratioPlac));

            deltaHSL[0] = params.dt * (double(eQ::data::parameters["hslProductionRate_C4"]) * conc[S])//PRODUCTION
                       - dHSL[0];                                               //MEMBRANE DIFFUSION

            deltaPROTEIN[tetR] = params.dt * ((gamma_d * ratio_H_tau)/(1.0  + ratio_H_tau));
            double ratio_tetR = pow(tPROTEIN[tetR]/0.5, hn);

            //READOUT:
            delta[FP] = params.dt * gamma_d * (
                    ratioPlac/(1.0  + ratioPlac)//LOCAL EXCITATION
//                    +    1.0/(1.0 + ratio_tetR)     //GLOBAL INHIBITION
                    +    1.0/(1.0 + ratio_H_tau)     //GLOBAL INHIBITION
            );

            conc[L] = ratioPlac;//use as temp storage for testing value
            conc[H] = ratio_H_tau;
            conc[LEGI_A] = ratio_tetR;

        }

        pushConcentrations();
        return dHSL;
}

double modulusUpdate(double x)
{
    //    a = v(4)/D;
    //    cxR = c0 / (1 - exp(a*L))  ...
    //        * (exp(a*x) - exp(a*L)) + offset;
    double cl = double(eQ::data::parameters["simulationChannelLengthLeft"]);//already scaled
    double cr = double(eQ::data::parameters["simulationChannelLengthRight"]);
    double tw = double(eQ::data::parameters["simulationTrapWidthMicrons"]);
    double xt = cl+cr+tw;

    double vx = double(eQ::data::parameters["simulationFlowRate"]);
    //um^2/min for IPTG (relative to C4HSL a la Pai and You, using molecular weight vs. C4HSL)
    double D = 3.0e4 * sqrt(159.0/238.0) * double(eQ::data::parameters["diffusionScaling"]);
    double a = vx/D;
    double expaL = exp(a*xt);//exp(ax) evaluated at x=L (right side of device)
    double expax = exp(a*(x+cl));//exp(ax) evaluated at cell position x (add left-side offset)

    double iptgOffset   = 0.0;//low side of gradient at media channel right (Modulus c0 offset of input gradient)
    double iptgLeft     = +100.0;//high side of gradient at media channel left (+ difference from right side)

    double thisC = (iptgLeft/(1.0 - expaL)) * (expax - expaL) + iptgOffset;

        if ("MODULUS_2" == eQ::data::parameters["simType"])
        {
//                        eQ::data::parameters["MODULUS_IPTG"] = thisC;
//                        const double iptgOffset = 0.0;
//                        const double iptgOffset = 10.0;
//                        const double iptgOffset = 20.0;
//                        const double iptgOffset = 30.0;
//                        const double iptgOffset = 40.0;
//                        const double iptgOffset = 50.0;
            const double iptgOffset = 60.0;
            double iptgMin = 0.0;
            double iptgMax = 30.0;
//                        double iptgMax = 00.0;//control input
//                        const double iptgMin = 20.0;
//                        const double iptgMax = 80.0;
//                        const double iptgMin = 10.0;
//                        const double iptgMax = 100.0;
//                        const double iptgMin = 100.0;
//                        const double iptgMax = 10.0;
//                        const double iptgMin = 30.0;
//                        const double iptgMax = 30.0;
            iptgMin += iptgOffset;
            iptgMax += iptgOffset;
            double iptg = (x/tw)*(iptgMax-iptgMin) + iptgMin;
            return iptg;
        }
        return  double(eQ::data::parameters["MODULUS_IPTG"]);
}


/*
void Strain::solveLinearSystem(double &H_i, double &H_e, double y1, double y2, double cellFraction)
{
    //MEMBRANE DIFFUSION OF HSL (degradation via spatial diffusion is done on grid)
    //                double newHi, newHe;
    //                solveLinearSystem(newHi, newHe, dt*double(eQ::data::parameters["hslProductionRate_C4"]) * conc[S] + conc[H], C4ext, rho_cell);
    //                delta[H] = newHi - conc[H];
    //                dC4HSL = newHe  - C4ext;
    double a,b,c,d;

    double ri = 0.5/cellFraction;
    double re = 0.5/(1.0 - cellFraction);
    double di = params->DH * ri;
    double de = params->DH * re;
    a = 1.0 + dt*(di);
    b = -dt*di;
    c = -dt*de;
    d = 1.0 + dt*de;

    H_e = (y1 - y2*a/c)/(b - a*d/c);
    H_i = (y1 - b*H_e)/a;
}
*/
////////////////////////////////////////////////////////////////////////////////
//				loadParams():
////////////////////////////////////////////////////////////////////////////////
void loadParams(eQ::Cell::strainType which, struct dso_parameters &params, struct rates &rates)
{//SELECT PROMOTER STRENGTHS HERE:
    if(eQ::Cell::strainType::ACTIVATOR == which)
	{//ACTIVATOR:
		params.eta.S0 	= rates.eta[R0][STRONG];
		params.eta.S1 	= rates.eta[R1][STRONG];
		params.eta.FP0 	= rates.eta[F0][STRONG];
		params.eta.FP1 	= rates.eta[F1][STRONG];

		params.K.H 		= rates.K[STRONG].H;
		params.K.I 		= rates.K[WEAK].I;
		params.K.L  	= rates.K[WEAK].L;

		params.phi 		= rates.phiH;
	}
    else if(eQ::Cell::strainType::REPRESSOR == which)
	{//REPRESSOR:
		params.eta.S0 	= rates.eta[C0][WEAK];
		params.eta.S1 	= rates.eta[C1][WEAK];
		params.eta.FP0 	= rates.eta[Y0][WEAK];
		params.eta.FP1 	= rates.eta[Y1][WEAK];

		params.K.H 		= rates.K[STRONG].H;
		params.K.I 		= rates.K[WEAK].I;
		params.K.L  	= rates.K[WEAK].L;

		params.phi 		= rates.phiI;
	}
	//THESE HAVE ONLY ONE VALUE: POPULATE params directly:
	params.eta.L0 	= rates.eta[L0][ONLY];
	params.eta.L1 	= rates.eta[L1][ONLY];
	params.eta.A0 	= rates.eta[A0][ONLY];
	params.eta.A1 	= rates.eta[A1][ONLY];

	params.K.A 			= rates.K[ONLY].A;
	params.K.C 			= rates.K[ONLY].C;

	params.nH 			= rates.nH;
	params.nL 			= rates.nL;
	params.nI 			= rates.nI;
	params.dC 			= rates.dC;
	params.dA 			= rates.dA;
	params.dilution 	= rates.dilution;
	params.mu_e 		= rates.mu_e;
	params.DH 			= rates.DH;
	params.DI 			= rates.DI;
	params.m 			= rates.m;
}
////////////////////////////////////////////////////////////////////////////////
//				initRates():
////////////////////////////////////////////////////////////////////////////////
void initRates(double *basalRate, struct rates &rates)
{
// MASTER TABLE OF RATES FROM CHEN (2015)

	// //basal rate chosen from histogram in Chen:
	// basalRate[SR] = 10.0;
	// basalRate[SC] = 10.0;
	// basalRate[SL] = 10.0;
	// basalRate[SA] = 10.0;
	// basalRate[SF] = 10.0;
	// basalRate[SY] = 2.0;
	// basalRate[ClpXP] = 1400.0;

	//values posted in Jae Kim's update paper:
	basalRate[SR] = 3.06;
	basalRate[SC] = 37.23;
	basalRate[SL] = 4.52;
//    basalRate[SA] = 258.0;
    basalRate[SA] = 9.54;
    basalRate[SF] = 113.0;
	basalRate[SY] = 6.8;
	basalRate[ClpXP] = 1820.0;

	rates.eta[R0][WEAK] 		= 1.0 * basalRate[SR];
	rates.eta[R0][MEDIUM] 		= 3.04 * basalRate[SR];
	rates.eta[R0][STRONG] 		= 20.13 * basalRate[SR];
	rates.eta[R0][LAC]			= 177.44 * basalRate[SR];

	rates.eta[R1][WEAK] 		= 624.44 * basalRate[SR];
	rates.eta[R1][MEDIUM] 		= 574.97 * basalRate[SR];
	rates.eta[R1][STRONG] 		= 367.48 * basalRate[SR];
	rates.eta[R1][LAC]			= 0.0 * basalRate[SR];

	rates.eta[C0][WEAK] 		= 1.0 * basalRate[SC];
	rates.eta[C0][MEDIUM] 		= 3.04 * basalRate[SC];
	rates.eta[C0][STRONG] 		= 20.13 * basalRate[SC];

	rates.eta[C1][WEAK] 		= 624.44 * basalRate[SC];
	rates.eta[C1][MEDIUM] 		= 574.97 * basalRate[SC];
	rates.eta[C1][STRONG] 		= 367.48 * basalRate[SC];

	rates.eta[F0][WEAK] 		= 1.0 * basalRate[SF];
	rates.eta[F0][MEDIUM] 		= 3.04 * basalRate[SF];
	rates.eta[F0][STRONG] 		= 20.13 * basalRate[SF];
	rates.eta[F0][LAC]			= 177.44 * basalRate[SF];

	rates.eta[F1][WEAK] 		= 624.44 * basalRate[SF];
	rates.eta[F1][MEDIUM] 		= 574.97 * basalRate[SF];
	rates.eta[F1][STRONG] 		= 367.48 * basalRate[SF];
	rates.eta[F1][LAC]			= 0.0 * basalRate[SF];

	rates.eta[Y0][WEAK] 		= 1.0 * basalRate[SY];
	rates.eta[Y0][MEDIUM] 		= 41.8 * basalRate[SY];
	rates.eta[Y1][WEAK] 		= 1713.0 * basalRate[SY];
	rates.eta[Y1][MEDIUM] 		= 197.49 * basalRate[SY];

	//THESE HAVE ONLY ONE VALUE: POPULATE params directly:
	rates.eta[L0][ONLY] 		= 1.0 * basalRate[SL]; 
	rates.eta[L1][ONLY] 		= 1735.47 * basalRate[SL]; 
	rates.eta[A0][ONLY] 		= 27.03 * basalRate[SA]; 
	rates.eta[A1][ONLY] 		= 141.61 * basalRate[SA];

	rates.K[WEAK].H				= 16599.0;
	rates.K[MEDIUM].H			= 10333.0;
	rates.K[STRONG].H			= 5937.0;

	rates.K[WEAK].L				= 47.7;
	rates.K[LAC].L				= 85.38;

	rates.K[WEAK].I				= 2357.3;
	rates.K[MEDIUM].I			= 594.23;

	//THESE HAVE ONLY ONE VALUE: POPULATE params directly:
    rates.K[ONLY].A				= 5110.0 * 1000.0;
    rates.K[ONLY].C 			= 1300.0;

//    rates.nH = 4;
//	rates.nL = 2;
//	rates.nI = 4;
    rates.nH = 4.;
    rates.nL = 2.;
    rates.nI = 4.;

    rates.dC = 1.8 * basalRate[ClpXP];
	rates.dA = 2257.0;
	rates.dilution = log(2.0)/25.0;
	rates.mu_e = 0.1;
	// rates.mu_e = 0.01;

	rates.DH = 3.0;
	rates.DI = 2.1;
	rates.phiH = 16.0;
	rates.phiI = 2.0;
	rates.m = log(2.0)/3.0;
}
