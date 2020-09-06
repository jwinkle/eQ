#include "Strain.h"


//static data member for each strain type (can be initialized via enum within scope of class):
std::vector<bool> aspectRatioInvasionStrain::inductionFlags = {aspectRatioInvasionStrain::NUM_INDUCTIONFLAGS, false};
std::vector<bool> sendRecvStrain::inductionFlags            = {sendRecvStrain::NUM_INDUCTIONFLAGS, false};
std::vector<bool> MODULUSmodule::inductionFlags             = {MODULUSmodule::NUM_INDUCTIONFLAGS, false};
std::vector<bool> aspectRatioOscillator::inductionFlags     = {aspectRatioOscillator::NUM_INDUCTIONFLAGS, false};
std::vector<bool> parB_MotherStrain::inductionFlags         = {parB_MotherStrain::NUM_INDUCTIONFLAGS, false};


std::vector<double>
sendRecvStrain::computeProteins
    (const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
{


    //===========================================================================
    //SENDER_RECEIVER
    //===========================================================================
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
                std::cout<<std::endl;
                std::cout<<"\t hit aspect ratio flag...changing to: "
                        <<aspectRatioScaling
                       <<" with strain: "<<getStrainType()
                       <<std::endl;
                std::cout<<std::endl;
            }
        }
//    }


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

double parB_MotherStrain::growthRateScaling()
{
    const double growthArrestThreshold = double(eQ::data::parameters["growthArrestThreshold"]);
//    const double tetRThreshold = double(eQ::data::parameters["tetRThreshold"]);
    const double lacIThreshold = double(eQ::data::parameters["lacIThreshold"]);

    return (tHSL[C14] > growthArrestThreshold) && (tPROTEIN[_lacI] < lacIThreshold)
//    return (tHSL[C14] > growthArrestThreshold)
            ? 0 : 1;
}
std::vector<double>
parB_MotherStrain::computeProteins
    (const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
{
    computeConcentrations(eHSL, membraneRate, lengthMicrons);

    const double hTet = 4;
    const double KTet = 0.25;
    double ratioTetR = pow(tPROTEIN[_tetR]/KTet, hTet);

    const double hCin = 4;
    const double KCin = 200;
    double ratioCin = pow(tHSL[C14]/KCin, hCin);

    const double parBThreshold = double(eQ::data::parameters["parBThreshold"]);

    if(eQ::Cell::strainType::ACTIVATOR == params.whichType)
    {
        deltaPROTEIN[_tetR] = params.dt * cellGrowthRate;
        deltaPROTEIN[_lacI] = params.dt * cellGrowthRate;
        deltaHSL[C4]        = params.dt * double(eQ::data::parameters["hslProductionRate_C4"]);
        if(inductionFlags[INDUCTION])
            parB_losePlasmid = (getDelayedHSL(Strain::hsl::C4) > parBThreshold);
    }

    if(inductionFlags[INDUCTION])
        deltaHSL[C14]  = params.dt * double(eQ::data::parameters["hslProductionRate_C14"])
                *  1/(1 + ratioTetR)//hill function repressor
                * ( 1
                    +  1 * ratioCin/(1 + ratioCin) //feedback
                );

    deltaPROTEIN[_tetR]  -=  params.dt * (iPROTEIN[_tetR] * cellGrowthRate);
    deltaPROTEIN[_lacI]  -=  params.dt * (iPROTEIN[_lacI] * cellGrowthRate);

    //LACTONASE DECAY OF HSL:
//    deltaHSL[C4]  -=  params.dt * (iHSL[C4] * aiiA_decayRate);
//    deltaHSL[C14]  -=  params.dt * (iHSL[C14] * aiiA_decayRate);

    //MEMBRANE DIFFUSION
    deltaHSL[C4]  -=  dHSL[C4];
    deltaHSL[C14]  -=  dHSL[C14];

    conc[tetR] = iPROTEIN[_tetR];//copy to old data structure for now (for recording)
    conc[L] = iPROTEIN[_lacI];//copy to old data structure for now (for recording)

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

