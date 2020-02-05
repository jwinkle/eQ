#include "Strain.h"

Strain::~Strain() {} //required to be defined outside of declaration to avoid obscure "undefined reference to vtable" linking errors

Strain::Strain(eQ::strainType strainType, const struct dso_parameters *p, const double timestep, const double nodesPerMicronScale)
    : whichType(strainType), params(p), dt(timestep), nodesPerMicron(nodesPerMicronScale), inductionFlag(false)
{
	//NOTE:  valid for fixed dt/delay times only
    double delay = double(eQ::parameters["promoterDelayTimeMinutes"]);
    queueDepth = size_t(ceil(delay/dt));
//    queueDepth = size_t(ceil(TAU_DELAY/dt));

	for(size_t i=0; i<queueDepth; i++)
	{
		qH_tau.push(0.0);
		qL_tau.push(0.0);
		qI_tau.push(0.0);
	}
	for(int i=0; i<numConcentrations; i++)
	{
		conc[i]=0.0;//don't omit this!
		delta[i]=0.0;
	}
}
void Strain::divideProteins(std::shared_ptr<Strain> daughter, double fracToDaughter)
{
    //NOTE: the size of all data structures is determined by constructor of derived class
	double protein;
	double rd = fracToDaughter;
	double rm = (1.0 - fracToDaughter);

    //OLD DATA STRUCTURES:
    for(size_t i=0; i<queueDepth; i++)
	{
        //NOTE:  DON'T DIVIDE SINCE QUEUE STORES CONCENTRATION:
        protein = qH_tau.front();
            qH_tau.pop();
            qH_tau.push(protein);
                daughter->qH_tau.pop();
                daughter->qH_tau.push(protein);
        protein = qL_tau.front();
            qL_tau.pop();
            qL_tau.push(protein);
                daughter->qL_tau.pop();
                daughter->qL_tau.push(protein);
        protein = qI_tau.front();
            qI_tau.pop();
            qI_tau.push(protein);
                daughter->qI_tau.pop();
                daughter->qI_tau.push(protein);
	}

	for(size_t i=0;i<numConcentrations;i++)
    {//DIVIDE PROTEINS BY FRACTION
        daughter->conc[i] =  rd*conc[i];
        conc[i] *= rm;
	}

    //NEW DATA STRUCTURES:
    for(size_t i(0); i< iHSL.size(); i++)
    {
        //DIVIDE HSL# BY FRACTION
        daughter->iHSL[i] = rd*iHSL[i];
        iHSL[i] *= rm;
        //NOTE:  DON'T DIVIDE SINCE QUEUE STORES CONCENTRATION:
        daughter->tHSL[i] = tHSL[i];
        for(size_t j=0; j<queueDepth; j++)
        {//pop-push j times for each i:
            protein = HSL_tau[i].front();
                HSL_tau[i].pop();
                HSL_tau[i].push(protein);
                    daughter->HSL_tau[i].pop();
                    daughter->HSL_tau[i].push(protein);
        }
    }

    for(size_t i(0); i< iPROTEIN.size(); i++)
    {
        //DIVIDE PROTEIN# BY FRACTION
        daughter->iPROTEIN[i] = rd*iPROTEIN[i];
        iPROTEIN[i] *= rm;
        //NOTE:  DON'T DIVIDE SINCE QUEUE STORES CONCENTRATION:
        daughter->tPROTEIN[i] = tPROTEIN[i];
        for(size_t j=0; j<queueDepth; j++)
        {//pop-push j times for each i:

            protein = PROTEIN_tau[i].front();
                PROTEIN_tau[i].pop();
                PROTEIN_tau[i].push(protein);
                    daughter->PROTEIN_tau[i].pop();
                    daughter->PROTEIN_tau[i].push(protein);
        }
    }
}
std::vector<double>
Strain::computeProteins
    (double C4ext, double C14ext, double lengthMicrons)
{

    //data:  conc[i] stores protein # for each protein between timesteps (covert to conc. first)
//    double explicitDilution = params->dilution;
    double explicitDilution         = 0.0;//set to zero when using ABM  with cell growth
    //we store protein count, convert to nanoMolar concentration: 0.602 protein/um^3 = 1nM    
    double oneProteinToNanoMolar = eQ::proteinNumberToNanoMolar(1.0, lengthMicrons);

    if(lengthMicrons == 0.0)
    {
        std::cout<<"ERROR: scale for computeProteins() set to 0!"<<std::endl;
        return std::vector<double>{};
    }
    for(int i=0; i<numConcentrations; i++)
    {
        conc[i] *= oneProteinToNanoMolar;//convert to proper concentration
    }
    //reactions on promoters use delayed concentrations (stored as t-tau concentrations, not #, thus no scaling):
    //promoter Hill functions "started" transcription/translation  tau time ago, so use those concentrations
    H_tau = qH_tau.front();  qH_tau.pop();
    L_tau = qL_tau.front();  qL_tau.pop();
    I_tau = qI_tau.front();  qI_tau.pop();
    //compute Hill function ratios once:
     ratio_H_tau = pow(H_tau/params->K.H, params->nH);
     ratio_L_tau = pow(L_tau/params->K.L, params->nL);
     ratio_I_tau = pow(I_tau/params->K.I, params->nI);

	totalProtein = 
            params->K.C + conc[S]+ conc[L] + conc[A] + conc[FP] + conc[M];
	totalHSL =
            params->K.A + conc[H] + conc[I];
	degradationProtein =
            (params->dC / totalProtein) + explicitDilution;//explicitly model dilution here
	degradationHSL =
            (params->dA / totalHSL) * conc[A] + explicitDilution;//explicitly model dilution here

    //C4 or C14 synthase:
    delta[S] = dt * ( //rhlI or cinI both have rhl/lac hybrid promoter
			(params->eta.S0 + params->eta.S1 * ratio_H_tau)/(1.0 + ratio_H_tau + ratio_L_tau)
			- degradationProtein * conc[S] );
	//LacI
    delta[L] = dt * ( //lacI has cin promoter
			(params->eta.L0 + params->eta.L1 * ratio_I_tau)/(1.0 + ratio_I_tau)
			- degradationProtein * conc[L] );
	//AiiA
    delta[A] = dt * ( //aiiA has cin promoter
			(params->eta.A0 + params->eta.A1 * ratio_I_tau)/(1.0 + ratio_I_tau)
			- degradationProtein * conc[A] );
	//Mature FP
	delta[M] = dt * (
			params->m * conc[FP] 
			- degradationProtein * conc[M] );
    //MEMBRANE DIFFUSION OF HSL (degradation via spatial diffusion is done on grid)
    dC4HSL = dt * (
            params->DH * (conc[H] - C4ext));    //MEMBRANE DIFFUSION
    dC14HSL = dt * (
            params->DI * (conc[I] - C14ext));   //MEMBRANE DIFFUSION

    if(eQ::strainType::ACTIVATOR == whichType)
	{
        //CFP --> tracks C4HSL synthase:
        delta[FP] = dt * ( //activator has rhl/lac hybrid promoter
				(params->eta.FP0 + params->eta.FP1 * ratio_H_tau)/(1.0 + ratio_H_tau + ratio_L_tau)
				- (degradationProtein + params->m) * conc[FP] );
        delta[H] = dt * (//C4HSL
				params->phi * conc[S]									//PRODUCTION
                - degradationHSL * conc[H] )							//DEGRADATION+DILUTION
                - dC4HSL;                                               //MEMBRANE DIFFUSION
        delta[I] = dt * (//C14HSL
                - degradationHSL * conc[I] ) 							//DEGRADATION+DILUTION
                - dC14HSL;                                               //MEMBRANE DIFFUSION
    }
    else if(eQ::strainType::REPRESSOR == whichType)
	{
        //YFP --> tracks C14HSL synthase:
        delta[FP] = dt * ( //repressor has cin/lac promoter
				(params->eta.FP0 + params->eta.FP1 * ratio_I_tau)/(1.0 + ratio_I_tau + ratio_L_tau)
				- (degradationProtein + params->m) * conc[FP] );
        delta[I] = dt * (//C14HSL
				params->phi * conc[S]									//PRODUCTION
                - degradationHSL * conc[I] )							//DEGRADATION+DILUTION
                - dC14HSL;                                               //MEMBRANE DIFFUSION
        delta[H] = dt * (//C4HSL
                - degradationHSL * conc[H] )							//DEGRADATION+DILUTION
                - dC4HSL;                                               //MEMBRANE DIFFUSION
    }

    for(int i=0;i<numConcentrations;i++)
    {
        conc[i] += delta[i];
        if (conc[i] < 0.0)  conc[i] = 0.0;
        else conc[i] /= oneProteinToNanoMolar;//store as protein #
    }
    //push new HSL and LacI values to queue: [STORE AS CONCENTRATION in nM]:
    qH_tau.push(conc[H] * oneProteinToNanoMolar);
    qL_tau.push(conc[L] * oneProteinToNanoMolar);
    qI_tau.push(conc[I] * oneProteinToNanoMolar);

    return std::vector<double> {dC4HSL/oneProteinToNanoMolar, dC14HSL/oneProteinToNanoMolar};//return HSL # changes
}

std::vector<double>
Strain::computeProteins
    (const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
{
    computeConcentrations(eHSL, membraneRate, lengthMicrons);

    //        //compute Hill function ratios once:
//    ratio_H_tau = pow(H_tau/params->K.H, params->nH);
//    ratio_L_tau = pow(L_tau/params->K.L, params->nL);
//    ratio_I_tau = pow(I_tau/params->K.I, params->nI);

}


void Strain::computeConcentrations(const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
{
    //COMMON CODE TO INITIALIZE DATA, COMPUTE CONVERSIONS FROM # TO CONC., POP DELAY QUEUES, COMPUTE MEMBRANE FLUX:

    if(lengthMicrons == 0.0)
    {
        std::cout<<"ERROR: scale for computeProteins() set to 0!"<<std::endl;
        return;
    }
    nanoMolarPerMolecule = eQ::proteinNumberToNanoMolar(1.0, lengthMicrons);
    //    double nanoMolarPerMolecule = eQ::proteinNumberToNanoMolar(1.0, computeAverageCellLength());

    //OLD PROTEIN ARRAYS:
    for(int i=0; i<numConcentrations; i++)
    {//concentrations are stored as protein # over time steps in the conc[] array
        conc[i] *= nanoMolarPerMolecule;//convert to proper concentration (units: nM)
        delta[i] = 0.0;
    }
    //OLD DELAY QUEUES:
    //        //reactions on promoters use delayed concentrations (stored as t-tau concentrations, not #, thus no scaling):
//    H_tau = qH_tau.front();  qH_tau.pop();
//    L_tau = qL_tau.front();  qL_tau.pop();
//    I_tau = qI_tau.front();  qI_tau.pop();

    //NEW DATA STRUCTURES:
    for(size_t i(0); i< iHSL.size(); i++)
    {
        deltaHSL[i] =0.0;//hsl generation deltas
        iHSL[i] *= nanoMolarPerMolecule;//convert to proper concentration (units: nM)
        //reactions on promoters use delayed concentrations (stored as t-tau concentrations, not #, thus no scaling):
        tHSL[i] = HSL_tau[i].front();  HSL_tau[i].pop();
        //MEMBRANE DIFFUSION OF HSL (degradation via spatial diffusion is done on grid)
        dHSL[i] = dt * (
                    membraneRate[i] * (iHSL[i] - eHSL[i]));    //MEMBRANE DIFFUSION out
    }
    for(size_t i(0); i< iPROTEIN.size(); i++)
    {
        deltaPROTEIN[i] = 0.0;// generation deltas
        iPROTEIN[i] *= nanoMolarPerMolecule;//convert to proper concentration (units: nM)
        //reactions on promoters use delayed concentrations (stored as t-tau concentrations, not #, thus no scaling):
        tPROTEIN[i] = PROTEIN_tau[i].front();  PROTEIN_tau[i].pop();
    }
}
void Strain::pushConcentrations()
{
    //OLD DATA STRUCTURES:
    for(int i=0;i<numConcentrations;i++)
    {
        conc[i] += delta[i];
        if (conc[i] < 0.0)  conc[i] = 0.0;
        else conc[i] /= nanoMolarPerMolecule;//store as protein #
    }
    //NEW DATA STRUCTURES:
    for(size_t i(0); i< iHSL.size(); i++)
    {
        iHSL[i] += deltaHSL[i];
        if (iHSL[i] < 0.0)
        {
            iHSL[i] = 0.0;
        }
        HSL_tau[i].push(iHSL[i]);//store as concentration
        iHSL[i] /= nanoMolarPerMolecule;//store as protein #
    }
    for(size_t i(0); i< iPROTEIN.size(); i++)
    {
        iPROTEIN[i] += deltaPROTEIN[i];
        if (iPROTEIN[i] < 0.0) iPROTEIN[i] = 0.0;
        PROTEIN_tau[i].push(iPROTEIN[i]);//store as concentration
        iPROTEIN[i] /= nanoMolarPerMolecule;//store as protein #
    }
}
std:: vector<double>
sendRecvStrain::computeProteins
    (const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
{
    computeConcentrations(eHSL, membraneRate, lengthMicrons);

    if(eQ::strainType::ACTIVATOR == whichType)
    {
        double offset = 0.0;
//        double offset = 0.25;
//        double offset = 0.5;
        double scale = double(eQ::parameters["MODULUS_IPTG"]) + offset;
        //over-ride
        scale = 1.0;
        double alpha = scale * double(eQ::parameters["hslProductionRate_C4"]);

        //Sender cells produce C4HSL at maximum rate:
//        deltaHSL[0]  = dt * (double(eQ::parameters["hslProductionRate_C4"]) * 1.0)//PRODUCTION maximum=1.0
        deltaHSL[0]  = dt * (alpha)//PRODUCTION maximum=1.0
                - dHSL[0];                                               //MEMBRANE DIFFUSION
    }
    else if(eQ::strainType::REPRESSOR == whichType)
    {
        if("DUAL_SENDER_RECEIVER" == eQ::parameters["simType"])
        {
            ratio_I_tau = pow(I_tau/params->K.H, params->nH);//change to match that of C4
            //C14 synthase:
            delta[S] = dt * (
                      (params->eta.S0 + params->eta.S1 * ratio_I_tau)/(1.0 + ratio_I_tau)
                    - degradationProtein * conc[S] );
            delta[I] = dt * (//C14HSL
                    params->phi * conc[S]									//PRODUCTION
                    - degradationHSL * conc[I] )							//DEGRADATION+DILUTION
                    - dC14HSL;                                               //MEMBRANE DIFFUSION
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
            delta[FP] = dt * (
                                gamma_d * ratio_H_tau/(1.0 + ratio_H_tau)
                        );
        }
    }

    pushConcentrations();
    return dHSL;
}

std:: vector<double>
aspectRatioInvasionStrain::computeProteins
    (const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
{
    computeConcentrations(eHSL, membraneRate, lengthMicrons);

    double gamma_d = (log(2)/20.0);

    double pHinThresh           = 1000.0;
    double pHControlThresh      = 1200.0;
    double pLacControlThresh    = 1000.0;
    double placThresh           = 0.25;

    double ratio_L = pow(L_tau/placThresh, 10.0);

    double responseScale = 30.0;

    if(eQ::strainType::ACTIVATOR == whichType)
    {
        //HSL production from other strain "priming signal" (want high sensitivity and output)
        double ratio_Hcontrol   = pow(tHSL[2]/pHControlThresh, 10.0);
        //intra-strain signal
        deltaHSL[3] = dt * responseScale * (double(eQ::parameters["hslProductionRate_C14"]) * ratio_Hcontrol/(1.0 + ratio_Hcontrol))
                - dHSL[3];
//        //LacI
        double ratio_Hin        = pow(tHSL[2]/pHinThresh, 10.0);
        delta[L] = dt * (gamma_d * ratio_Hin/(1.0 + ratio_Hin));
//        double ratio_Hlac = pow(tHSL[3]/pLacControlThresh, 10.0);  //want just sensitive enough to switch on to ensure switching off
//        delta[L] = dt * (gamma_d * ratio_Hlac/(1.0 + ratio_Hlac));

        //HSL signal to other strain (inhibited by other strain via toggle switch topology)
        deltaHSL[0] = dt * (double(eQ::parameters["hslProductionRate_C4"]))
                *  1.0/(1.0 + ratio_L)
                - dHSL[0];

        //report HSL from receiver cells:
        deltaHSL[1] = - dHSL[1];                                         //MEMBRANE DIFFUSION
        deltaHSL[2] = - dHSL[2];                                         //MEMBRANE DIFFUSION
    }
    else if(eQ::strainType::REPRESSOR == whichType)
    {
        double ratio_Hcontrol = pow(tHSL[0]/pHControlThresh, 10.0);
        //intra-strain signal
        deltaHSL[1] = dt * responseScale * (double(eQ::parameters["hslProductionRate_C14"]) * ratio_Hcontrol/(1.0 + ratio_Hcontrol))
                - dHSL[1];
        //LacI
        double ratio_Hin        = pow(tHSL[0]/pHinThresh, 10.0);
        delta[L] = dt * (gamma_d * ratio_Hin/(1.0 + ratio_Hin));
//        double ratio_Hlac = pow(tHSL[1]/pLacControlThresh, 10.0);  //want just sensitive enough to switch on to ensure switching off
//        delta[L] = dt * (gamma_d * ratio_Hlac/(1.0 + ratio_Hlac));


        //priming signal output to other strain:
        deltaHSL[2] = dt * (double(eQ::parameters["hslProductionRate_C4"]))
                *  1.0/(1.0 + ratio_L)
                - dHSL[2];

        //report HSL from receiver cells:
        deltaHSL[0] = - dHSL[0];                                         //MEMBRANE DIFFUSION
        deltaHSL[3] = - dHSL[3];                                         //MEMBRANE DIFFUSION
    }

    //for recording into existing data files:
    conc[H] = iHSL[0];
    conc[I] = iHSL[1];

    for(int i=0;i<numConcentrations;i++)
    {
        conc[i] += delta[i];
        if (conc[i] < 0.0)  conc[i] = 0.0;
        conc[i] /= nanoMolarPerMolecule;//store as protein #
    }
    for(size_t i(0); i< iHSL.size(); i++)
    {
        iHSL[i] += deltaHSL[i];
        if (iHSL[i] < 0.0)
        {
            iHSL[i] = 0.0;
        }
        HSL_tau[i].push(iHSL[i]);//store as concentration
        iHSL[i] /= nanoMolarPerMolecule;//store as protein #
    }
    qL_tau.push(conc[L] * nanoMolarPerMolecule);

    return dHSL;
}

std:: vector<double>
MODULUSmodule::computeProteins
    (const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
{
    computeConcentrations(eHSL, membraneRate, lengthMicrons);

    const double hn = 2.0;//hill exponent
        double gamma_d = (log(2)/20.0);


        if("MODULUS_1" == eQ::parameters["simType"])
        {
            //convert signal "c" to fixed activating signal for Plac
            iptg = double(eQ::parameters["MODULUS_IPTG"]);
            double Klac = 30.0; //nanoMolar Plac K50 via Lutz et al.

            //FIXED K50:
//            double ratioPlac = pow(iptg/Klac, hn);
            //use lognormal Plac K50:
//            double ratioPlac = pow(iptg/lnKL, hn);
            //use lognormal INPUT: lnKL is now the input IPTG value (set on cell birth/division):
            double ratioPlac = pow(lnKL/Klac, hn);

            //FIXED K50 for output promoter:
            ratio_H_tau = pow(H_tau/double(eQ::parameters["lnmean"]), hn);
//            //FIXED K50 = 0.5 for feedback promoter:
//            double ratio_H_tau2 = pow(H_tau/0.5, hn);

            //use generated KH values from lognormal distributions (set on cell birth):
//            ratio_H_tau = pow(H_tau/lnKH, hn);
            //for feedback circuit:
//            double ratio_H_tau2 = pow(H_tau/lnKH2, hn);

            double fbalpha = double(eQ::parameters["MODULUS_FEEDBACK_FRACTION"]);//convex combination version
            double fbscale = double(eQ::parameters["MODULUS_FEEDBACK_STRENGTH"]);//multiple for feedback production

            //WITH DIFFUSION:
            if( (eQ::parameters["MODULUS_option"] == "+D+F") || (eQ::parameters["MODULUS_option"] == "+D-F") )
            {
                if(eQ::parameters["MODULUS_option"] == "+D-F")
                {//NO FEEDBACK:
                    //C4 synthase: [S] \in [0, 1]
                    delta[S] = dt * ( //rhlI under PL-lac promoter (direct activation by "c"=iptg, ratio set explicitly above)
                             (gamma_d * ratioPlac)/(1.0  + ratioPlac));
                }
                else
                {//WITH FEEDBACK:
                    //rhlI under PL-lac promoter (direct activation by "c"=iptg, ratio set explicitly above)
                    delta[S] = dt * (
                           gamma_d * ratioPlac/(1.0  + ratioPlac)
                        +  fbscale * gamma_d * ratio_H_tau/(1.0 + ratio_H_tau));
//                    (1.0 - fbalpha)*gamma_d * ratioPlac/(1.0  + ratioPlac)
//                 +  fbalpha * gamma_d * ratio_H_tau/(1.0 + ratio_H_tau));
                }

                //READOUT:
                delta[FP] = dt * ( //Prhl promoter
                           gamma_d * (ratio_H_tau)/(1.0 + ratio_H_tau));

//                dC4HSL = dt * (//positive direction is outwards from cell
//                        params->DH * (conc[H] - eHSL[0]));    //MEMBRANE DIFFUSION

                delta[H] = dt * (double(eQ::parameters["hslProductionRate_C4"]) * conc[S])//PRODUCTION
                        - dHSL[0];                                               //MEMBRANE DIFFUSION
//                        - dC4HSL;                                               //MEMBRANE DIFFUSION
            }
            //NO DIFFUSION: use AraC
            else if( (eQ::parameters["MODULUS_option"] == "-D+F") || (eQ::parameters["MODULUS_option"] == "-D-F") )
            {
                //for recording only:
                delta[L] = dt * ( gamma_d * iptg );

                //AraC production:
                if(eQ::parameters["MODULUS_option"] == "-D-F")
                {//NO FEEDBACK:
                    delta[S] = dt * ( //S=AraC under PL-lac promoter (H is identically S here)
                             (gamma_d * ratioPlac)/(1.0  + ratioPlac));
                }
                //with positive feedback, add AraC under Pbad promoter
                else if(eQ::parameters["MODULUS_option"] == "-D+F")
                {//WITH FEEDBACK:
                    delta[S] = dt * (
                            gamma_d * ratioPlac/(1.0  + ratioPlac)
                         +  fbscale * gamma_d * ratio_H_tau/(1.0 + ratio_H_tau));
//                           (1.0 - fbalpha)*gamma_d * ratioPlac/(1.0  + ratioPlac)
//                        +  fbalpha * gamma_d * ratio_H_tau/(1.0 + ratio_H_tau));
                }
                //READOUT
                delta[FP] = dt * ( //Pbad promoter for AraC
                    //normalize by the expected division time (to cancel dilution due to growth)
                    gamma_d * (ratio_H_tau)/(1.0 + ratio_H_tau));

                //MEMBRANE DIFFUSION OF HSL (degradation via spatial diffusion is done on grid)
                dHSL[0] = 0.0;    //MEMBRANE DIFFUSION
//                dC4HSL = 0.0;    //MEMBRANE DIFFUSION
                //H=AraC here:
                double gammaT = double(eQ::parameters["gammaT_C4"]);//need to scale back out the HSL decay rate:
                delta[H] = dt * (gamma_d * double(eQ::parameters["hslProductionRate_C4"])/gammaT * conc[S]);//PRODUCTION
            }
            else ;//bad option
        }
        else if("MODULUS_2" == eQ::parameters["simType"])
        {
            //convert signal "c" to fixed activating signal for Plac
            //passed-in value set in eQabm.cpp:
            iptg = double(eQ::parameters["MODULUS_IPTG"]);
            double Klac = 30.0; //nanoMolar Plac K50 via Lutz et al.
            Klac = 60.0;
            double ratioPlac = pow(iptg/Klac, hn);
            //FIXED K50 for output promoter:
//            double HK = double(eQ::parameters["lnmean"]) * 8.0;
//            double HK = 1.8e4;
            double HK = 3.0e4;
            ratio_H_tau = pow(tHSL[0]/HK, hn);

            //C4 synthase:
            delta[S] = dt * ( //rhlI under PL-lac promoter (direct activation by "c"=iptg, ratio set explicitly above)
                     (gamma_d * ratioPlac)/(1.0  + ratioPlac));

            deltaHSL[0] = dt * (double(eQ::parameters["hslProductionRate_C4"]) * conc[S])//PRODUCTION
                       - dHSL[0];                                               //MEMBRANE DIFFUSION

            deltaPROTEIN[tetR] = dt * ((gamma_d * ratio_H_tau)/(1.0  + ratio_H_tau));
            double ratio_tetR = pow(tPROTEIN[tetR]/0.5, hn);

            //READOUT:
            delta[FP] = dt * gamma_d * (
                    ratioPlac/(1.0  + ratioPlac)//LOCAL EXCITATION
//                    +    1.0/(1.0 + ratio_tetR)     //GLOBAL INHIBITION
                    +    1.0/(1.0 + ratio_H_tau)     //GLOBAL INHIBITION
            );

            conc[L] = ratioPlac;//use as temp storage for testing value
            conc[H] = ratio_H_tau;
            conc[LEGI_A] = ratio_tetR;
        }

        for(int i=0;i<numConcentrations;i++)
        {
            conc[i] += delta[i];
            if (conc[i] < 0.0)  conc[i] = 0.0;
            else conc[i] /= nanoMolarPerMolecule;//store as protein #
        }
        for(size_t i(0); i< iHSL.size(); i++)
        {
            iHSL[i] += deltaHSL[i];
            if (iHSL[i] < 0.0)
            {
                iHSL[i] = 0.0;
            }
            HSL_tau[i].push(iHSL[i]);//store as concentration
            iHSL[i] /= nanoMolarPerMolecule;//store as protein #
        }
        for(size_t i(0); i< iPROTEIN.size(); i++)
        {
            iPROTEIN[i] += deltaPROTEIN[i];
            if (iPROTEIN[i] < 0.0) iPROTEIN[i] = 0.0;
            PROTEIN_tau[i].push(iPROTEIN[i]);//store as concentration
            iPROTEIN[i] /= nanoMolarPerMolecule;//store as protein #
        }
            //push new HSL and LacI values to queue: [STORE AS CONCENTRATION in nM]:
//            qH_tau.push(conc[H] * nanoMolarPerMolecule);
//            qL_tau.push(conc[L] * nanoMolarPerMolecule);
//            qI_tau.push(conc[I] * nanoMolarPerMolecule);

            pushCellLengthValue(lengthMicrons);

        return dHSL;
}

void Strain::solveLinearSystem(double &H_i, double &H_e, double y1, double y2, double cellFraction)
{
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
/*
//MEMBRANE DIFFUSION OF HSL (degradation via spatial diffusion is done on grid)
//                double newHi, newHe;
//                solveLinearSystem(newHi, newHe, dt*double(eQ::parameters["hslProductionRate_C4"]) * conc[S] + conc[H], C4ext, rho_cell);
//                delta[H] = newHi - conc[H];
//                dC4HSL = newHe  - C4ext;
//{
//    //data:  conc[i] stores protein # for each protein between timesteps (covert to conc. first)
////    double explicitDilution = params->dilution;
//    double explicitDilution         = 0.0;//set to zero when using ABM  with cell growth
//    //we store protein count, convert to nanoMolar concentration:
////    double volumeScaling = lengthMicrons*cellVolumePerLength*0.602;
//    double oneProteinToNanoMolar = eQ::proteinNumberToNanoMolar(1.0, lengthMicrons);

//    if(lengthMicrons == 0.0)
//    {
//        std::cout<<"ERROR: scale for computeProteins() set to 0!"<<std::endl;
//        return std::vector<double> {};
//    }
//    for(int i=0; i<numConcentrations; i++)
//    {
//        conc[i] *= oneProteinToNanoMolar;//convert to proper concentration (units: nM)
//    }
//    //reactions on promoters use delayed concentrations (stored as t-tau concentrations, not #, thus no scaling):
//    H_tau = qH_tau.front();  qH_tau.pop();
//    L_tau = qL_tau.front();  qL_tau.pop();
//    I_tau = qI_tau.front();  qI_tau.pop();
//    //compute Hill function ratios once:
//    ratio_H_tau = pow(H_tau/params->K.H, params->nH);
//    ratio_L_tau = pow(L_tau/params->K.L, params->nL);
//    ratio_I_tau = pow(I_tau/params->K.I, params->nI);

//    totalProtein =
//            params->K.C + conc[S]+ conc[L]+ conc[A]+ conc[FP]+ conc[M];
//    totalHSL =
//            params->K.A + conc[H]+ conc[I];
//    degradationProtein =
//            (params->dC / totalProtein) + explicitDilution;//explicitly model dilution here
//    degradationHSL =
//            (params->dA / totalHSL) * conc[A] + explicitDilution;//explicitly model dilution here

//    //AiiA
//    delta[A] = dt * ( //aiiA has cin promoter
//              (params->eta.A0)
////              (params->eta.A0 + params->eta.A1 * ratio_I_tau)/(1.0 + ratio_I_tau)
//            - degradationProtein * conc[A] );
//    //Mature FP
//    delta[M] = dt * (
//            params->m * conc[FP]
//            - degradationProtein * conc[M] );

//    //MEMBRANE DIFFUSION OF HSL (degradation via spatial diffusion is done on grid)
//    dC4HSL = dt * (
//            params->DH * (conc[H] - C4ext));    //MEMBRANE DIFFUSION
//    dC14HSL = dt * (
//            params->DI * (conc[I] - C14ext));   //MEMBRANE DIFFUSION

//    if(eQ::strainType::ACTIVATOR == whichType)
//    {
//        //hack-in a pressure value set by model, use as Hill function parameter:
////        ratio_I_tau = pow(pressureValue/pressureK50, 2.0);
//        //LacI
//        delta[L] = dt * ( //lacI has cin promoter
//                (params->eta.L0 + params->eta.L1 * ratio_I_tau)/(1.0 + ratio_I_tau)
//                - degradationProtein * conc[L] );
//        //C4 synthase:
//        delta[S] = dt * ( //rhlI or cinI both have rhl/lac hybrid promoter
//                (params->eta.S0 + params->eta.S1 * ratio_H_tau)/(1.0 + ratio_H_tau + ratio_L_tau)
//                - degradationProtein * conc[S] );
//        //CFP tracks  C4HSL synthase:
//        delta[FP] = 0.0;
//        delta[H] = dt * (
//                params->phi * conc[S]									//PRODUCTION
//                - degradationHSL * conc[H] )							//DEGRADATION+DILUTION
//                - dC4HSL;                                               //MEMBRANE DIFFUSION
//        delta[I] = 0.0;//no repressor
//    }
//    else if(eQ::strainType::REPRESSOR == whichType)
//    {
//        //C14 synthase:
//        delta[S] = 0.0;
//        //YFP tracks C14HSL synthase:
//        delta[FP] = dt * ( //repressor has cin/lac promoter
//                (params->eta.FP0 + params->eta.FP1 * ratio_H_tau)/(1.0 + ratio_H_tau)
////                    (params->eta.FP0 + params->eta.FP1 * ratio_I_tau)/(1.0 + ratio_I_tau + ratio_L_tau)
//                - (degradationProtein + params->m) * conc[FP] );
//        delta[I] = 0.0;
//        delta[H] = dt * (- degradationHSL * conc[H] )                   //DEGRADATION+DILUTION
//                - dC4HSL;                                               //MEMBRANE DIFFUSION
//    }

//    for(int i=0;i<numConcentrations;i++)
//    {
//        conc[i] += delta[i];
//        if (conc[i] < 0.0)  conc[i] = 0.0;
//        else conc[i] /= oneProteinToNanoMolar;//store as protein #
//    }
//    //push new HSL and LacI values to queue: [STORE AS CONCENTRATION in nM]:
//    qH_tau.push(conc[H] * oneProteinToNanoMolar);
//    qL_tau.push(conc[L] * oneProteinToNanoMolar);
//    qI_tau.push(conc[I] * oneProteinToNanoMolar);

//    return std::vector<double> {dC4HSL/oneProteinToNanoMolar, dC14HSL/oneProteinToNanoMolar};//return HSL # changes
//}

*/
////////////////////////////////////////////////////////////////////////////////
//				loadParams():
////////////////////////////////////////////////////////////////////////////////
void loadParams(eQ::strainType which, struct dso_parameters &params, struct rates &rates)
{//SELECT PROMOTER STRENGTHS HERE:
	if(eQ::strainType::ACTIVATOR == which)
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
	else if(eQ::strainType::REPRESSOR == which)
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


/*

    if(false){
    //Aspect ratio gene circuit:
    //ACTIVATOR produces C4HSL under induction, diffuses out of cell
    //REPRESSOR diffuses in C4 and produces ftsZ
    //snip: old code
    if(eQ::strainType::ACTIVATOR == whichType)
    {
        //CFP <-- C4HSL synthase:
        delta[FP] = dt * (
                (params->eta.FP0 + params->eta.FP1 * ratio_H_tau)/(1.0 + ratio_H_tau + ratio_L_tau)
                - (degradationProtein + params->m) * conc[FP] );

        double hsl_ext,hsl_int;
        solveLinearSystem(hsl_int, hsl_ext, conc[H] + dt*params->phi * conc[S], C4ext);
        delta[H] = hsl_int - conc[H];
        dC4HSL = hsl_ext - C4ext;

        solveLinearSystem(hsl_int, hsl_ext, conc[I], C14ext);
        delta[I] = hsl_int - conc[I];
        dC14HSL = hsl_ext - C14ext;

        // delta[H] = dt * (
        // 		params->phi * conc[S]									//PRODUCTION
        // 		- params->DH * (conc[H] - C4ext)						//MEMBRANE DIFFUSION
        // 		- degradationHSL * conc[H] );							//DEGRADATION+DILUTION
        // delta[I] = dt * (
        // 		- params->DI * (conc[I] - C14ext)						//MEMBRANE DIFFUSION
        // 		- degradationHSL * conc[I] );							//DEGRADATION+DILUTION
    }
    else if(eQ::strainType::REPRESSOR == whichType)
    {
        //YFP <-- C14HSL synthase:
        delta[FP] = dt * (
                (params->eta.FP0 + params->eta.FP1 * ratio_I_tau)/(1.0 + ratio_I_tau + ratio_L_tau)
                - (degradationProtein + params->m) * conc[FP] );
        double hsl_ext,hsl_int;

        solveLinearSystem(hsl_int, hsl_ext, conc[I] + dt*params->phi * conc[S], C14ext);
        delta[I] = hsl_int - conc[I];
        dC14HSL = hsl_ext - C14ext;
            solveLinearSystem(hsl_int, hsl_ext, conc[H], C4ext);
            delta[H] = hsl_int - conc[H];
            dC4HSL = hsl_ext - C4ext;
        // delta[I] = dt * (
        // 		params->phi * conc[S]									//PRODUCTION
        // 		- params->DI * (conc[I] - C14ext)						//MEMBRANE DIFFUSION
        // 		- degradationHSL * conc[I] );							//DEGRADATION+DILUTION
        // delta[H] = dt * (
        // 		- params->DH * (conc[H] - C4ext)						//MEMBRANE DIFFUSION
        // 		- degradationHSL * conc[H] );							//DEGRADATION+DILUTION
    }

    // //MEMBRANE DIFFUSION OF HSL (degradation is done on grid)
    // dC4HSL = dt * (
    // 		// params->DH * (conc[H] - C4ext));
    // 		params->DH * (conc[H] - C4ext*scale));
    dC4HSL *= scale;
    C4ext += dC4HSL;
        if(C4ext < 0.0) C4ext = 0.0;
    // dC14HSL = dt * (
    // 		// params->DH * (conc[I] - C14ext));
    // 		params->DH * (conc[I] - C14ext*scale));
    dC14HSL *= scale;
    C14ext += dC14HSL;
        if(C14ext < 0.0) C14ext = 0.0;

    //set new concentrations for use in diffusion:
    for(int i=0;i<numConcentrations;i++)
    {
        conc[i] += delta[i];
        conc[i] *= scale;//store as protein #
        if (conc[i] < 0.0)  conc[i] = 0.0;
    }
    //push new HSL and LacI values to queue:
    qH_tau.push(conc[H]);//push as protein #
    qL_tau.push(conc[L]);//push as protein #
    qI_tau.push(conc[I]);//push as protein #
  }
*/
