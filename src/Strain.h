#ifndef STRAIN_H
#define STRAIN_H

#include <iostream>
#include <array>
#include <queue>
#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <memory>

#include "eQ.h"
#include "eQcell.h"
#include "eQrandom.h"

//PARAMETER SET FOR PRODUCTION RATES = "eta" (BASAL, ACTIVE)
enum etas
{
	R0,R1,
	C0,C1,
	F0,F1,
	Y0,Y1,
	L0,L1,
	A0,A1,
	numEtas
};
enum promoterStrengths
{
	WEAK, MEDIUM, STRONG, LAC,
	numStrengths,
	ONLY=0
};
enum basalRates
{//Indeterminant parameter values (SR, SC, SF, SY, SL, SA, and ClpXP)
	SR, SC, SL, SA, SF, SY, ClpXP, 
	numBasalRates
};


using etas_t = struct eta
{
	double 
	S0,S1,
	FP0,FP1,
	L0,L1,
	A0,A1;
};
using K50_t = struct K50
{
	double H,L,I,C,A;
};


struct rates//master table of values published in Chen et al (2015)
{
    double eta
        [etas::numEtas]
        [promoterStrengths::numStrengths];
    K50_t K
        [promoterStrengths::numStrengths];
    double nH,nL,nI;					//Hill coeffs.
	double dC,dA,dilution,mu_e;				//degradation, dilution constants
	double DH, DI, phiH, phiI, m;			//membrane diffusion, production, maturation rates
};

//data struct to hold actual populated values of all constants per-strain
struct dso_parameters
{
	etas_t eta;		//production rates structure
	K50_t K;		//x50 (half-max values) structure
	double dC,dA,dilution,mu_e;		//degradation, dilution constants
	double DH, DI, phi, m;			//diffusion, production, maturation rates
    double nH,nL,nI;			//Hill coeffs.
};

////////////////////////////////////////////////////////////////////////////////
//							class Strain
////////////////////////////////////////////////////////////////////////////////
class Strain
{
public:
    //per-cell values
    //S = "synthase" activator or repressor dep. on strain type
    //FP = "fluorescent protein" CFP or YFP dep. on strain type
    enum concentrations
    {
        S, L, A, FP, M,
        H_T, HR, R_T, R,
        H, I, LEGI_C, LEGI_A, LEGI_R,
        AraC, tetR,
        CFP,YFP,RFP,GFP,
        NUM_CONCENTRATIONS
    };
    enum hsl
    {
        C4, C14,
        NUM_HSL
    };


    struct Params
    {
        eQ::Cell::strainType	whichType;
        double					dt;
        double					nodesPerMicronScale;
        size_t					numHSL;
        double					promoterDelayTimeMins;
        eQ::Cell::Params *      baseData;//populated on cell creation; NB: cannot be a reference
    };
    Strain::Params				params;

    //virtual destructor needed in base class so that derived classes can be deleted properly:
    virtual ~Strain() = default;

	//clone will call default copy constructor (duplicates all objects by value)
    virtual std::shared_ptr<Strain>
        clone()  {return std::make_shared<Strain>(*this);}

    virtual void    init(){}
    virtual double  growthRateScaling(){return 1.0;}

    //Note: must declare methods virtual if you want to over-ride them in subclasses
    virtual double	getProteinNumber(concentrations which)	const {return conc[which];}
    virtual double	getHSL(hsl which)						const {return iHSL[which];}
    virtual double	getDelayedHSL(hsl which)				const {return tHSL[which];}
    virtual const std::vector<double>&	getHSLvec()			const {return iHSL;}

    eQ::Cell::strainType  getStrainType(void)               const {return params.whichType;}
    bool			inductionActive(void)                   const {return inductionFlag;}
    void			induce()                                {inductionFlag = true;}
    void			setPressureValue(double pressure)       {pressureValue = pressure;}
    void			setPressureK50(double pressure)         {pressureK50 = pressure;}


  protected:



    //ALL VARIABLES WILL BE COPIED BY THE DEFAULT COPY CONSTRUCTOR:
    size_t	queueDepth;
	bool 	inductionFlag;
	double  pressureValue, pressureK50;

    //THE FOLLOWING MUST BE RE-INITIALIZED ON DIVISION OF DAUGHTER CELL:
    double								conc[NUM_CONCENTRATIONS];
    std::vector<double>                 iHSL;
    std::vector<double>                 tHSL;
    std::vector<std::queue<double>>     HSL_tau;

    std::vector<double>                 iPROTEIN;
    std::vector<double>                 tPROTEIN;
    std::vector<std::queue<double>>     PROTEIN_tau;
    //THE FOLLOWING NEED NOT BE RE-INITIALIZED:
	std::vector<double>                 dHSL;//membrane diffusion value
	std::vector<double>                 deltaHSL;//per-timestep production value
	std::vector<double>                 deltaPROTEIN;//per-timestep production value

	
	//THE FOLLOWING NEED NOT BE RE-INITIALIZED:
    double	ratio_H_tau, ratio_L_tau, ratio_I_tau;
	double 	totalProtein, totalHSL;
    double 	delta[NUM_CONCENTRATIONS];
	double  dC4HSL, dC14HSL;
	double  degradationProtein, degradationHSL;
    double	nanoMolarPerMolecule;

    double	R_total;
    double	H_tau, L_tau, I_tau;

    friend std::ostream &operator<<(std::ostream &os, const eQ::Cell::strainType &type);

public:
    //Strain::Strain(const eQ::Cell::strainType strainType, const dso_parameters *p, const double timestep, const double nodesPerMicronScale)
    //    : whichType(strainType), params(p), dt(timestep), nodesPerMicron(nodesPerMicronScale), inductionFlag(false)
    Strain(const Strain::Params &p)
        : params(p)
    {

        //NOTE:  valid for fixed dt/delay times only
        queueDepth = size_t(ceil(params.promoterDelayTimeMins/params.dt));

    //    //old data structures:
    //	for(size_t i=0; i<queueDepth; i++)
    //	{
    //		qH_tau.push(0.0);
    //		qL_tau.push(0.0);
    //		qI_tau.push(0.0);
    //	}
        for(int i=0; i<NUM_CONCENTRATIONS; i++)
        {
            conc[i]=0.0;//don't omit this!
            delta[i]=0.0;
        }
    }
    void initializeDataStructures(const size_t numHSL, const size_t numProteins)
    {
        iHSL.assign(numHSL, 0.0);
        tHSL.assign(numHSL, 0.0);//per timestep value of popped queue value for convenience
        dHSL.assign(numHSL, 0.0);
        deltaHSL.assign(numHSL, 0.0);

        iPROTEIN.assign(numProteins, 0.0);
        tPROTEIN.assign(numProteins, 0.0);//per timestep value of popped queue value for convenience
        deltaPROTEIN.assign(numProteins, 0.0);

        std::queue<double> q;
        for(size_t j=0; j<queueDepth; j++) {q.push(0.0);}

        HSL_tau.assign(numHSL, q);
        PROTEIN_tau.assign(numProteins, q);

    }
    virtual void divideProteins(std::shared_ptr<Strain> daughter, double fracToDaughter)
    {
        //NOTE: the size of all data structures is determined by constructor of derived class
        double protein;
        double rd = fracToDaughter;
        double rm = (1.0 - fracToDaughter);

        for(size_t i=0;i<NUM_CONCENTRATIONS;i++)
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


    virtual std::vector<double>
    computeProteins
        (const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
    {
        computeConcentrations(eHSL, membraneRate, lengthMicrons);

        //        //compute Hill function ratios once:
    //    ratio_H_tau = pow(H_tau/params->K.H, params->nH);
    //    ratio_L_tau = pow(L_tau/params->K.L, params->nL);
    //    ratio_I_tau = pow(I_tau/params->K.I, params->nI);
        return std::vector<double> {0,0};//return HSL # changes

    }


    void computeConcentrations(const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons)
    {
        //COMMON CODE TO INITIALIZE DATA, COMPUTE CONVERSIONS FROM # TO CONC., POP DELAY QUEUES, COMPUTE MEMBRANE FLUX:

        if(lengthMicrons == 0.0)
        {
            std::cout<<"ERROR: scale for computeProteins() set to 0!"<<std::endl;
            return;
        }
        nanoMolarPerMolecule = eQ::Cell::moleculeNumberToNanoMolar(1.0, lengthMicrons);

        //OLD PROTEIN ARRAYS:
        for(int i=0; i<NUM_CONCENTRATIONS; i++)
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
            dHSL[i] = params.dt * (
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
    void pushConcentrations()
    {
        //OLD DATA STRUCTURES:
        for(int i=0;i<NUM_CONCENTRATIONS;i++)
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
};

class sendRecvStrain : public Strain
{
private:
    enum  qProteins
    {
        FP = 0,
        NUM_QUEUES
    };
public:
    enum hslTypes
    {
        C4HSL   = 0,
        C14HSL  = 1,
        NUM_HSLTYPES
    };
    enum inductionFlag
    {
        INDUCTION,
        NUM_INDUCTIONFLAGS
    };

    static std::vector<bool> inductionFlags;

    static void setFlag(inductionFlag which) {inductionFlags[which] = true;}

    //CONSTRUCTOR: (called for initial seed cells only)
//	sendRecvStrain(eQ::Cell::strainType strainType, const double timestep, const double nodesPerMicron, const size_t numHSL)
//        : Strain(strainType, nullptr, timestep, nodesPerMicron)
    sendRecvStrain(const Strain::Params &p)
        : Strain(p)
    {
        initializeDataStructures(params.numHSL, qProteins::NUM_QUEUES);
        inductionFlags.assign(NUM_INDUCTIONFLAGS, false);
    }
    //calls default copy constructor for derived class:
    std::shared_ptr<Strain>
        clone()  override {return std::make_shared<sendRecvStrain>(*this);}  // requires C++ 14
    std::vector<double>
        computeProteins(
            const std::vector<double> &eHSL, const std::vector<double> &membraneD, const double lengthMicrons) override;
};
class aspectRatioInvasionStrain : public Strain
{
private:
	enum  qProteins
	{
		ftsZ = 0,
		NUM_QUEUES
	};

public:
    enum hslType
    {
        C4HSL   = 0,
        C14HSL  = 1,
        NUM_HSLTYPES
    };
    enum inductionFlag
    {
        ASPECTRATIO_INDUCTION,
        NUM_INDUCTIONFLAGS
    };

    static std::vector<bool> inductionFlags;

    static void setFlag(inductionFlag which) {inductionFlags[which] = true;}

    //CONSTRUCTOR: (called for initial seed cells only)
//	aspectRatioInvasionStrain(eQ::Cell::strainType strainType, const double timestep, const double nodesPerMicronScale, const size_t numHSL)
//        : Strain(strainType, nullptr, timestep, nodesPerMicronScale)
//    {
//        initializeDataStructures(numHSL, qProteins::NUM_QUEUES);
//    }
    aspectRatioInvasionStrain(const Strain::Params &p) : Strain(p)
    {
        initializeDataStructures(params.numHSL, qProteins::NUM_QUEUES);
        inductionFlags.assign(NUM_INDUCTIONFLAGS, false);
    }
    //calls default copy constructor for derived class:
    std::shared_ptr<Strain>
        clone()  override {return std::make_shared<aspectRatioInvasionStrain>(*this);}  // requires C++ 14

    std::vector<double>
        computeProteins(const std::vector<double> &eHSL, const std::vector<double> &membraneD, const double lengthMicrons) override;

};
class aspectRatioOscillator : public Strain
{
private:
    enum  qProteins
    {
        ftsZ = 0,
        NUM_QUEUES
    };
public:
    enum hslType
    {
        C4HSL   = 0,
        C14HSL  = 1,
        NUM_HSLTYPES
    };
    enum inductionFlag
    {
        ASPECTRATIO_INDUCTION,
        NUM_INDUCTIONFLAGS
    };
    static std::vector<bool> inductionFlags;
    static void setFlag(inductionFlag which) {inductionFlags[which] = true;}

    aspectRatioOscillator(const Strain::Params &p) : Strain(p)
    {
        initializeDataStructures(params.numHSL, qProteins::NUM_QUEUES);
        inductionFlags.assign(NUM_INDUCTIONFLAGS, false);
    }
    //calls default copy constructor for derived class:
    std::shared_ptr<Strain>
        clone()  override {return std::make_shared<aspectRatioOscillator>(*this);}  // requires C++ 14
    std::vector<double>
        computeProteins(const std::vector<double> &eHSL, const std::vector<double> &membraneD, const double lengthMicrons) override;
};
class parB_MotherStrain : public Strain
{
private:
    enum  qProteins
    {
        _tetR = 0,
        NUM_QUEUES
    };
public:
    enum hslType
    {
        C4HSL   = 0,
        C14HSL  = 1,
        NUM_HSLTYPES
    };
    enum inductionFlag
    {
        INDUCTION,
        NUM_INDUCTIONFLAGS
    };
    static std::vector<bool> inductionFlags;
    static void setFlag(inductionFlag which) {inductionFlags[which] = true;}

    parB_MotherStrain(const Strain::Params &p) : Strain(p)
    {
        initializeDataStructures(params.numHSL, qProteins::NUM_QUEUES);
        inductionFlags.assign(NUM_INDUCTIONFLAGS, false);
        std::random_device rdev;
        rn = std::make_shared<eQ::uniformRandomNumber>(rdev);
    }
    //calls default copy constructor for derived class:
    std::shared_ptr<Strain> clone()  override  // requires C++ 14
    {
        auto clone = std::make_shared<parB_MotherStrain>(*this);
        if(parB_losePlasmid)
        {
            clone->parB_losePlasmid = false;//reset for daughter cell
            clone->params.whichType = eQ::Cell::strainType::REPRESSOR;
//            if(rn->successCoinFlip())   clone->params.whichType = eQ::Cell::strainType::REPRESSOR;
//            else                        params.whichType = eQ::Cell::strainType::REPRESSOR;
        }
        return clone;
    }
    void init () override
    {//called after seeding creation to set params dependent on the base data (otherwise invalid pointer)
        tetR_productionRate = log(2)/params.baseData->doublingPeriodMinutes;
    }
    std::vector<double>
        computeProteins(const std::vector<double> &eHSL, const std::vector<double> &membraneD, const double lengthMicrons) override;

    double growthRateScaling() override;
private:
    bool parB_losePlasmid = false;
    std::shared_ptr<eQ::uniformRandomNumber> rn;
    double tetR_productionRate;
};
class MODULUSmodule : public Strain
{
private:
    enum  qProteins
    {
        tetR = 0,
        NUM_QUEUES
    };
public:
    enum hslTypes
    {
        C4HSL   = 0,
        C14HSL  = 1,
        NUM_HSLTYPES
    };
    enum inductionFlags
    {
        INDUCTION,
        NUM_INDUCTIONFLAGS
    };

    static std::vector<bool> inductionFlags;

    //CONSTRUCTOR: (called for initial seed cells only)
//    MODULUSmodule(const double timestep, const double nodesPerMicron, const size_t numHSL)
//		: Strain(eQ::Cell::strainType::ACTIVATOR, nullptr, timestep, nodesPerMicron)
    MODULUSmodule(const Strain::Params &p) : Strain(p)
    {
        initializeDataStructures(params.numHSL, qProteins::NUM_QUEUES);
        inductionFlags.assign(NUM_INDUCTIONFLAGS, false);

        //sample the K50 parameters for promoters
        //note:  will be reset for daughter cells
//        setKH();

//        //create a queue for time-averaging the cell length to reduce noise:
//        double minsQ = double(eQ::parameters["MODULUS_TIME_AVERAGE_MINS"]);
//        auto qdepth = size_t(ceil(minsQ/timestep));
//        for(size_t i=0; i<qdepth; i++)
//        {
//            qReadout.push_front(0.0);
//        }
//        debugCounter=0;

    }
    //calls default copy constructor for derived class:
    std::shared_ptr<Strain>
        clone()  override {return std::make_shared<MODULUSmodule>(*this);}  // requires C++ 14
    std::vector<double>
        computeProteins(const std::vector<double> &eHSL, const std::vector<double> &membraneD, const double lengthMicrons) override;


////    void divideProteins(std::shared_ptr<MODULUSmodule> daughter, double fracToDaughter)
//    void divideProteins(std::shared_ptr<Strain> pdaughter, double fracToDaughter) override
//    {
//        //must cast the base-class pointer parameter to most easily override divideProteins() virtual in base class:
//        auto daughter = std::dynamic_pointer_cast<MODULUSmodule>(pdaughter);

////        //copy parent values to daughter before re-sampling:
////        daughter->lnKL = lnKL;
////		daughter->lnKH = lnKH;
////		daughter->lnKH2 = lnKH2;
////        reSampleK50s();
////        daughter->reSampleK50s();

//		Strain::divideProteins(daughter, fracToDaughter);

////        //divide cell length averaging queue
////        double rm = 1.0-fracToDaughter;
////        for(size_t i=0;i<qReadout.size();i++)
////        {//DIVIDE BY FRACTION
////            auto data = qReadout.back();
////                qReadout.push_front(rm*data);
////                qReadout.pop_back();
////					daughter->qReadout.push_front(fracToDaughter*data);
////					daughter->qReadout.pop_back();
////        }
////        computeAverageCellLength();
//	}

    //	void setKH()
//	{//for initial seeds:
//		//TODO: verify log normal generator was initialized
//		lnKL =   logNormals.getLogNormalRandomNumber(0);//PL_lac promoter K50
//		lnKH =   logNormals.getLogNormalRandomNumber(1);//Prhl or Pbad promoter K50
//		lnKH2 =   logNormals.getLogNormalRandomNumber(1);//Prhl or Pbad promoter K50
//	}
//    void reSampleK50s()
//	{//for correlating values parent/daughter:
//		//TODO: verify log normal generator was initialized
//		auto correlateValue = [](double parent, double sample)
//		{
//            const double a = double(eQ::parameters["K50_correlationScale"]);//relative geometric weight of parent value
////            const double a = 0.9;//relative geometric weight of parent value
//            return pow(parent, a)*pow(sample, 1-a);
//		};
//		lnKL =  correlateValue(lnKL, logNormals.getLogNormalRandomNumber(0));//PL_lac promoter K50
//		lnKH =  correlateValue(lnKH, logNormals.getLogNormalRandomNumber(1));//Prhl or Pbad promoter K50
//		lnKH2 =  correlateValue(lnKH2, logNormals.getLogNormalRandomNumber(1));//Prhl or Pbad promoter K50
//	}

//	eQ::logNormalRandomNumber logNormals;
    double lnKH, lnKL, lnKH2;

//    std::deque<double> qReadout;
//    double filteredCellLength;

//    double computeAverageCellLength()
//    {
//        double sum = std::accumulate(qReadout.begin(), qReadout.end(), 0.0);
//        filteredCellLength = sum / qReadout.size();
//        return filteredCellLength;
//    }
//    void pushCellLengthValue(double value)
//    {
//        qReadout.push_front(value);
//        qReadout.pop_back();
//        computeAverageCellLength();
//    }

};

//void initRates(double *, struct rates &);
//void loadParams(eQ::Cell::strainType which, struct dso_parameters &params, struct rates &rates);
	

#endif
