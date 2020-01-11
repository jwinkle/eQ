#ifndef STRAIN_H
#define STRAIN_H

#include <iostream>
#include <array>
#include <queue>
#include <vector>
#include <string>
#include <iostream>
#include <math.h>

#include "eQ.h"

// #define TAU_DELAY 	7500 //7.5 minutes -> timestep=0.001
#define TAU_DELAY 	7.5 //7.5 minutes -> timestep=0.001

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
//per-cell values 
//S = "synthase" activator or repressor dep. on strain type
//FP = "fluorescent protein" CFP or YFP dep. on strain type
enum concentrations
{
    S, L, A, FP, M,
    H_T, HR, R_T, R,
    H, I, LEGI_C, LEGI_A, LEGI_R,
    AraC,
	numConcentrations
};


typedef struct etas_t
{
	double 
	S0,S1,
	FP0,FP1,
	L0,L1,
	A0,A1;
}etas_t;
typedef struct K50
{
	double H,L,I,C,A;
}K50_t;


struct rates//master table of values published in Chen et al (2015)
{
	double eta
		[static_cast<size_t>(etas::numEtas)]
		[static_cast<size_t>(promoterStrengths::numStrengths)];
	K50_t K
		[static_cast<size_t>(promoterStrengths::numStrengths)];

//    unsigned int nH,nL,nI;					//Hill coeffs.
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
//    unsigned int nH,nL,nI;			//Hill coeffs.
    double nH,nL,nI;			//Hill coeffs.
};

////////////////////////////////////////////////////////////////////////////////
//							class Strain
////////////////////////////////////////////////////////////////////////////////
class Strain
{
  public:
	Strain(eQ::strainType strainType, const struct dso_parameters *p, const double timestep, const double nodesPerMicronScale);

	virtual ~Strain();

	//clone will call default copy constructor (duplicates all objects by value)
	virtual std::shared_ptr<Strain> clone() const {return std::make_shared<Strain>(*this);}

    //Note: must declare methods virtual if you want to over-ride them in subclasses
    virtual std::vector<double>
		computeProteins( double C4ext, double C14ext, double scale);
    virtual std::vector<double>
        computeProteins(const std::vector<double> &eHSL, const std::vector<double> &membraneD, const double lengthMicrons);

	virtual void
		divideProteins(std::shared_ptr<Strain> daughter, double fracToDaughter);

    void computeConcentrations(const std::vector<double> &eHSL, const std::vector<double> &membraneRate, const double lengthMicrons);
    void pushConcentrations();

	eQ::strainType  getStrainType(void){return whichType;}
	void			loadParams(eQ::strainType which);
	void			induce(){inductionFlag = true;}
	bool			inductionActive(void){return inductionFlag;}
	void			setPressureValue(double pressure) {pressureValue = pressure;}
	void			setPressureK50(double pressure) {pressureK50 = pressure;}
	double			getProteinNumber(enum concentrations which) {return conc[which];}
	void			solveLinearSystem(double &H_i, double &H_e, double y1, double y2, double rho_cell);

	//ALL VARIABLES WILL WILL BE COPIED BY THE DEFAULT COPY CONSTRUCTOR:
	eQ::strainType					whichType;
	size_t							queueDepth;
	const struct dso_parameters 	*params;

	//THE FOLLOWING MUST BE RE-INITIALIZED ON DIVISION OF DAUGHTER CELL:
	double 	conc[numConcentrations];
	std::vector<double> iHSL;
	std::vector<double> tHSL;
	std::vector<std::queue<double>> HSL_tau;

	std::vector<double> iPROTEIN;
	std::vector<double> tPROTEIN;
	std::vector<std::queue<double>> PROTEIN_tau;

	//THE FOLLOWING NEED NOT BE RE-INITIALIZED:
	std::vector<double> dHSL;//membrand diffusion value
	std::vector<double> deltaHSL;//per-timestep production value
	std::vector<double> deltaPROTEIN;//per-timestep production value

	double	R_total;
    double	H_tau, L_tau, I_tau;


  protected:
	//ALL VARIABLES WILL WILL BE COPIED BY THE DEFAULT COPY CONSTRUCTOR:
	double  dt;//note: should not be changed by user due to delay buffer size
	double  nodesPerMicron;
	bool 	inductionFlag;
	double  pressureValue, pressureK50;

	//THE FOLLOWING WILL BE INITIALIZED ON DIVISION OF DAUGHTER CELL:
	std::queue<double> qH_tau;
	std::queue<double> qI_tau;
	std::queue<double> qL_tau;
	
	//THE FOLLOWING NEED NOT BE RE-INITIALIZED:
	double	ratio_H_tau, ratio_L_tau, ratio_I_tau;
	double 	totalProtein, totalHSL;
	double 	delta[numConcentrations];
	double  dC4HSL, dC14HSL;
	double  degradationProtein, degradationHSL;
    double nanoMolarPerMolecule;

    friend std::ostream &operator<<(std::ostream &os, const eQ::strainType &type);
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
    double			getProteinNumber(sendRecvStrain::qProteins which) {return iPROTEIN[which];}

//    sendRecvStrain(eQ::strainType strainType, const struct dso_parameters *p,
//                   const double timestep, const double nodesPerMicron)
//        : Strain(strainType, p, timestep, nodesPerMicron) {}

    //CONSTRUCTOR: (called for initial seed cells only)
    sendRecvStrain(eQ::strainType strainType, const double timestep, const double nodesPerMicron, const size_t numHSL)
        : Strain(strainType, nullptr, timestep, nodesPerMicron)
    {
        iHSL.assign(numHSL, 0.0);
        tHSL.assign(numHSL, 0.0);//per timestep value of popped queue value for convenience
        dHSL.assign(numHSL, 0.0);
        deltaHSL.assign(numHSL, 0.0);
            iPROTEIN.assign(qProteins::NUM_QUEUES, 0.0);
            tPROTEIN.assign(qProteins::NUM_QUEUES, 0.0);//per timestep value of popped queue value for convenience
            deltaPROTEIN.assign(qProteins::NUM_QUEUES, 0.0);

        std::queue<double> q;
        for(size_t j=0; j<queueDepth; j++) {q.push(0.0);}
        HSL_tau.assign(numHSL, q);
        PROTEIN_tau.assign(sendRecvStrain::qProteins::NUM_QUEUES, q);

    }
    //calls default copy constructor for derived class:
    std::shared_ptr<Strain>
        clone() const {return std::make_shared<sendRecvStrain>(*this);}  // requires C++ 14
    std::vector<double>
        computeProteins(const std::vector<double> &eHSL, const std::vector<double> &membraneD, const double lengthMicrons);

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
	//CONSTRUCTOR: (called for initial seed cells only)
	aspectRatioInvasionStrain(eQ::strainType strainType, const double timestep, const double nodesPerMicronScale, const size_t numHSL)
        : Strain(strainType, nullptr, timestep, nodesPerMicronScale)
    {
        iHSL.assign(numHSL, 0.0);
        tHSL.assign(numHSL, 0.0);//per timestep value of popped queue value for convenience
		dHSL.assign(numHSL, 0.0);
		deltaHSL.assign(numHSL, 0.0);

		iPROTEIN.assign(qProteins::NUM_QUEUES, 0.0);
		tPROTEIN.assign(qProteins::NUM_QUEUES, 0.0);//per timestep value of popped queue value for convenience
		deltaPROTEIN.assign(qProteins::NUM_QUEUES, 0.0);

		std::queue<double> q;
		for(size_t j=0; j<queueDepth; j++) {q.push(0.0);}

		HSL_tau.assign(numHSL, q);
		PROTEIN_tau.assign(qProteins::NUM_QUEUES, q);

    }

	//calls default copy constructor for derived class:
	std::shared_ptr<Strain>
		clone() const {return std::make_shared<aspectRatioInvasionStrain>(*this);}  // requires C++ 14

    std::vector<double>
        computeProteins(const std::vector<double> &eHSL, const std::vector<double> &membraneD, const double lengthMicrons);

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
    MODULUSmodule(const struct dso_parameters *p,
                  const double timestep, const double nodesPerMicron)
		: Strain(eQ::strainType::ACTIVATOR, p, timestep, nodesPerMicron)
	{
        //sample the K50 parameters for promoters
        //note:  will be reset for daughter cells
		setKH();

        //create a queue for time-averaging the cell length to reduce noise:
        double minsQ = double(eQ::parameters["MODULUS_TIME_AVERAGE_MINS"]);
        auto qdepth = size_t(ceil(minsQ/timestep));
        for(size_t i=0; i<qdepth; i++)
        {
            qReadout.push_front(0.0);
        }
        debugCounter=0;
	}

    //CONSTRUCTOR: (called for initial seed cells only)
    MODULUSmodule(const double timestep, const double nodesPerMicron, const size_t numHSL)
          : MODULUSmodule(nullptr, timestep, nodesPerMicron)
    {
        iHSL.assign(numHSL, 0.0);
        tHSL.assign(numHSL, 0.0);//per timestep value of popped queue value for convenience
        dHSL.assign(numHSL, 0.0);
        deltaHSL.assign(numHSL, 0.0);
            iPROTEIN.assign(qProteins::NUM_QUEUES, 0.0);
            tPROTEIN.assign(qProteins::NUM_QUEUES, 0.0);//per timestep value of popped queue value for convenience
            deltaPROTEIN.assign(qProteins::NUM_QUEUES, 0.0);

        std::queue<double> q;
        for(size_t j=0; j<queueDepth; j++) {q.push(0.0);}
        HSL_tau.assign(numHSL, q);
        PROTEIN_tau.assign(MODULUSmodule::qProteins::NUM_QUEUES, q);

    }
    //calls default copy constructor for derived class:
    std::shared_ptr<Strain>
        clone() const {return std::make_shared<MODULUSmodule>(*this);}  // requires C++ 14
    std::vector<double>
        computeProteins(const std::vector<double> &eHSL, const std::vector<double> &membraneD, const double lengthMicrons);

    void divideProteins(std::shared_ptr<MODULUSmodule> daughter, double fracToDaughter)
    {
		//copy parent values to daughter before re-sampling:
//        auto pmodulus = std::dynamic_pointer_cast<MODULUSmodule>(daughter);
		daughter->lnKL = lnKL;
		daughter->lnKH = lnKH;
		daughter->lnKH2 = lnKH2;
        reSampleK50s();
        daughter->reSampleK50s();

		Strain::divideProteins(daughter, fracToDaughter);

        //divide cell length averaging queue
        double rm = 1.0-fracToDaughter;
        for(size_t i=0;i<qReadout.size();i++)
        {//DIVIDE BY FRACTION
            auto data = qReadout.back();
                qReadout.push_front(rm*data);
                qReadout.pop_back();
					daughter->qReadout.push_front(fracToDaughter*data);
					daughter->qReadout.pop_back();
        }
        computeAverageCellLength();
	}
	void setKH()
	{//for initial seeds:
		//TODO: verify log normal generator was initialized
		lnKL =   eQ::getLogNormalRandomNumber(0);//PL_lac promoter K50
		lnKH =   eQ::getLogNormalRandomNumber(1);//Prhl or Pbad promoter K50
		lnKH2 =   eQ::getLogNormalRandomNumber(1);//Prhl or Pbad promoter K50
	}
    void reSampleK50s()
	{//for correlating values parent/daughter:
		//TODO: verify log normal generator was initialized
		auto correlateValue = [](double parent, double sample)
		{
            const double a = double(eQ::parameters["K50_correlationScale"]);//relative geometric weight of parent value
//            const double a = 0.9;//relative geometric weight of parent value
            return pow(parent, a)*pow(sample, 1-a);
		};
		lnKL =  correlateValue(lnKL, eQ::getLogNormalRandomNumber(0));//PL_lac promoter K50
		lnKH =  correlateValue(lnKH, eQ::getLogNormalRandomNumber(1));//Prhl or Pbad promoter K50
		lnKH2 =  correlateValue(lnKH2, eQ::getLogNormalRandomNumber(1));//Prhl or Pbad promoter K50
	}



    double lnKH, lnKL, lnKH2;
    double iptg;
    std::deque<double> qReadout;
    double filteredCellLength;

    int debugCounter;


    double computeAverageCellLength()
    {
        double sum = std::accumulate(qReadout.begin(), qReadout.end(), 0.0);
        filteredCellLength = sum / qReadout.size();
        return filteredCellLength;
    }
    void pushCellLengthValue(double value)
    {
        qReadout.push_front(value);
        qReadout.pop_back();
        computeAverageCellLength();
    }

};

void initRates(double *, struct rates &);
void loadParams(eQ::strainType which, struct dso_parameters &params, struct rates &rates);
	

#endif
