#ifndef EQ_H
#define EQ_H

#include <random>
#include <cmath>
#include <fstream>
#include <csignal>
#include <iostream>
#include <chrono>
#include <thread>
#include <pthread.h>
#include <memory>
#include <vector>
#include <list>
#include <forward_list>
#include <mpi.h>

#include <petscsys.h>
#include <petsc.h>
#include <petsclog.h>

//https://github.com/nlohmann/json
#include "../nlohmann/json.hpp"
// for convenience
using json = nlohmann::json;

class eQ
{
public:

    class tensorDataSource;

    //global cell parameters (fixed diameter here...may change cell width then put in cell model)
    //assume cell cylinder, radius 1/2 um;
    static constexpr double cylindricalCellVolumePerLength = M_PI*(0.5*0.5);//area of circle
    static constexpr double poleVolume = 4.0/3.0 * M_PI*(0.5*0.5*0.5);//both poles total

    static double computeCellVolumeActual(double cellLength)
    {//uses cylindrical body and hemispherical poles; models actual (ideal) cell volume
        return (cellLength - 1.0) * cylindricalCellVolumePerLength + poleVolume;
    }

    static double computeCellVolumeForConcentrations(double cellLength)
    {//09Oct.2019: using pole volumes leads to ~10% jump in conc. at division!
        return cellLength * cylindricalCellVolumePerLength;
    }

    static double proteinNumberToNanoMolar(double p, double cellLength)
    {//CONVERSION FROM PROTEIN # TO NANOMOLAR:
            // 1 nanoMolar [nM] = 1e-9mol/L    =    1e-9 * 0.602e24#/L  = 1e15*0.602 [protein#/liter]
            //1 um^3 = 1e-18m^3    =     1e-18m^3 * 1L/(0.1m)^3 = 1e-15L
            //==> 1 nM =  1e15*0.602 [protein#/L] * 1e-15 [L/um^3] = 0.602 [protein#/um^3]
            //==> to nanomolar: (p#/um^3)/0.602
//        return (p/computeCellVolumeActual(cellLength)) / 0.602;
        return (p/computeCellVolumeForConcentrations(cellLength)) / 0.602;
    }

    static double computeIntraCellularVolumeFraction(double cellLength)
    {
        //compute volume fraction of cell vs. rectangular region 1um x 1um x cellLength
//        return computeCellVolumeActual(cellLength) / (cellLength * 1.0 * 1.0);
        //use blunt-end cell volume to avoid division discontinuity:
        return computeCellVolumeForConcentrations(cellLength) / (cellLength * 1.0 * 1.0);
    }
    static double computeExtraCellularVolumeFraction(double cellLength)
    {
        //compute volume fraction of cell vs. rectangular region 1um x 1um x cellLength
        return 1.0 - computeIntraCellularVolumeFraction(cellLength);
    }
//    static double computeExtraCellularVolume(double cellLength)
//    {
//        //compute volume of outside of cell vs. rectangular region 1um x 1um x cellLength
//        return (cellLength * 1.0 * 1.0) -  computeCellVolumeActual(cellLength);
//    }
    static double computeVolumeRatio_ExtraToIntra(double cellLength)
    {
        return computeExtraCellularVolumeFraction(cellLength)/computeIntraCellularVolumeFraction(cellLength);
    }

//=======================================================================================
	typedef json			parametersType;
	typedef std::string		dataStringsType;

	static void initDefaultParameters()
	{
		eQ::initializedParameters = true;
        eQ::parameters["simulationTrapHeightMicrons"] = 20;
        eQ::parameters["simulationTrapWidthMicrons"] = 40;
        eQ::parameters["lengthScaling"] = 1.0;

        eQ::parameters["trapType"]      = "NOWALLED";
        eQ::parameters["boundaryType"]  = "DIRICHLET_0";
        eQ::parameters["cellInitType"]  = "RANDOM";
        eQ::parameters["simType"]		= "NO_SIGNALING";
        eQ::parameters["nodesPerMicronSignaling"] = 2;
        eQ::parameters["nodesPerMicronData"]       = 1;

        eQ::parameters["defaultAspectRatioFactor"]     = 1.0;
        eQ::parameters["mutantAspectRatioScale"]       = 1.0;
        eQ::parameters["aspectRatioThresholdHSL"]       = 1000.0;

        eQ::parameters["trapChannelLinearFlowRate"] = 100.0 * 60.0;//microns/sec * 60sec/min;
        eQ::parameters["channelSolverNumberIterations"] = 1;

        eQ::parameters["hslSignaling"]  = false;
        eQ::parameters["D_HSL"]			= {};

        eQ::parameters["numberSeedCells"] = 8;
        eQ::parameters["cellInitType"] = "RANDOM";

        eQ::parameters["openWalledDirichlet0"]     = true;//preset
        eQ::parameters["boundaries"] =
        {
            {"left", {"Dirichlet" , 0.0}}, //boundary value
//            {"left", {"Neumann" , 0.0}}, //normal derivative
//			{"left", {"Robin" , {0.0, 1.0, 0.0}}}, //alpha (Dirichlet), beta (normal derivative), gamma (boundary value)
            {"right", {"Dirichlet" , 0.0}},
            {"top", {"Dirichlet" , 0.0}},
            {"bottom", {"Dirichlet" , 0.0}}
        };
        //fraction = +/- 0.5*x
        eQ::parameters["divisionNoiseScale"] = 0.05;// = +/- 0.025
        eQ::parameters["MODULUS_TIME_AVERAGE_MINS"] = 1.0;
        eQ::parameters["K50_correlationScale"] = 0.0;

        eQ::parameters["channelLengthMicronsLeft"] = 100;
        eQ::parameters["channelLengthMicronsRight"] = 100;

        eQ::parameters["promoterDelayTimeMinutes"] = 8.0;

    }

    //declaration;  must be defined before main();
    static void populateFileStrings();

    static void init()
    {
        eQ::initDefaultParameters();
        eQ::populateFileStrings();
    }

	enum dataParameter
    {//NOTE: WHEN UPDATING, ADD THE FILENAMES to populateFileStrings() defined above main.
		//the data is written in  eQabm::updateCell()
		SPRING_COMPRESSION,
		CELL_ANGLE,
		CELL_VELOCITY,
		C4, C14,
		FP,
		CFP,
		YFP,
		SYN,
		LACI,
		AIIA,
		MFP,
		QTENSOR,
		GRADVEL_ISOTROPIC,
        LEGI_DATA_H,
        LEGI_DATA_A,
        C4RHL,
        RHL_T,
        MODULUS_S,
        MODULUS_H,
        MODULUS_R,

		DTENSOR_11,
		DTENSOR_22,
		DTENSOR_12,

		NUM_DATAPARAMETERS
	};


	typedef std::vector<
				std::pair<std::shared_ptr<eQ::tensorDataSource>, eQ::dataParameter>>
			dataFiles_t;

	static dataFiles_t initDataRecording(std::vector<eQ::dataParameter> &dataToRecord)
	{
		size_t n = parameters["nodesPerMicronData"];
		size_t y = parameters["simulationTrapHeightMicrons"];
		size_t x = parameters["simulationTrapWidthMicrons"];

		dataFiles_t dataRecords;
		size_t		rank;

		for(auto data : dataToRecord)
		{
			if(eQ::dataParameter::CELL_VELOCITY == data)
				rank=1;
			else
				rank=0;

			dataRecords.push_back(std::make_pair(
								  std::make_shared<tensorDataSource>(n,y,x,rank),
								  data));
		}
		return dataRecords;
	}



    class diffusionSolver
    {
    public:
        struct params
        {
            size_t                      uniqueID;
            MPI_Comm                    comm;
            double                      dt;
            double                      D_HSL;
            std::string                 filePath;
            eQ::dataFiles_t             dataFiles;
            double                      trapHeightMicrons;
            double                      trapWidthMicrons;
            double                      nodesPerMicron;
            double                      trapChannelVelocity;
        };
//        virtual void initDiffusion(MPI_Comm comm, std::vector<std::string> filePaths, int argc, char* argv[]) =0;// {}
//        virtual void initDiffusion(size_t id, MPI_Comm comm, std::string filePath, double D, double dt, int argc, char* argv[]) =0;// {}
        virtual void initDiffusion(eQ::diffusionSolver::params &) =0;// {}
        virtual void stepDiffusion() {}
        virtual void setBoundaryValues(const eQ::parametersType &bvals) =0;
        virtual eQ::parametersType getBoundaryFlux(void) =0; //must at least define value of "totalFlux"
        virtual void writeDiffusionFiles(double timestamp) =0;
        virtual void finalize(void) {}
    };

//	enum trapType
//	{
//		NOTRAP,
//		NOWALLED,
//		TWOWALLED,
//		THREEWALLED,
//		ONEWALLED_TOP,
//		ONEWALLED_LEFT,
//		LEFTCORNER_HALF,
//		OPPOSITE_CORNERS,
//		NUM_TRAPTYPES
//	};
//    enum boundaryType
//    {
//        DIRICHLET_0,
//        NEUMANN_0,
//        ROBIN,
//        NEUMANN_3WALLED_TEST,
//        NUM_BOUNDARYTYPES
//    };
//    enum cellInitType
//    {
//        RANDOM,
//        AB_HALF,
//        ABA_THIRDS,
//        NUM_CELLINITTYPES
//    };
//    enum simType
//	{
//		NO_SIGNALING,
//		SENDER_RECEIVER,
//        TOGGLE_SWITCH,
//        CONSENSUS_TOGGLE,
//		DUALSTRAIN_OSCILLATOR,
//        STATIC_ASPECTRATIO,
//        DSO_LATTICE,
//        SINGLEMUTANT_ASPECTRATIO,
//        INDUCED_DYNAMIC_ASPECTRATIO,
//		NUM_SIMTYPES
//	};
    enum HSLType
	{
		C4HSL,
		C14HSL,
		NUM_HSLTYPES
	};
	enum class strainType
	{
	    ACTIVATOR, REPRESSOR,
	    X,Y,Z,
		NUM_STRAINTYPES
	};


//=======================================================================================
	static inline std::pair<double, double> xy_from_ij(size_t i, size_t j, double n)
	{//Note: j is horizontal (x) direction, i is vertical (y)
		return std::make_pair(double(j)/n,double(i)/n);
	}
	static inline std::pair<size_t, size_t> ij_from_xy(double x, double y, double n)
	{//Note: j is horizontal (x) direction, i is vertical (y)
		size_t j = size_t(round(x * n));
		size_t i = size_t(round(y * n));
		return std::make_pair(i,j);
	}
	static inline std::pair<size_t, size_t> ij_from_xy(std::pair<double, double> xy, double n)
    {//Note: j is horizontal (x) direction, i is vertical (y)
        size_t j = size_t(round(xy.first * n));
        size_t i = size_t(round(xy.second * n));
        return std::make_pair(i,j);
    }
	static inline size_t index_from_ij(size_t i, size_t j, size_t nh, size_t nw)
    {//Note: j is horizontal (x) direction, i is vertical (y)
        return ((i < nh) && (j < nw)) ? j+i*nw : 0;
    }
//=======================================================================================
	template<class T>
    class gridFunction
    {
    public:
        gridFunction(size_t nh, size_t nw)
                    : nodesHigh(nh), nodesWide(nw)
        {
            for(size_t i(0); i<nodesHigh; i++)
                grid.push_back(std::vector<T>(nodesWide, T(0)));
        }
		void clear()
		{
			if(grid.empty())
			{
				std::cout<<"grid empty error!"<<std::endl;
				return;
			}
			for(size_t i(0); i<nodesHigh; i++)
				for(size_t j(0); j<nodesWide; j++)
					grid[i][j] = T(0);
		}
		void assign(T value)
		{
			if(grid.empty())
			{
				std::cout<<"grid empty error!"<<std::endl;
				return;
			}
			for(size_t i(0); i<nodesHigh; i++)
				for(size_t j(0); j<nodesWide; j++)
					grid[i][j] = T(value);
		}
		bool isValidIndex(std::pair<size_t, size_t> point)
        {
            return (point.first < nodesHigh) && (point.second < nodesWide);
        }
        std::vector<std::vector<T>> grid;
    private:
        size_t nodesHigh, nodesWide;
    };
//=======================================================================================
	template<class T>
	static void linearize2Dgrid(std::shared_ptr<gridFunction<T>> g, std::vector<T> &dataArray)
	{
		auto nrows = g->grid.size();
		auto ncols = g->grid[0].size();
		size_t index=0;
		for(size_t i(0); i<nrows; i++)
			for(size_t j(0); j<ncols; j++)
				dataArray.at(index++) = g->grid[i][j];
	}

//=======================================================================================
    class tensorDataSource
    {
	public:
		tensorDataSource(size_t n, size_t y, size_t x, size_t thisRank)
			: nodesPerMicron(n),//use the Data version of per-micron
			  rank(thisRank)//scalar data source default
		{
			auto nh = y*n + 1;
			auto nw = x*n + 1;
			if(0 == rank)
			{
				dataGrid = std::make_shared<eQ::gridFunction<double>>(nh,nw);
			}
			else
			{
				xdataGrid = std::make_shared<eQ::gridFunction<double>>(nh,nw);
				ydataGrid = std::make_shared<eQ::gridFunction<double>>(nh,nw);
			}
		}

        tensorDataSource(const tensorDataSource &)
            {std::cout<<"tensorDataSource copyConstructor not defined!\n";}

        size_t getRank(){return rank;}
        double eval(double x, double y)
        {
//			auto ij = getIndices(x,y);
			auto ij = ij_from_xy(x,y,nodesPerMicron);
			return dataGrid->grid[ij.first][ij.second];
        }
        std::pair<double,double> evalVector(double x, double y)
        {
//			auto ij = getIndices(x,y);
			auto ij = ij_from_xy(x,y,nodesPerMicron);
			return std::make_pair(
                xdataGrid->grid[ij.first][ij.second],
                ydataGrid->grid[ij.first][ij.second]);
        }
        //todo: make these a vector of pointers to grid, auto-sized for tensor rank
        std::shared_ptr<eQ::gridFunction<double>> dataGrid;
        std::shared_ptr<eQ::gridFunction<double>> xdataGrid;
        std::shared_ptr<eQ::gridFunction<double>> ydataGrid;

    private:
//        std::pair<size_t, size_t> getIndices(double x, double y)
//        {
//            size_t j = size_t(round(x * nodesPerMicron));
//            size_t i = size_t(round(y * nodesPerMicron));
//            return std::make_pair(i,j);
//        }
        double nodesPerMicron;
        size_t rank;
    };
//=======================================================================================
//=======================================================================================
//=======================================================================================
//=======================================================================================
	class uniformRandomNumber
    {
    public:
        uniformRandomNumber()
        {
    //      use a random device generator to generate seed:
//            std::random_device rdev{};
//            generator.seed(rdev());
//            distribution =  std::make_shared<std::uniform_real_distribution<double>>(0.0,1.0);
            // construct a trivial random generator engine from a time-based seed:
            auto seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::cout<<"Called uniformRandomNumber() which generated seed="<<seed<<std::endl;
            uniformRandomNumber(size_t(seed));
        }
        uniformRandomNumber(size_t seed)
        {
            generator.seed(unsigned(seed));
//            generator.seed(0);
            distribution =  std::make_shared<std::uniform_real_distribution<double>>(0.0,1.0);
            std::cout<<"Initializing uniformRandomNumber(size_t seed) with seed="<<seed<<std::endl;
        }
        double randomNumber()
        {
            return (*distribution)(generator);
        }
    private:
        std::default_random_engine generator;
        std::shared_ptr<std::uniform_real_distribution<double>> distribution;
    };
//=======================================================================================

    //static class members which must be defined outside of main
    static eQ::dataStringsType dataStrings[];
    static bool initializedParameters;
    static eQ::parametersType parameters;
    static std::default_random_engine generator;
    static std::vector<std::shared_ptr<std::lognormal_distribution<double>>>
        lndistributions;

	static std::vector<double> solutionVector;
    static eQ::uniformRandomNumber zeroOne;

    static void initLogNormalGenerators(std::vector<std::pair<double, double>> &meansStdevs)
    {
//      use a random device generator to generate seed:
//        std::random_device rdev{};
//        generator.seed(rdev());
        // construct a trivial random generator engine from a time-based seed:
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        generator.seed(static_cast<unsigned long>(seed));
        for(auto &init : meansStdevs)
        {
            //convert to (lognormal) mean and variance of the log of target distribution
            double lnmu = log(init.first/sqrt(1.0 + init.second/(init.first*init.first)));
            double lns2 = log(1.0 + init.second/(init.first*init.first));

            lndistributions.push_back(std::make_shared<std::lognormal_distribution<double>>(lnmu, lns2));
        }
    }
//    static double getLogNormalMean(size_t which)
//    {
//        return (lndistributions.size() > which) ? lndistributions[which]->m(): 0.0;
//    }
    static double getLogNormalRandomNumber(size_t which)
    {
        return (lndistributions.size() > which) ? (*(lndistributions[which]))(generator) : 0.0;
    }
//    static double getLogNormalRandomNumber(size_t which, double mean, double var)
//    {
//        if(lndistributions.size() > which)
//        {
//            //update the mean, var
//            lndistributions[which]->param(std::lognormal_distribution<double>::param_type(mean, var));
//            return (*(lndistributions[which]))(generator);
//        }
//        return mean;
//    }

    friend std::ostream &operator<<(std::ostream &os, const eQ::strainType &type);

    class databaseClass
    {
    public:
        std::vector<eQ::dataStringsType> dataStrings;
        std::default_random_engine generator;
        std::vector<std::shared_ptr<std::lognormal_distribution<double>>>
            lndistributions;
        std::vector<std::shared_ptr<std::uniform_real_distribution<double>>>
            uniformRandomNumbers;
        bool initializedParameters;
        eQ::parametersType parameters;
        void init_default_random_engine()
        {
            // construct a trivial random generator engine from a time-based seed:
            auto seed = std::chrono::system_clock::now().time_since_epoch().count();
            generator.seed(static_cast<unsigned long>(seed));
        }
        void initLogNormalGenerators(std::vector<std::pair<double, double>> &meansStdevs)
        {
            for(auto &init : meansStdevs)
            {
                lndistributions.push_back(std::make_shared<std::lognormal_distribution<double>>(init.first, init.second));
            }
        }
        double getLogNormalRandomNumber(size_t which)
        {
            return (lndistributions.size() > which) ? (*(lndistributions[which]))(generator) : 0.0;
        }
        void initUniformRandomNumberGenerators(std::vector<std::pair<double, double>> &lowHighs)
        {
            for(auto &init : lowHighs)
            {
                uniformRandomNumbers.push_back(std::make_shared<std::uniform_real_distribution<double>>(init.first, init.second));
            }
        }
        double getUniformRandomNumber(size_t which)
        {
            return (uniformRandomNumbers.size() > which) ? (*(uniformRandomNumbers[which]))(generator) : 0.0;
        }
    };

    static eQ::databaseClass database;
};


#endif // EQ_H
