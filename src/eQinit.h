#ifndef EQINIT_H
#define EQINIT_H

#include "eQ.h"

namespace eQ {

void initDiffusionRates()
{
	eQ::data::physicalDiffusionRates =
	{
		{"C4", std::vector<double> {3.0e4, 3.0}},
		{"C14", std::vector<double>{1.6e4, 2.1}}
	};
}

void initDefaultParameters()
	{
		eQ::data::initializedParameters = true;
		eQ::data::parameters["simulationTrapHeightMicrons"] = 20;
		eQ::data::parameters["simulationTrapWidthMicrons"] = 40;
		eQ::data::parameters["lengthScaling"] = 1.0;

		eQ::data::parameters["trapType"]      = "NOWALLED";
		eQ::data::parameters["boundaryType"]  = "DIRICHLET_0";
		eQ::data::parameters["cellInitType"]  = "RANDOM";
		eQ::data::parameters["simType"]		= "NO_SIGNALING";
		eQ::data::parameters["nodesPerMicronSignaling"] = 2;
		eQ::data::parameters["nodesPerMicronData"]       = 1;

		eQ::data::parameters["defaultAspectRatioFactor"]     = 1.0;
		eQ::data::parameters["mutantAspectRatioScale"]       = 1.0;
		eQ::data::parameters["aspectRatioThresholdHSL"]       = 1000.0;

		eQ::data::parameters["trapChannelLinearFlowRate"] = 100.0 * 60.0;//microns/sec * 60sec/min;
		eQ::data::parameters["channelSolverNumberIterations"] = 1;

		eQ::data::parameters["hslSignaling"]  = false;
		eQ::data::parameters["D_HSL"]			= {};

		eQ::data::parameters["numberSeedCells"] = 8;
		eQ::data::parameters["cellInitType"] = "RANDOM";

		eQ::data::parameters["PETSC_SIMULATION"] = false;
		eQ::data::parameters["boundaries"] =
		{
			{"left", {"Dirichlet" , 0.0}}, //boundary value
//            {"left", {"Neumann" , 0.0}}, //normal derivative
//			{"left", {"Robin" , {0.0, 1.0, 0.0}}}, //alpha (Dirichlet), beta (normal derivative), gamma (boundary value)
			{"right", {"Dirichlet" , 0.0}},
			{"top", {"Dirichlet" , 0.0}},
			{"bottom", {"Dirichlet" , 0.0}}
		};
		//fraction = +/- 0.5*x
		eQ::data::parameters["divisionNoiseScale"] = 0.05;// = +/- 0.025
		eQ::data::parameters["MODULUS_TIME_AVERAGE_MINS"] = 1.0;
		eQ::data::parameters["K50_correlationScale"] = 0.0;

		eQ::data::parameters["mediaChannelMicronsLeft"] = 100;
		eQ::data::parameters["mediaChannelMicronsRight"] = 100;
		eQ::data::parameters["simulationChannelLengthLeft"] = 100;
		eQ::data::parameters["simulationChannelLengthRight"] = 100;

		eQ::data::parameters["AnisotropicDiffusion_Axial"] = 1.0;
		eQ::data::parameters["AnisotropicDiffusion_Transverse"] = 1.0;

	}
}

#endif // EQINIT_H
