#ifndef EQCELL_H
#define EQCELL_H

#include <cmath>
#include <iostream>

namespace eQ {
class Cell
{
public:
	//virtual destructor needed in base class so that derived classes can be deleted properly:
	virtual ~Cell() = default;

	enum class strainType
	{
		ACTIVATOR, REPRESSOR,
		X,Y,Z,
		NUM_STRAINTYPES
	};
	struct Params
	{
		double mass, moment, length, width;
		double x, y, angle, vx, vy, av;
		double divisionLength;
		double meanDivisionLength;
		double doublingPeriodMinutes;
	};
	virtual void       setCellID(long i) {ID = i;}
	virtual long       getCellID() {return ID;}
	virtual void       setParentID(long i) {parentID = i;}
	virtual long       getParentID() {return parentID;}
	virtual double     getLengthMicrons()=0;//length in microns
	virtual double     getWidthMicrons()=0;//length in microns
	virtual double     getCenter_x()=0;
	virtual double     getCenter_y()=0;
	virtual double     getAngle()=0;

    virtual double     getDivisionLength()                  { return baseData.divisionLength; }
    virtual void       setDivisionLength(double dl)         { baseData.divisionLength = dl; }
    virtual double     getMeanDivisionLength()              { return baseData.meanDivisionLength; }
    virtual void       setMeanDivisionLength(double mdl)    { baseData.meanDivisionLength = mdl; }


	static constexpr double  DEFAULT_CELL_MASS                      =             1.0;
	static constexpr double  DEFAULT_CELL_MOMENT                    =             100.0;
	static constexpr double  DEFAULT_CELL_WIDTH_MICRONS             =             1.0;
	static constexpr double  DEFAULT_DIVISION_LENGTH_MICRONS        =             4.2;
	static constexpr double  DEFAULT_CELL_DOUBLING_PERIOD_MINUTES   =             20.0;
	static constexpr double  DEFAULT_PROMOTER_DELAY_TIME_MINUTES	=             7.5;//see Chen et al. Science 2015


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

protected:
	virtual double  computeDivisionLength(double lengthAtBirth, double meanLengthAtBirth, double alpha)
	{
		//15Feb.2019:  cell size regulation using ~"incremental model" by Ariel Amir (PRL 112, 2014)
		//"sizer" regulation <== alpha=1
		return 2.0 * pow(lengthAtBirth, 1.0-alpha) * pow(meanLengthAtBirth, alpha);
	}
	virtual void    computeDivisionLength(double birthLength)
	{
		double alpha = 0.5;
		baseData.divisionLength = computeDivisionLength(birthLength, 0.5*baseData.meanDivisionLength, alpha);
	};

	Cell::Params     baseData;
	double          _dt;
	long            ID;
	long            parentID;

	friend std::ostream &operator<<(std::ostream &os, const eQ::Cell::strainType &type);

};

}
#endif // EQCELL_H
