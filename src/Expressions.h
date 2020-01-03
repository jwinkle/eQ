#ifndef EXPRESSIONS_H
#define EXPRESSIONS_H

#include <dolfin.h>
#include "eQ.h"
#include "./abm/eColi.h"


using namespace dolfin;
#include <mshr.h>
using namespace mshr;

class updatingDirchletBoundary: public Expression
{
  public:
    updatingDirchletBoundary(Function & funcRef)
    {
        pFunction = &funcRef;
    }

    updatingDirchletBoundary(double bval) : boundaryValue(bval), pFunction(nullptr) {}

	void eval(Array<double>& values, const Array<double>& x) const
	{
        Array<double> myval(1);

        if(nullptr == pFunction)
        {
            values[0] = boundaryValue;
        }
        else
        {
            pFunction->eval(values, x);
        }

	}
	void updateBoundary(double bval)
	{
		boundaryValue = bval;
	}

	double boundaryValue;
    Function *pFunction;
};

class scalarDataExpression : public Expression
{
public:
//    scalarDataExpression(eQ::tensorDataSource *source)
	scalarDataExpression(std::shared_ptr<eQ::tensorDataSource> source)
		: dataSource(source) {}
	void eval(Array<double>& values, const Array<double>& x) const
	{
		values[0] = dataSource->eval(x[0], x[1]);
	}
private:
	std::shared_ptr<eQ::tensorDataSource> dataSource;
};
class vectorDataExpression : public Expression
{
public:
//    vectorDataExpression(eQ::tensorDataSource *source)
	vectorDataExpression(std::shared_ptr<eQ::tensorDataSource> source)
		: Expression(2), dataSource(source) {}
	void eval(Array<double>& values, const Array<double>& x) const
	{
		auto vec = dataSource->evalVector(x[0], x[1]);
		values[0] = vec.first;
		values[1] = vec.second;
	}
private:
	std::shared_ptr<eQ::tensorDataSource> dataSource;
};
class zeroInitScalar : public Expression
{
  public:
  void eval(Array<double>& values, const Array<double>& x) const
  {
	  values[0]=0.0;
  }

};
class zeroInitVector : public Expression
{
  public:
  zeroInitVector():Expression(2) {}
  void eval(Array<double>& values, const Array<double>& x) const
  {
	  values[0]=0.0;
	  values[1]=0.0;
  }

};

class AnisotropicDiffusionTensor: public Expression
{
  public:
	AnisotropicDiffusionTensor(std::shared_ptr<std::vector<double>> D,
			size_t nodeDensity, size_t nodesHigh, size_t nodesWide, double test)
		: D(D), n(nodeDensity), nh(nodesHigh), nw(nodesWide), testValue(test)
	{
		maxNodes = nh*nw;
//        std::cout<<"size of D: "<<D->size()<<std::endl;
//        std::cout<<"size of nodeDensity: "<<nodeDensity<<std::endl;
//        std::cout<<"size of nodesHigh: "<<nodesHigh<<std::endl;
//        std::cout<<"size of nodesWide: "<<nodesWide<<std::endl;
	}
	void printSize(void)
	{
		std::cout<<"size of D: "<<D->size()<<std::endl;
	}
	void eval(Array<double>& values, const Array<double>& x) const
	{
		auto index = eQ::ij_from_xy(x[0], x[1], n);
		auto entry = index.first*nw + index.second;
		values[0] = (entry < maxNodes) ? D->at(entry) : testValue;
//        values[0] = testValue;
	}
	std::shared_ptr<std::vector<double>> D;
	size_t n, nh, nw, maxNodes;
	double testValue;
private:
	bool sizeFlag;
};


class ThetaClass : public Expression
{
  public:
  void eval(Array<double>& values, const Array<double>& x) const
  {
	double theta;
	// theta = 2.0*DOLFIN_PI*(x[0]/gridWidth);
	// values[0] = (x[0] < gridWidth/2.0) ? 0.0: DOLFIN_PI/2.0;
	// theta = DOLFIN_PI/4.0;
	// theta = 0.0;
	theta = 1.0;
	// double plusminus_point5 = 0.5 - dolfin::rand();
	// theta =  2.0*DOLFIN_PI*dolfin::rand();

	// theta =  1.0*DOLFIN_PI*dolfin::rand();
	// values[0] = theta;

	//testing:
	auto dx = x[0];
	auto dy = x[1];
	auto r = sqrt( pow(0.25-dy, 2.0) + pow(1.5-dx, 2.0));

	//hijack theta class for eta:
	// values[0] = (r < 0.1) ? 1.0 : 0.0;
	values[0] = (r < 0.1) ? DOLFIN_PI/4.0 : 0.0;
	// values[0] = (r < 0.1) ? DOLFIN_PI + DOLFIN_PI/4.0 : 0.0;
	// values[0] = 0.0;

	// values[0] = DOLFIN_PI/4.0;
	// values[0] = DOLFIN_PI + DOLFIN_PI/4.0;
	// values[0] = -DOLFIN_PI/2.0;

	// //BAND-LIMITED INIT:
	//   auto scale = 0.125;
	//   auto a1=DOLFIN_PI*5.0;
	//   auto a2=30.0;
	//   auto a3=a2+DOLFIN_PI/3.0;
	//   auto a4=a2-DOLFIN_PI/3.0;
	//   values[0] = 0.5 + scale*(
	//     sin(a1*DOLFIN_PI*x[0])*sin(a2*DOLFIN_PI*x[0])
	//     *sin(a3*DOLFIN_PI*x[1])*sin(a4*DOLFIN_PI*x[1])
	//     );
  }
};
class GammaClass : public Expression
{
  public:
  void eval(Array<double>& values, const Array<double>& x) const
  {
//    OcclusionBoundary boundary;
//    bool onBC = boundary.inside(x);
	// values[0] = onBC ? 10.0:1.0;
	// values[0] = onBC ? (1.0 + 2.0*boundary.scaleValue):1.0;
	//occlusion OFF:
	values[0] = 1.0;
  }
  private:
};
class RhoClass : public Expression
{
  public:
  void eval(Array<double>& values, const Array<double>& x) const
  {
	double a = 0.2;
	double dx = x[0]-0.4;
	double dy = x[1]-0.4;
	double r = sqrt( pow(dx,2.0) + pow(dy,2.0));
	// values[0] = r<a ? 1.0:0.0;
	values[0] = 1.0;
  }
};
class AlphaGammaRhoClass : public Expression
{
 public:
  AlphaGammaRhoClass(GammaClass &fg, RhoClass &fr)//sets aplpha to 1.0
  {
	falpha = nullptr; constantAlpha = true; alpha = 1.0;
	fgamma = &fg;
	frho = &fr;
  }
  AlphaGammaRhoClass(Function &fa, GammaClass &fg, RhoClass &fr)
  {
	constantAlpha = false;
	falpha = &fa;
	fgamma = &fg;
	frho = &fr;
  }
  AlphaGammaRhoClass(double a, GammaClass &fg, RhoClass &fr)
  {
	falpha = nullptr; constantAlpha = true;
	alpha = a;
	fgamma = &fg;
	frho = &fr;
  }
  void setGenerator(std::shared_ptr<std::default_random_engine> gen)
  {
	  generator = gen;
	  noiseFlag = true;
  }
 void eval(Array<double>& values, const Array<double>& x) const
  {
	//NOTE: output value mus be positive
	double agr;
	Array<double> myval(1);
	if(constantAlpha)
	{
	  agr = alpha;
	}
	else
	{
	  falpha->eval(myval, x);
	  agr = myval[0];
	}
	fgamma->eval(myval, x);
	agr *= myval[0];
	  frho->eval(myval, x);
	  agr *= myval[0];

	values[0] = agr;
 }
  private:
  Function *falpha;
  GammaClass *fgamma;
  RhoClass *frho;
  std::shared_ptr<std::default_random_engine>  generator;
  bool noiseFlag = false;

  bool constantAlpha = false;
  double alpha = 0.0;
};
class QScalar : public Expression
  {
  public:
	QScalar(Function &fref1, Function &fref2, int which,
			  std::shared_ptr<std::default_random_engine> gen, double scale)
	{
	  f1 = &fref1;//Qa scalar
	  f2 = &fref2;//Qb scalar
	  _which = which;
	  generator = gen;
	  noiseScale = scale;
	  // QScalar(fref1,fref2,which);
	}
	void eval(Array<double>& values, const Array<double>& x) const
	{
	  double angle1,angle2;
	  Array<double> myval(1);

	  f1->eval(myval, x);
	  auto ua = myval[0];//ua = cos(2theta)
		f2->eval(myval, x);
		auto ub = myval[0];//ub = sin(2theta)

	  auto qnorm = sqrt(ua*ua + ub*ub);//norm = 1/2 q-scalar = eigenvalue

	  // if(qnorm == 0.0)
	  if(qnorm < DOLFIN_EPS)
	  {
		 angle1 = 0.;
		 angle2 = 0.;
	  }
	  else
	  {
		 angle1 = acos(ua/qnorm);//acos() in [0,pi]
		 angle2 = asin(ub/qnorm);//asin() in [-pi/2, +pi/2]
	  }

		// if(ub>0.0)
		if(ub<0.0)// y<0, quadrants III, IV
		{//acos conversion
		  angle1 = 2.0*DOLFIN_PI - angle1;
		}
		// if(ua>0.0)
		  // angle2 = -angle2;//y<0, quadrants III, IV
		// else
		//   angle2 = DOLFIN_PI + angle2;
		if(ua<0.0) //x<0, quadrants II, III
		{
		  angle2 = DOLFIN_PI - angle2;
		}
		if(angle2 < 0.0) angle2 += 2.0*DOLFIN_PI;

		angle1 *= 0.5;
		angle2 *= 0.5;

	  auto averageAngle = 0.5*(angle1+angle2);

	  //these values are for writing to file (and to throttle via 1-q^2 term)
	  if(0 == _which)
	  {
		values[0] = 2.0*qnorm;//returns the q-scalar order parameter
		  // values[0] = ua;
		  // values[0] = sqrt(ua*ua);
	  }
	  else if(1 == _which)
		values[0] = averageAngle;//returns Q director
		  // values[0] = ub;
		  // values[0] = angle1;
		  // values[0] = angle2;
	  else
	  {
		// const double scale = 64.0;//pi/32 = ~ 6deg
		double noisevar;
		noisevar = (noiseScale < DOLFIN_EPS) ? 0.0:DOLFIN_PI/noiseScale;
		std::normal_distribution<double> distribution(0.0,noisevar);
		double rn = distribution(*generator);
		auto an = 2.0*(averageAngle + rn);
		// auto an = averageAngle;
		if(2 == _which)
		{
		  values[0] = cos(an)*qnorm;
		  // std::cout<<ua<<", "<<values[0]<<std::endl;
		  // values[0] = ua;
		}
		else if(3 == _which)
		{
		  // if(ua>0.0) an = -an;
		  values[0] = sin(an)*qnorm;
		  // std::cout<<ub<<", "<<values[0]<<std::endl;
		  // values[0] = ub;
		}
	  }
	}
  private:
	double noiseScale;
	Function *f1, *f2;
	int _which;
	std::shared_ptr<std::default_random_engine>  generator;
	int debugCounter=0;
  };
class StrainRate : public Expression
{//forms the symmetric, traceless strain rate, accessed by component, from grad(v) tensor
	//used soley for Q-dot.
	public:
	StrainRate(Function &fref, size_t w)
	{
	  f = &fref;//tensor strain rate
	  which = w;//only if tensor, which entry to track: 1,1 or 1,2
	}
	void eval(Array<double>& values, const Array<double>& x) const
	{
	  Array<double> myval(4);//2x2 tensor
	  f->eval(myval, x);

	  double k11 = myval[0];
	  double k12 = myval[1];
	  double norm = sqrt(k11*k11 + k12*k12);
	  //compute trace/2
	  double traceBy2 = 0.5*(myval[0]+myval[3]);

	  //kinetic rate parameter
	  // double B = 1.0;
	  double B = 0.5;

		if(0 == which)
		  values[0] = B*(k11-traceBy2);
		  // values[0] = myval[1];
		else if(1 == which)
		  values[0] = B*k12;
		  // values[0] = myval[2];
	}
	private:
	Function *f;
	size_t which;
};

class LEGI_boundaryValues : public Expression
{
  public:
	LEGI_boundaryValues(double gridHeight, double gridWidth, double c0val, double gval)
		: trapH(gridHeight), trapW(gridWidth), c0(c0val), g(gval) {}
	void eval(Array<double>& values, const Array<double>& x) const
	{
	  Array<double> myval(1);
	  values[0] = c0 + g*(x[0] - trapW/2.0);
	}
	double trapH, trapW, c0, g;
};

#endif // EXPRESSIONS_H
