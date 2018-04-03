#ifndef P2H2D_IFUNCTIONAL_H
#define P2H2D_IFUNCTIONAL_H

#include <function.h>
#include <gradient.h>
#include <projection.h>
#include <printer.h>
#include "iproblem2h2d.h"

using namespace IProblem2H;

namespace IProblem2H
{

class IFunctional : public RnFunction, public IGradient, public IProjection, public IPrinter
{
public:
    IFunctional();

    virtual double fx(const DoubleVector &prms) const;
    virtual double integral(const DoubleMatrix &u, const DoubleMatrix &ut) const;
    virtual double integral1(const DoubleMatrix &u, const DoubleMatrix &ut) const;
    virtual double integral2(const DoubleMatrix &u, const DoubleMatrix &ut) const;
    virtual double norm(const EquationParameter& eprm, const OptimizeParameter &oprm, const OptimizeParameter &oprm0) const;
    virtual double penalty(const vector<ExtendedSpaceNode2DH> &info, const OptimizeParameter &o_prm) const;
    virtual double gpi(unsigned int i, unsigned int layer, const vector<ExtendedSpaceNode2DH> &info, const OptimizeParameter &o_prm) const;
    virtual double g0i(unsigned int i, unsigned int layer, const vector<ExtendedSpaceNode2DH> &info, const OptimizeParameter &o_prm) const;
    //virtual double mu(double x, double y) const;

    virtual void gradient(const DoubleVector &prms, DoubleVector &g) const;

    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, GradientMethod::MethodResult result) const;
    virtual void project(DoubleVector &x, unsigned int index);

    //void setGridParameters(Dimension timeDimension, Dimension spaceDimensionX, Dimension spaceDimensionY);
    //void setGridTimeDimension(Dimension timeDimension);
    //void setEquationParameters(double a, double lambda0, double lambda);
    //void setIntTemperature(double fi);
    //void setEnvTemperature(double theta);
    //void setRegEpsilon(double regEpsilon);
    //void setPenaltyCoefficient(double r);
    //void setPenaltyLimits(const DoubleVector &vmin, const DoubleVector &vmax);
    //void setParameter(const Parameter &parameter);
    //void setParameter0(const Parameter &parameter0);

    void toVector(const OptimizeParameter &prm, DoubleVector &x) const;
    void fromVector(const DoubleVector &x, OptimizeParameter &prm) const;

public:
    OptimizeParameter mOptParameter0;
    EquationParameter mEquParameter;

    Dimension mTimeDimension;
    Dimension mSpaceDimensionX;
    Dimension mSpaceDimensionY;

    DoubleVector vmin;
    DoubleVector vmax;

    double regEpsilon;
    //double r;
public:
    DoubleMatrix V0;
    double alpha0;
    DoubleMatrix V1;
    double alpha1;

    bool optimizeK;
    bool optimizeZ;
    bool optimizeC;
    bool optimizeO;

    double r;
};

}

#endif // P2H2D_IFUNCTIONAL_H
