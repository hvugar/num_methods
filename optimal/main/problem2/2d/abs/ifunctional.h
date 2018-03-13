#ifndef IFUNCTIONAL_H
#define IFUNCTIONAL_H

#include <function.h>
#include <gradient.h>
#include <projection.h>
#include <printer.h>
#include "problem2forward2d.h"
#include "problem2backward2d.h"

class IFunctional : public RnFunction, public IGradient, public IProjection, public IPrinter
{
public:
    IFunctional();

    virtual double fx(const DoubleVector &prms) const;
    virtual double integral(const DoubleMatrix &u) const;
    virtual double norm() const;
    virtual double penalty(vector<ExtendedSpaceNode2D> &info) const;
    virtual double gpi(unsigned int i, unsigned int layer, const vector<ExtendedSpaceNode2D> &info) const;
    virtual double g0i(unsigned int i, unsigned int layer, const vector<ExtendedSpaceNode2D> &info) const;
    virtual double mu(double x, double y) const;

    virtual void gradient(const DoubleVector &prms, DoubleVector &g) const;

    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, GradientMethod::MethodResult result) const;
    virtual void project(DoubleVector &x, unsigned int index);

    void setGridParameters(Dimension timeDimension, Dimension spaceDimensionX, Dimension spaceDimensionY);
    void setGridTimeDimension(Dimension timeDimension);
    void setEquationParameters(double a, double lambda0, double lambda);
    void setIntTemperature(double fi);
    void setEnvTemperature(double theta);
    void setRegEpsilon(double regEpsilon);
    void setPenaltyCoefficient(double r);
    void setPenaltyLimits(const DoubleVector &vmin, const DoubleVector &vmax);
    void setParameter(const Parameter &parameter);
    void setParameter0(const Parameter &parameter0);

    void toVector(const Parameter &prm, DoubleVector &x) const;
    void fromVector(const DoubleVector &x, Parameter &prm);

private:
    double sgn(double x) const;

protected:
    Parameter mParameter;
    Parameter mParameter0;

    Dimension mTimeDimension;
    Dimension mSpaceDimensionX;
    Dimension mSpaceDimensionY;

    DoubleVector vmin;
    DoubleVector vmax;

    double regEpsilon;
    double r;
public:
    Problem2Forward2D *forward;
    Problem2Backward2D *backward;
    DoubleMatrix U;

    bool optimizeK;
    bool optimizeZ;
    bool optimizeC;
    bool optimizeO;
};

#endif // IFUNCTIONAL_H
