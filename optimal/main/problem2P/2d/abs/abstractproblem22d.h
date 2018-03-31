#ifndef ABSTRACTPROBLEM22D_H
#define ABSTRACTPROBLEM22D_H

#include <function.h>
#include <gradient.h>
#include <gradient_cjt.h>
#include <projection.h>
#include <printer.h>
#include <utils/random.h>
#include "problem2forward2d.h"
#include "problem2backward2d.h"

#include <imaging.h>

class AbstactProblem22D : public RnFunction, public IGradient, public IProjection, public IPrinter
{
public:
    AbstactProblem22D();
    virtual ~AbstactProblem22D();

    void setGridParameters(Dimension timeDimension, Dimension spaceDimensionX, Dimension spaceDimensionY);
    void setEquationParameters(double a, double lambda0, double lambda);
    void setIntTemperature(double fi);
    void setEnvTemperature(double theta);
    void setEpsilon(double epsilon);
    void setPenaltyCoefficient(double r);
    void setPenaltyLimits(const DoubleVector &vmin, const DoubleVector &vmax);
    void setParameter(const Parameter &parameter);
    void setParameter0(const Parameter &parameter0);

    virtual double fx(const DoubleVector &prms) const;
    virtual double integral(const DoubleMatrix &u) const;
    virtual double norm() const;
    virtual double penalty(const vector<ExtendedSpaceNode2D> &info) const;
    virtual double g0i(unsigned int i, unsigned int layer, const vector<ExtendedSpaceNode2D> &info) const;
    virtual double gpi(unsigned int i, unsigned int layer, const vector<ExtendedSpaceNode2D> &info) const;

    virtual void gradient(const DoubleVector &prms, DoubleVector &g);

    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, GradientMethod::MethodResult result) const;
    virtual void project(DoubleVector &x, unsigned int index);

    virtual double mu(double x, double y) const;

    void optimization(DoubleVector &prm0);

    void calculateU();

    double sgn(double x) const;

protected:
    double epsilon;
    double r;

    DoubleVector vmin;
    DoubleVector vmax;

    Dimension mTimeDimension;
    Dimension mSpaceDimensionX;
    Dimension mSpaceDimensionY;

    Parameter mParameter;
    Parameter mParameter0;

    Problem2Forward2D *forward;
    Problem2Backward2D *backward;

public:
    DoubleMatrix U;
};

#endif // ABSTRACTPROBLEM22D_H
