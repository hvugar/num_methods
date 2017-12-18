#ifndef ABSTRACTPROBLEM22D_H
#define ABSTRACTPROBLEM22D_H

#include <function.h>
#include <gradient.h>
#include <gradient_cjt.h>
#include <projection.h>
#include <printer.h>
#include <utils/random.h>
#include "problem2forward2dex4.h"
#include "problem2backward2dex4.h"

#include <imaging.h>

class AbstactProblem22D : public RnFunction, public IGradient, public IProjection, public IPrinter
{
public:
    AbstactProblem22D();
    virtual ~AbstactProblem22D();

    void setGridParameters(Dimension timeDimension, Dimension spaceDimensionX, Dimension spaceDimensionY);

    virtual double fx(const DoubleVector &prms) const;
    virtual void gradient(const DoubleVector &prms, DoubleVector &g);

    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, GradientMethod::MethodResult result) const;
    virtual void project(DoubleVector &x, unsigned int index);

    const Parameter& parameter() const;
    void setParameter(const Parameter& parameter);

    virtual double integral(const DoubleMatrix &u) const;
    virtual double norm() const;
    virtual double mu(double x, double y) const;

    void optimization(DoubleVector &prm0);
    void calculateU(const Parameter &prm0);

public:
    double espilon;

protected:
    Dimension mTimeDimension;
    Dimension mSpaceDimensionX;
    Dimension mSpaceDimensionY;

    Parameter mParameter;
    Parameter mParameter0;

    double alpha0;

    Problem2Forward2DEx4 *forward;
    Problem2Backward2DEx4 *backward;

public:
    DoubleMatrix U;
};

#endif // ABSTRACTPROBLEM22D_H
