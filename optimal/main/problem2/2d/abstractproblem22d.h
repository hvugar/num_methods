#ifndef ABSTRACTPROBLEM22D_H
#define ABSTRACTPROBLEM22D_H

#include "iproblem2forward2d.h"
#include "iproblem2backward2d.h"

#include <function.h>
#include <gradient.h>
#include <utils/random.h>

#include <imaging.h>

class Problem2Forward2DEx4;
class Problem2Backward2DEx4;
class AbstactProblem22D;

//---------------------------------------------------------------------------------------------------------------//

class Problem2Forward2DEx4 : public IProblem2Forward2D
{
public:
    virtual ~Problem2Forward2DEx4() {}
    double fi;

protected:
    virtual double initial(const SpaceNodePDE &) const;
    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const;
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const;

    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const;

    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;
};

//---------------------------------------------------------------------------------------------------------------//


class Problem2Backward2DEx4 : public IProblem2Backward2D
{
public:
    virtual ~Problem2Backward2DEx4() {}

    AbstactProblem22D *ap22d;
    DoubleMatrix *U;
    DoubleMatrix *u;

protected:
    virtual double initial(const SpaceNodePDE &) const;
    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const;
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const;

    virtual double h(const SpaceNodePDE &) const;
    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const;

    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;
};

//---------------------------------------------------------------------------------------------------------------//

class AbstactProblem22D : public RnFunction, public IGradient
{
public:
    AbstactProblem22D();
    virtual ~AbstactProblem22D();

    void setGridParameters(Dimension timeDimension, Dimension spaceDimensionX, Dimension spaceDimensionY);

    virtual double fx(const DoubleVector &prms) const;
    virtual void gradient(const DoubleVector &prms, DoubleVector &g);

    const Parameter& parameter() const;
    void setParameter(const Parameter& parameter);

    virtual double integral(const DoubleMatrix &u) const;
    virtual double mu(double x, double y) const;

protected:
    Dimension mTimeDimension;
    Dimension mSpaceDimensionX;
    Dimension mSpaceDimensionY;

    Parameter mParameter;

    double alpha0;

    Problem2Forward2DEx4 *forward;
    Problem2Backward2DEx4 *backward;

public:
    DoubleMatrix U;
};


#endif // ABSTRACTPROBLEM22D_H
