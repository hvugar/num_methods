#ifndef ABSTRACTPROBLEM22D_H
#define ABSTRACTPROBLEM22D_H

#include "problem2setting.h"
#include "iproblem2forward2d.h"
#include "iproblem2backward2d.h"

#include <function.h>
#include <gradient.h>
#include <utils/random.h>

class AProblem2Forward2D;
class AProblem2Backward2D;
class AbstactProblem22D;

class AProblem2Forward2D : public IProblem2Forward2D
{
public:
    AbstactProblem22D *p22d;

protected:
    virtual double initial(const SpaceNodePDE &) const { return 0.0; }
    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const { return NAN; }
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

    //virtual void layerInfo(const DoubleMatrix &u, unsigned int layerNumber) const;
};

class AProblem2Backward2D : public IProblem2Backward2D
{
public:
    AbstactProblem22D *p22d;

protected:
    virtual double initial(const SpaceNodePDE &) const { return 0.0; }
    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const { return NAN; }
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

    virtual double h(const SpaceNodePDE &) const { return 0.0; }
    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

    //virtual void layerInfo(const DoubleMatrix &p, unsigned int layerNumber) const;
};

class AbstactProblem22D : public RnFunction, public IGradient
{
public:
    void setGridParameters(Dimension timeDimension, Dimension spaceDimensionX, Dimension spaceDimensionY);

    virtual double fx(const DoubleVector &prms) const;
    virtual void gradient(const DoubleVector &prms, DoubleVector &g);

    const P2Setting& setting() const;
    void setP2Setting(const P2Setting& setting);

    virtual double integral(const DoubleMatrix &u) const;
    virtual double mu(double x, double y) const;

    void setForward(IProblem2Forward2D *);
    void setBackward(IProblem2Backward2D *);
private:
    Dimension mTimeDimension;
    Dimension mSpaceDimensionX;
    Dimension mSpaceDimensionY;

    P2Setting msetting;

    double alpha0;

    IProblem2Forward2D *forward;
    IProblem2Backward2D *backward;

public:
    DoubleMatrix U;
    DoubleVector *grad;
};


#endif // ABSTRACTPROBLEM22D_H
