#ifndef ABSTRACTPROBLEM22D_H
#define ABSTRACTPROBLEM22D_H

#include "problem2setting.h"
#include "iproblem2forward2d.h"
#include "iproblem2backward2d.h"

#include <function.h>
#include <gradient.h>
#include <utils/random.h>

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

    DoubleMatrix U;
};

#endif // ABSTRACTPROBLEM22D_H
