#ifndef PROBLEM22D_H
#define PROBLEM22D_H

#include "problem2setting.h"
#include "iproblem2forward2d.h"
#include "iproblem2backward2d.h"

#include <function.h>
#include <gradient.h>
#include <utils/random.h>

class Problem2Forward2D : public IProblem2Forward2D
{
protected:
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double initial(const SpaceNodePDE &) const { return 0.0; }
};

class Problem2Backward2D : public IProblem2Backward2D
{
protected:
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double h(const SpaceNodePDE &) const { return 0.0; }
    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
};

class Problem22D : public RnFunction, public IGradient
{
public:
    static void Main(int argc, char* argv[]);

    void setGridParameters(Dimension timeDimension, Dimension spaceDimensionX, Dimension spaceDimensionY);
    virtual ~Problem22D() {}

    virtual double fx(const DoubleVector &x) const;
    virtual void gradient(const DoubleVector &x, DoubleVector &g);

    virtual double integral(const DoubleMatrix &u) const;
    virtual double mu(double x, double y) const;

    DoubleMatrix U;
    P2Setting setting;

    static void testExample1();

private:
    Dimension mTimeDimension;
    Dimension mSpaceDimensionX;
    Dimension mSpaceDimensionY;

    double alpha0;

    void testForwardEquation(const P2Setting &setting) const;
    void testBackwardEquation(const P2Setting &setting) const;
};

#endif // PROBLEM22D_H
