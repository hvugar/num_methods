#ifndef PROBLEM0H_SOLVER_H
#define PROBLEM0H_SOLVER_H

#include <grid/hibvp.h>
#include <benchmark.h>
#include <function.h>
#include "problem0h_global.h"

#include <functional>

//int count1 = 0;
//int count2 = 0;

class Problem0HCommon;
class Problem0HForward;
class Problem0HBckward;
class Problem0HFunctional;

class Problem0HCommon
{
public:
    Problem0HCommon();
    virtual ~Problem0HCommon();

protected:
    inline auto virtual mu(const SpaceNodePDE &) const -> double { return 0.0; }
    inline auto virtual mu(unsigned int, unsigned int) const -> double { return 1.0; }

    double a = 1.0;
    double gamma = 0.0;
    double alpha0 = 1.0;
    double alpha1 = 1.0;
    double alpha2 = 1.0;
    double alpha3 = 1.0;
    SpacePoint ksi;

    DoubleMatrix U1;
    DoubleMatrix U2;
    DoubleMatrix u1;
    DoubleMatrix u2;

    struct Psi {
        double vl;
        double dx;
        double dy;
    };

    std::vector<double> v1;
    SpacePoint p1;
    std::vector<double> v2;
    SpacePoint p2;
    std::vector<Psi> ps1;
    std::vector<Psi> ps2;
};

class Problem0HForward : public CdIHyperbolicIBVP, public virtual Problem0HCommon
{
public:
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

private:
    double p(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
};

class Problem0HBckward : public ConjugateCdIHyperbolicIBVP, public virtual Problem0HCommon
{
public:
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
};

class Problem0HFunctional : public RnFunction, public IGradient,
        public virtual Problem0HForward, public virtual Problem0HBckward
{
public:
    static void Main(int argc, char** argv);

public:
    virtual auto fx(const DoubleVector &x) const -> double;
    virtual auto gradient(const DoubleVector &x, DoubleVector &g) const -> void;

    virtual auto integral1(const DoubleMatrix &u) const -> double;
    virtual auto integral2(const DoubleMatrix &u) const -> double;
    virtual auto norm() const -> double;
    virtual auto penalty() const -> double;

    void setDimension(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY);

private:
    Problem0HFunctional* const_this = nullptr;
};


#endif // PROBLEM0H_SOLVER_H
