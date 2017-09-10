#ifndef HEATCONTROL1_H
#define HEATCONTROL1_H

#include <math.h>
#include <stdlib.h>
#include <function.h>
#include <pde_old/parabolicequation.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

#include <grid/pibvp.h>
#include <grid/bpibvp.h>

/**
 * @brief The HeatControl struct
 * du/dt = a(d^2u/dx^2) + f(x,t);
 * u(x,0) = fi(x);
 * u(0,t) = m1(t);
 * u(l,t) = m2(t);
 */

class MINIMUMSHARED_EXPORT HeatControl1 : public RnFunction, public IGradient, public IPrinter
{
public:
    HeatControl1();
    virtual ~HeatControl1() {}

    class CParabolicIBVP : public ParabolicIBVP
    {
    protected:
        virtual double initial(const SpaceNodePDE &sn) const;
        virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
        virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
        virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    public:
        const DoubleVector *pf;
    } forward;

    class CBackwardParabolicIBVP : public BackwardParabolicIBVP
    {
    protected:
        virtual double initial(const SpaceNodePDE &sn) const;
        virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
        virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
        virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

        virtual void layerInfo(const DoubleVector &, unsigned int) const;

    public:
        const DoubleVector *pu;
        const DoubleVector *pU;
        const DoubleMatrix *pp;
    } backward;

    virtual double fx(const DoubleVector& f) const;
    virtual void gradient(const DoubleVector& f, DoubleVector &g);

    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &g, double fx, GradientMethod::MethodResult result) const;

private:
    double ht;
    double hx;
    unsigned int N;
    unsigned int M;

    double u(double x, double t) const;
    double fxt(double x, double t) const;
    DoubleVector U;

public:
    static void Main(int argc, char *argv[]);
};

#endif // HEATCONTROL1_H
