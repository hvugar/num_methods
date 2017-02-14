#ifndef HEATCONTROL1_H
#define HEATCONTROL1_H

#include <math.h>
#include <stdlib.h>
#include <function.h>
#include <parabolicequation.h>
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

class MINIMUMSHARED_EXPORT HeatControl1 : public RnFunction, public IGradient,
        //public ParabolicIBVP, public BackwardParabolicIBVP,
        //public IParabolicEquation, public IBackwardParabolicEquation,
        public IPrinter
{
public:
    HeatControl1();
    virtual ~HeatControl1() {}

    class CParabolicIBVP : public ParabolicIBVP
    {
    public:
        virtual double initial(const SpaceNode &sn) const;
        virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const;
        virtual double f(const SpaceNode &sn, const TimeNode &tn) const;
        virtual double a(const SpaceNode &sn, const TimeNode &tn) const;
        const DoubleVector *pf;
    } pibvp;

    class CBackwardParabolicIBVP : public BackwardParabolicIBVP
    {
    public:
        virtual double initial(const SpaceNode &sn) const;
        virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const;
        virtual double f(const SpaceNode &sn, const TimeNode &tn) const;
        virtual double a(const SpaceNode &sn, const TimeNode &tn) const;
        const DoubleVector *pu;
        const DoubleVector *pU;
    } bpibvp;

    virtual double fx(const DoubleVector& f);
    virtual void gradient(const DoubleVector& f, DoubleVector &g);

    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &s, double a, RnFunction* f) const;

private:
    double ht;
    double hx;
    unsigned int N;
    unsigned int M;

    double u(double x, double t) const;
    double fxt(double x, double t);
    DoubleVector U;

public:
    static void Main(int argc, char *argv[]);
};

#endif // HEATCONTROL_H
