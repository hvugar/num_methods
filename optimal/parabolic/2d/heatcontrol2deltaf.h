#ifndef HEATCONTROL2DELTAF_H
#define HEATCONTROL2DELTAF_H

#include <math.h>
#include <stdlib.h>

#include <function.h>
#include <pde/parabolicequation.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

class HeatControl2DeltaF : public RnFunction, public IGradient, public IParabolicEquation2D, public IBackwardParabolicEquation2D, public IPrinter, public IProjection
{
public:
    HeatControl2DeltaF(unsigned int m, unsigned int n2, unsigned int n1);
    virtual ~HeatControl2DeltaF() {}

    virtual double fx(const DoubleVector &v) const;
    virtual void gradient(const DoubleVector &v, DoubleVector &g);
    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &g, double fx, GradientMethod::MethodResult result) const;
    virtual void project(DoubleVector &v, int index);

    virtual double initial(unsigned int i, unsigned int j) const;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    virtual double binitial(unsigned int i, unsigned int j) const;
    virtual double bboundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const;

    void calculateGF(const DoubleVector &v, const DoubleMatrix& psi, DoubleVector& g, unsigned int k);

    static void main(int argc, char *argv[]);
private:
    double u(double x1, double x2, double t) const { return x1*x1 + x2*x2 + t*t; }
    double norm(const DoubleVector& v) const;

    inline double v1(double t) const { return 10*t; }
    inline double v2(double t) const { return 20*t; }
    inline double v3(double t) const { return 30*t; }

    inline void psiDerivative(double &psiX1, double &psiX2, double e1, double e2, const DoubleMatrix &psi);

    double t0;
    double t1;

    double x10;
    double x11;
    double x20;
    double x21;

    unsigned int M;
    unsigned int N1;
    unsigned int N2;
    unsigned int L;

    double a1;
    double a2;

    double h1;
    double h2;
    double ht;

    double alpha;
    double gause_a;
    double gause_b;

    DoubleMatrix U;
    const DoubleVector *pv;
    const DoubleMatrix *pu;
    DoubleVector E;
};

#endif // HEATCONTROL2DELTAF_H
