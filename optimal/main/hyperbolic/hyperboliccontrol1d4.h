#ifndef HYPERBOLICCONTROL1D4_H
#define HYPERBOLICCONTROL1D4_H

#include <function.h>
#include <projection.h>
#include <printer.h>
#include <doublevector.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <hyperbolicequation.h>
#include <tomasmethod.h>
#include <r1minimize.h>
#include <stdlib.h>

class HyperbolicControl1D4 : public RnFunction, public R1Function, public Printer, public Projection
{
public:
    HyperbolicControl1D4();
    virtual ~HyperbolicControl1D4() {}

    void calculateSettings();

    virtual double fx(const DoubleVector& v);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step=0.000001);

    virtual double fx(double x);

    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;
    virtual void project(DoubleVector &v, int j);

    virtual double fi1(unsigned int i) const;
    virtual double fi2(unsigned int i) const;
    virtual double m1(unsigned int j) const;
    virtual double m2(unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    double pfi1(unsigned int i) const { return ((*pu)[M+D][i]-U)*hx; }
    double pfi2(unsigned int i) const { return 0.0; }
    double pmu1(unsigned int t) const { return 0.0; }
    double pmu2(unsigned int t) const { return 0.0; }

    void calculateU(const DoubleVector& v, DoubleMatrix &u);
    void calculateP(const DoubleMatrix &u, DoubleVector &g);
    void calculateG(const DoubleVector& psi, DoubleVector& g, unsigned int j);

    //double f1(double t) const { return U; }
    //double f2(double t) const { return 2.0*t*t; }
    //double f3(double t) const { return 4.0*t*t; }

    static void main();

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double ht;
    double hx;
    unsigned int M;
    unsigned int N;
    unsigned int D;
    unsigned int L;

    double U;
    DoubleVector e;

    double a;
    double lamda;
    double R;

    const DoubleVector *pv;
    const DoubleMatrix *pu;

    unsigned int count;
};

#endif // HYPERBOLICCONTROL1D4_H
