#ifndef HYPERBOLICCONTROL2D1_H
#define HYPERBOLICCONTROL2D1_H

//#define ONLY_POWER
//#define ONLY_COORDINATE
#define POWER_COORDINATE

#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>

class MINIMUMSHARED_EXPORT HyperbolicControl2D1 : public R1Function, public RnFunction, public IGradient, public IHyperbolicEquation2D,
        public IBackwardHyperbolicEquation2D, public IPrinter, public IProjection
{
public:
    HyperbolicControl2D1();
    virtual ~HyperbolicControl2D1() {}

    virtual double fx(double x) const;
    virtual double fx(const DoubleVector &x) const;
    virtual void gradient(const DoubleVector &v, DoubleVector &g);

    virtual double initial1(unsigned int i, unsigned int j) const;
    virtual double initial2(unsigned int i, unsigned int j) const;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    virtual double binitial1(unsigned int i, unsigned int j) const;
    virtual double binitial2(unsigned int i, unsigned int j) const;
    virtual double bboundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const;

    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &g, double fx) const;
    virtual void project(DoubleVector &x, int index);


public:
    static void Main(int argc, char* argv[]);

    double u(double i, double j, double k) const;
    void psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi);
    double fxt(unsigned int i, unsigned int j, unsigned int k) const;
    double norm(const DoubleVector& v) const;

    double v1(double t) const;
    double v2(double t) const;

    double t0;
    double t1;
    double x10;
    double x11;
    double x20;
    double x21;
    double h1;
    double h2;
    double ht;
    unsigned int L;
    unsigned int N1;
    unsigned int N2;
    unsigned int M;

    double alpha0;
    double alpha1;
    double alpha2;
    double alpha3;
    double qamma;
    double a1;
    double a2;

    DoubleMatrix U0;
    DoubleMatrix U1;

    DoubleVector e;
    DoubleVector c;
    const DoubleVector *px;
    const DoubleCube *pu;
    double h;

    double vd;
    double vu;
    unsigned int count;
};

#endif // HYPERBOLICCONTROL2D1_H
