#ifndef HYPERBOLICCONTROL2D24_H
#define HYPERBOLICCONTROL2D24_H

#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>

class HyperbolicControl2D24 : public R1Function, public RnFunction,
        public IHyperbolicEquation2D, public IBackwardHyperbolicEquation2D,
        public IGradient, public IPrinter, public IProjection
{
public:
    HyperbolicControl2D24();
    virtual ~HyperbolicControl2D24() {}

    virtual double fx(double x);
    virtual double fx(const DoubleVector &x);

    virtual void gradient(const DoubleVector &v, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;
    virtual void project(DoubleVector &x, int index);

    virtual double initial1(unsigned int i, unsigned int j) const;
    virtual double initial2(unsigned int i, unsigned int j) const;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    virtual double binitial1(unsigned int i, unsigned int j) const;
    virtual double binitial2(unsigned int i, unsigned int j) const;
    virtual double bboundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const;

public:
    static void main(int argc, char ** argv);

    double fxt(unsigned int i, unsigned int j, unsigned int k) const;
    void psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi);
    void printGradients(const DoubleVector &x, unsigned int i, FILE* f) const;

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

    double U0;
    double U1;

    DoubleVector e;
    const DoubleVector *px;
    const DoubleCube *pu;
    FILE *file;

    DoubleVector X;
};

#endif
