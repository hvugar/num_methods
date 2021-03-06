#ifndef HYPERBOLICCONTROL2D22_H
#define HYPERBOLICCONTROL2D22_H

#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <projection.h>

class MINIMUMSHARED_EXPORT HyperbolicControl2D22 : public R1Function, public RnFunction,
        public IHyperbolicEquation2D, public IBackwardHyperbolicEquation2D,
        public IGradient, public IPrinter//, public IProjection
{
public:
    HyperbolicControl2D22();
    virtual ~HyperbolicControl2D22() {}

    virtual double fx(double x) const;
    virtual double fx(const DoubleVector &x) const;

    virtual void gradient(const DoubleVector &v, DoubleVector &g);
    virtual void print(unsigned int i, const DoubleVector &v, const DoubleVector &g, double fx) const;

    virtual double initial1(unsigned int i, unsigned int j) const;
    virtual double initial2(unsigned int i, unsigned int j) const;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    virtual double binitial1(unsigned int i, unsigned int j) const;
    virtual double binitial2(unsigned int i, unsigned int j) const;
    virtual double bboundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const;

public:
    static void Main(int argc, char *argv[]);

    double fxt(unsigned int i, unsigned int j, unsigned int k) const;

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
    DoubleVector x;
    const DoubleVector *pv;
    const DoubleCube *pu;
    DoubleVector v0;

    FILE *file;
};

#endif // HYPERBOLICCONTROL2D22_H
