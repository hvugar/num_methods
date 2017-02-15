#ifndef HYPERBOLICCONTROL2D_H
#define HYPERBOLICCONTROL2D_H

#include <hyperbolicequation.h>
#include <function.h>
#include <gradient_cjt.h>
#include <printer.h>

class MINIMUMSHARED_EXPORT HyperbolicControl2D : public R1Function, public RnFunction, public IGradient, public IHyperbolicEquation2D, public IBackwardHyperbolicEquation2D, public IPrinter
{
public:
    HyperbolicControl2D();
    virtual ~HyperbolicControl2D();

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

    double v1(double t) const;
    double v2(double t) const;
    double v3(double t) const;
    double fxt(unsigned int i, unsigned int j, unsigned int k) const;

    static void Main(int argc, char* argv[]);

private:
    double ht;
    double h1;
    double h2;
    unsigned int M;
    unsigned int N1;
    unsigned int N2;
    unsigned int L;
    double a1;
    double a2;
    double alpha0;
    double alpha1;
    double qamma;

    DoubleMatrix U0;
    DoubleMatrix U1;
    DoubleVector E;
    const DoubleVector *pv;
    const DoubleCube *pu;

    inline double u(unsigned int i, unsigned int j, unsigned int k) const;
    double norm(const DoubleVector& v) const;
};

#endif // HYPERBOLICCONTROL2D_H
