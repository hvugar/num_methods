#ifndef HYPERBOLICCONTROL2DM_H
#define HYPERBOLICCONTROL2DM_H

#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>

class HyperbolicControl2DM : public R1Function, public RnFunction,
        public IGradient, public IHyperbolicEquation2D, public IBackwardHyperbolicEquation2D,
        public IPrinter, public Projection
{
public:
    HyperbolicControl2DM();
    virtual ~HyperbolicControl2DM();

    virtual double fx(double x);
    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &v, DoubleVector &e);

    virtual double fi1(unsigned int i, unsigned int j) const;
    virtual double fi2(unsigned int i, unsigned int j) const;
    virtual double m1(unsigned int j, unsigned int k) const;
    virtual double m2(unsigned int j, unsigned int k) const;
    virtual double m3(unsigned int i, unsigned int k) const;
    virtual double m4(unsigned int i, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    virtual double bfi1(unsigned int i, unsigned int j) const;
    virtual double bfi2(unsigned int i, unsigned int j) const;
    virtual double bm1(unsigned int j, unsigned int k) const;
    virtual double bm2(unsigned int j, unsigned int k) const;
    virtual double bm3(unsigned int i, unsigned int k) const;
    virtual double bm4(unsigned int i, unsigned int k) const;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const;

    virtual void print(unsigned int i, const DoubleVector& x, const DoubleVector &e, double alpha, RnFunction* fn) const;
    virtual void project(DoubleVector &x, int index);
    double fxt(unsigned int i, unsigned int j, unsigned int k) const;
    double norm(const DoubleVector& v) const;

    void calculateGX(const DoubleVector& x, const DoubleMatrix& psi, DoubleVector& g, unsigned int k);
    void psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi);

    static void main();

private:
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
    double a;

    double U0;
    double U1;

    //DoubleVector d;
    DoubleVector e;
    const DoubleVector *px;
    const DoubleCube *pu;

    double vd;
    double vu;
};

#endif // HYPERBOLICCONTROL2DM_H
