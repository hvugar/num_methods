#ifndef HYPERBOLICCONTROLX_H
#define HYPERBOLICCONTROLX_H

#include <stdio.h>
#include <math.h>
#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <tomasmethod.h>

class HyperbolicControlX : public R1Function, public RnFunction, public IGradient, public IHyperbolicEquation, public IBackwardHyperbolicEquation, public IPrinter
{
public:
    HyperbolicControlX();
    virtual ~HyperbolicControlX() {}

    virtual double fx(double t);
    virtual double fx(const DoubleVector& v);
    virtual void gradient(const DoubleVector &v, DoubleVector &g);

    virtual double fi1(unsigned int i) const;
    virtual double fi2(unsigned int i) const;
    virtual double m1(unsigned int j) const;
    virtual double m2(unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    virtual double bfi1(unsigned int i) const { return 0.0; }
    virtual double bfi2(unsigned int i) const { return 0.0; }
    virtual double bm1(unsigned int j) const { return 0.0; }
    virtual double bm2(unsigned int j) const { return 0.0; }
    virtual double bf(unsigned int i, unsigned int j) const
    {
        const DoubleMatrix &u = *pu;
        if (M<=j && j<=M+D-2)
        {
            return -(2.0*(u[j][i]-U));
        }
        return 0.0;
    }

    virtual void print(unsigned int iteration, const DoubleVector& v, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

//    void calculateP(const DoubleMatrix &u, const DoubleVector &v, DoubleVector &g);
//    void calculateG(const DoubleVector &p, DoubleVector &g,  unsigned int j);

    static void main();
    void gradient1(const DoubleVector &v, DoubleVector &g);

//    inline virtual double pfi1(unsigned int i) const { return 0.0; }
//    inline virtual double pfi2(unsigned int i) const { return 0.0; }
//    inline virtual double pm1(unsigned int j) const { return 0.0; }
//    inline virtual double pm2(unsigned int j) const { return 0.0; }

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double hx;
    double ht;
    double a;
    unsigned int N;
    unsigned int M;
    unsigned int D;
    unsigned int L;
    unsigned int Xi;
    double xi;
    const DoubleVector *pv;
    const DoubleMatrix *pu;
    double U;
    double lamda;
};

#endif // HYPERBOLICCONTROLX_H
