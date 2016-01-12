#ifndef DISCRETEHYPERBOLIC_H
#define DISCRETEHYPERBOLIC_H

#include <function.h>
#include <doublevector.h>
#include <hyperbolicequation.h>
#include <printer.h>

class DiscreteHyperbolic : public RnFunction, public IGradient, public IHyperbolicEquation, public Printer
{
public:
    DiscreteHyperbolic();
    virtual ~DiscreteHyperbolic() {}

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g);

    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    virtual double fi1(unsigned int i) const;
    virtual double fi2(unsigned int i) const;
    virtual double m1(unsigned int j) const;
    virtual double m2(unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    double F(unsigned int i, unsigned int j) const;
    void calculateP(const DoubleVector& f0, const DoubleMatrix &u, DoubleMatrix &psi, DoubleVector &g);

    //virtual void calculateP();

    static void main();

private:
    double t0;
    double t1;
    double x0;
    double x1;
    unsigned int N;
    unsigned int M;
    double ht;
    double hx;
    double a;
    double lamda;
    DoubleVector U;
    const DoubleVector *pf;
};

#endif // DISCRETEHYPERBOLIC_H
