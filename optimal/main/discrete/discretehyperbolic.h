#ifndef DISCRETEHYPERBOLIC_H
#define DISCRETEHYPERBOLIC_H

#include <function.h>
#include <doublevector.h>
#include <hyperbolicequation.h>
#include <printer.h>

class Discretehyperbolic : public RnFunction, public IHyperbolicEquation, public Printer
{
public:
    Discretehyperbolic();
    virtual ~Discretehyperbolic() {}

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step=0.000001);

    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    virtual double fi1(unsigned int i) const;
    virtual double fi2(unsigned int i) const;
    virtual double m1(unsigned int j) const;
    virtual double m2(unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    double F(unsigned int i, unsigned int j) const;
    void calculateP(const DoubleVector &u, DoubleMatrix &psi);

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
