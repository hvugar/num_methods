#ifndef DISCRETEHYPERBOLIC1_H
#define DISCRETEHYPERBOLIC1_H

#include <function.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <hyperbolicequation.h>
#include <printer.h>

class DiscreteHyperbolic1 : public RnFunction, public R1Function, public IGradient, public IHyperbolicEquation, public IPrinter
{
public:
    DiscreteHyperbolic1() {}
    virtual ~DiscreteHyperbolic1() {}

    virtual double fx(double t);

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &v, DoubleVector &g);

    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    virtual double initial1(unsigned int i) const;
    virtual double initial2(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    //double F(unsigned int i, unsigned int j) const;
    void calculateP(const DoubleVector& v, const DoubleMatrix &u, DoubleMatrix &psi, DoubleVector &g);

    static void main();

private:
    double t0;
    double t1;
    double x0;
    double x1;
    unsigned int N;
    unsigned int M;
    unsigned int D;
    unsigned int L;
    unsigned int X1;
    double ht;
    double hx;
    double a;
    double lamda;
    double U;
    const DoubleVector *pv;
};

#endif // DISCRETEHYPERBOLIC_H
