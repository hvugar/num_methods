#ifndef DISCRETEHEAT_H
#define DISCRETEHEAT_H

#include <function.h>
#include <doublevector.h>
#include <gradient_cjt.h>
#include <parabolicequation.h>
#include <printer.h>

class DiscreteHeat : public RnFunction, public IGradient, public IPrinter, public IParabolicEquation
{
public:
    DiscreteHeat();
    virtual ~DiscreteHeat() {}

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g);

    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    virtual double initial(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    double fxt(double x, double t);

//    virtual double pfi(unsigned int i) const;
//    virtual double pm1(unsigned int j) const;
//    virtual double pm2(unsigned int j) const;
//    virtual double pf(unsigned int i, unsigned int j) const;

    void calculateP(const DoubleVector &f, const DoubleVector &u, DoubleMatrix &psi, DoubleVector &g);

    static void main();

private:
    double x0;
    double x1;
    double t1;
    double t0;
    unsigned int N;
    unsigned int M;
    double hx;
    double ht;
    double a;

    const DoubleVector *pf;
    DoubleVector U;
};

#endif // DISCRETEHEAT_H
