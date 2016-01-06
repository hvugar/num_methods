#ifndef DISCRETEHEAT_H
#define DISCRETEHEAT_H

#include <function.h>
#include <doublevector.h>
#include <gradient_cjt.h>
#include <parabolicequation.h>
#include <printer.h>

class DiscreteHeat : public RnFunction, public Printer, public IParabolicEquation, public IBackwardParabolicEquation
{
public:
    DiscreteHeat();
    virtual ~DiscreteHeat() {}

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step=0.000001);

    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    virtual double fi(unsigned int i) const;
    virtual double m1(unsigned int j) const;
    virtual double m2(unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    double fxt(double x, double t);

//    virtual double pfi(unsigned int i) const;
//    virtual double pm1(unsigned int j) const;
//    virtual double pm2(unsigned int j) const;
//    virtual double pf(unsigned int i, unsigned int j) const;

    void calculateP(const DoubleVector &f, const DoubleVector &u, DoubleMatrix &psi, DoubleVector &g);

    static void main();
    double hx;
private:
    double x0;
    double x1;
    double t1;
    double t0;
    unsigned int N;
    unsigned int M;

    double ht;
    double a;

    const DoubleVector *pf;
    DoubleVector U;
};

#endif // DISCRETEHEAT_H
