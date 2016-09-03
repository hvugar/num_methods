#ifndef POINTCONTROL11_H
#define POINTCONTROL11_H

#include <function.h>
#include <printer.h>
#include <gradient_cjt.h>

class PointControl11 : public RnFunction, public IPrinter
{
public:
    PointControl11();
    virtual ~PointControl11();

    virtual double fx(const DoubleVector &q);
    virtual void gradient(const DoubleVector &q, DoubleVector &g, double gradient_step);
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    void calculateX(const DoubleVector &q, DoubleVector &x);
    void calculateP(const DoubleVector &q, const DoubleVector& x, DoubleVector& p);

    static void main(int argc, char ** argv);

    double f(double t, double x) const;
    double pf(double t, double p, double x) const;

private:
    unsigned int N;

    double t0;
    double t1;
    double x0;
    double x1;
    double ht;
    double hx;

};

#endif // POINTCONTROL11_H
