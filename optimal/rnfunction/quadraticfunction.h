#ifndef QUADRATICFUNCTION_H
#define QUADRATICFUNCTION_H

#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>
#include <printer.h>
#include <cmath>

class MINIMUMSHARED_EXPORT QuadraticFunction : public RnFunction, public IGradient, public IPrinter
{
public:
    QuadraticFunction();
    virtual ~QuadraticFunction() {}
    //RnFunction
    virtual double fx(const DoubleVector& x) const;
    //IGradient
    virtual void gradient(const DoubleVector &x, DoubleVector &g) const;
    //IPrinter
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, GradientMethod::MethodResult result) const;

    static void Main(int argc, char *argv[]);
private:
    double grad_step;
    unsigned int count;
};

#endif // QUADRATICFUNCTION_H
