#ifndef ROSENBROCK_H
#define ROSENBROCK_H

#include <function.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>

class MINIMUMSHARED_EXPORT Rosenbrock : public RnFunction, public IGradient, public IPrinter
{
public:
    Rosenbrock();
    virtual ~Rosenbrock() {}
    //RnFunction
    virtual double fx(const DoubleVector& x) const;
    //IGradient
    virtual void gradient(const DoubleVector &x, DoubleVector &g) const;
    //Printer
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double, GradientMethod::MethodResult result) const;

    static void Main(int argc, char *argv[]);
private:
    double grad_step;
    unsigned int count;
};

#endif // ROSENBROCK_H
