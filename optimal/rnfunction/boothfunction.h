#ifndef BOOTHFUNCTION_H
#define BOOTHFUNCTION_H

#include <function.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>

/**
 * @brief The Booth's function. Range -10.0 <= x,y <= +10.0. Optimal f(1, 3)=0;
 */
class MINIMUMSHARED_EXPORT BoothFunction : public RnFunction, public IGradient, public IPrinter, public IProjection
{
public:
    virtual ~BoothFunction() {}
    //RnFunction
    virtual double fx(const DoubleVector& x) const;
    //IGradient
    virtual void gradient(const DoubleVector &x, DoubleVector &g) const;
    //Printer
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, GradientMethod::MethodResult result) const;
    //IProjection
    virtual void project(DoubleVector &x, unsigned int index);

    static void main(int argc, char ** argv);

private:
    double grad_step;
    double a;
    double b;
};

#endif // BOOTHFUNCTION_H
