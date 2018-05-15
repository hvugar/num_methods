#ifndef BEALESFUNCTION_H
#define BEALESFUNCTION_H

#include <function.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>

/**
 * @brief The Beales Function. Range -4.5 <= x,y <= 4.5. Optimal f(3, 0.5)=0;
 */
class MINIMUMSHARED_EXPORT BealesFunction : public RnFunction, public IGradient, public IPrinter, public IProjection
{
public:
    virtual ~BealesFunction() {}
    //RnFunction
    virtual double fx(const DoubleVector &x) const;
    //IGradient
    virtual void gradient(const DoubleVector &x, DoubleVector &g) const;
    //Printer
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, GradientMethod::MethodResult result) const;
    //Projection
    virtual void project(DoubleVector &x, unsigned int index);

    static void main(int argc, char ** argv);

    double a;
    double b;

private:
    double grad_step;
};

#endif // BEALESFUNCTION_H
