#ifndef GRADIENT_H
#define GRADIENT_H

#include "global.h"

#include <vector>
#include <math.h>
#include <stdio.h>

#include "function.h"
#include "r1minimize.h"
#include "doublevector.h"
#include "printer.h"

class MINIMUMSHARED_EXPORT GradientMethod
{
public:
    GradientMethod();
    virtual ~GradientMethod();

    virtual RnFunction* function() const;
    virtual void setFunction(RnFunction* function);

    virtual const DoubleVector& x() const;
    virtual void setX(const DoubleVector& x);

    virtual void calculate() = 0;

    double epsilon() const;
    void setEpsilon(double epsilon);

    void setR1MinimizeEpsilon(double step, double epsilon);
    void setGradientStep(double step);
    int count() const;

    void setPrinter(GrPrinter* printer);

protected:
    virtual double minimize() = 0;
    virtual void calculateGradient();
    virtual double gradientNorm() const;
    virtual double distance() const;

    RnFunction *m_fn;
    DoubleVector m_x;
    DoubleVector m_g;
    double m_alpha;
    double m_epsilon;
    double grad_step;
    double min_epsilon;
    double min_step;
    int iterationCount;
    int M;
    GrPrinter* printer;
};

#endif // GRADIENT_H
