#ifndef GRADIENT_H
#define GRADIENT_H

#include "global.h"

#include <vector>
#include <math.h>

#include "function.h"
#include "r1minimize.h"
#include "doublevector.h"
#include "printer.h"

/**
 * @brief The GradientMethod class
 */
class MINIMUMSHARED_EXPORT GradientMethod
{
public:
    GradientMethod();
    virtual ~GradientMethod();

    virtual void calculate(DoubleVector& x) = 0;

    virtual RnFunction* function() const;
    virtual void setFunction(RnFunction* function);

    double epsilon() const;
    void setEpsilon(double epsilon);

    void setR1MinimizeEpsilon(double step, double epsilon);
    void setGradientStep(double step);
    int count() const;

    void setPrinter(GrPrinter* printer);
    void setNormalize(bool normalize);

protected:
    virtual double minimize(const DoubleVector &x, const DoubleVector &g) = 0;

    RnFunction *m_fn;
    double m_epsilon;
    double grad_step;
    double min_epsilon;
    double min_step;
    int iterationCount;
    bool normalize;
    GrPrinter* printer;
};

#endif // GRADIENT_H
