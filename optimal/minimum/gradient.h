#ifndef GRADIENT_H
#define GRADIENT_H

#include "global.h"

#include <vector>
#include <math.h>
#include <stdio.h>

#include "function.h"
#include "r1minimize.h"

using namespace std;

class MINIMUMSHARED_EXPORT Gradient
{
public:
    Gradient();
    virtual ~Gradient();

    virtual RnFunction* function() const;
    virtual void setFunction(RnFunction* function);

    virtual const std::vector<double>& x() const;
    virtual void setX(const std::vector<double>& x);

    virtual void calculate() = 0;

    double epsilon() const;
    void setEpsilon(double epsilon);

    void setR1MinimizeEpsilon(double step, double epsilon);
    void setGradientStep(double step);
    int count() const;

protected:
    virtual double minimize() = 0;
    virtual void calculateGradient();
    virtual double gradientNorm() const;
    virtual double distance() const;

    RnFunction *m_fn;
    std::vector<double> m_x;
    std::vector<double> m_g;
    double m_alpha;
    double m_epsilon;
    double grad_step;
    double min_epsilon;
    double min_step;
    int k;
    int M;
};

#endif // GRADIENT_H
