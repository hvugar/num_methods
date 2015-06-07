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

    virtual RnFunction* f() const;
    virtual void setF(RnFunction* f);

    const std::vector<double>& x() const;
    void setX(const std::vector<double>& x);

    double epsilon() const;
    void setEpsilon(double epsilon);

    virtual void calcGradient();
    virtual double minimize() = 0;
    virtual void calculate() = 0;

    void setR1MinimizeEpsilon(double step, double epsilon);
    void setGradientStep(double step);
    int count() const;

protected:
    virtual double gradientNorm() const;
    virtual double distance() const;

    RnFunction *mfn;
    std::vector<double> mx;
    std::vector<double> mg;
    double malpha;
    double mepsilon;
    double grad_step;
    double min_epsilon;
    double min_step;
    int mcount;
};

#endif // GRADIENT_H
