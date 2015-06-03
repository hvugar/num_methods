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
    virtual ~Gradient();

    const std::vector<double>& x() const;
    void setX(const std::vector<double>& x);
    double epsilon() const;
    void setEpsilon(double epsilon);

    virtual RnFunction* f() const = 0;
    virtual void setF(RnFunction* f) = 0;

protected:
    virtual void gradient() = 0;
    virtual double argmin(double alpha) = 0;
    virtual double minimize() = 0;
    virtual void calculate() = 0;

    std::vector<double> mx;
    double mepsilon;
};

#endif // GRADIENT_H
