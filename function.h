#ifndef FUNCTION_H
#define FUNCTION_H

#include "doublevector.h"

using namespace std;

struct R1Function
{
    virtual double fx(double x) = 0;
};

class RnFunction
{
public:
    virtual double fx(const DoubleVector &x) = 0;
};

class IGradient
{
public:
    virtual void gradient(const DoubleVector &x, DoubleVector &g) = 0;

protected:
    static void Gradient(RnFunction *f, double step, const DoubleVector &x, DoubleVector &g);
};

#endif // FUNCTION_H
