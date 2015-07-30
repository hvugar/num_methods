#ifndef R1FUNCTION_H
#define R1FUNCTION_H

#include "doublevector.h"
#include "global.h"

using namespace std;

struct MINIMUMSHARED_EXPORT R1Function
{
    virtual double fx(double x) const = 0;
};

struct MINIMUMSHARED_EXPORT R2Function
{
    virtual double fx(double x, double y) const = 0;
};

struct MINIMUMSHARED_EXPORT RnFunction
{
    virtual double fx(const DoubleVector& x) const = 0;
    virtual void gradient(DoubleVector& g) const = 0;
};

typedef std::vector<RnFunction> RnFunctionList;
typedef std::vector<RnFunction*> PRnFunctionList;

#endif // R1FUNCTION_H
