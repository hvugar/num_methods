#ifndef FUNCTION_H
#define FUNCTION_H

//#include "doublevector.h"
#include "global.h"
#include "vector2d.h"

using namespace std;

struct MINIMUMSHARED_EXPORT R1Function
{
    virtual double fx(double x) = 0;
};

struct MINIMUMSHARED_EXPORT R2Function
{
    virtual double fx(double x, double y) = 0;
};

struct MINIMUMSHARED_EXPORT R3Function
{
    virtual double fx(double x, double y, double z) = 0;
};

class MINIMUMSHARED_EXPORT RnFunction
{
public:
    virtual double fx(const DoubleVector &x) = 0;
};

class MINIMUMSHARED_EXPORT IGradient
{
public:
    virtual void gradient(const DoubleVector &x, DoubleVector &g) = 0;

protected:
    static void Gradient(RnFunction *f, double step, const DoubleVector &x, DoubleVector &g);
};

struct MINIMUMSHARED_EXPORT OrdDifEquation
{
    virtual double fx(double x, double y) const = 0;
};

//typedef std::vector<RnFunction> RnFunctionList;
//typedef std::vector<RnFunction*> PRnFunctionList;

#endif // FUNCTION_H
