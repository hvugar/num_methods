#ifndef FUNCTION_H
#define FUNCTION_H

//#include "doublevector.h"
#include "global.h"
#include "vector2d.h"

struct MINIMUMSHARED_EXPORT R1Function
{
    virtual double fx(double x) const = 0;
};

struct MINIMUMSHARED_EXPORT R2Function
{
    virtual double fx(double x, double y) const = 0;
};

struct MINIMUMSHARED_EXPORT R3Function
{
    virtual double fx(double x, double y, double z) const = 0;
};

class MINIMUMSHARED_EXPORT RnFunction
{
public:
    virtual double fx(const DoubleVector &x) const = 0;
};

class MINIMUMSHARED_EXPORT IGradient
{
public:
    virtual void gradient(const DoubleVector &x, DoubleVector &g) = 0;

protected:
    static void Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g);
    static void Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g, unsigned int *inx, unsigned int size);
};

struct MINIMUMSHARED_EXPORT OrdDifEquation
{
    virtual double fx(double x, double y) const = 0;
};

//typedef std::vector<RnFunction> RnFunctionList;
//typedef std::vector<RnFunction*> PRnFunctionList;

#endif // FUNCTION_H
