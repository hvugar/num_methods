#ifndef FUNCTION_H
#define FUNCTION_H

#include "global.h"
#include "vector2d.h"

#ifdef __cplusplus
extern "C" {
#endif

double MINIMUMSHARED_EXPORT sgn(double x);

#ifdef __cplusplus
}
#endif

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

class MINIMUMSHARED_EXPORT VectorFunction
{
public:
    virtual double fx(double x, unsigned int num) const = 0;
};

class MINIMUMSHARED_EXPORT VectorRnFunction
{
public:
    virtual double fx(const DoubleVector &x, unsigned int num = 0) const = 0;
};

class MINIMUMSHARED_EXPORT MatrixFunction
{
public:
    virtual double fx(double x, unsigned int row, unsigned int col) const = 0;
};

class MINIMUMSHARED_EXPORT MatrixRnFunction
{
public:
    virtual double fx(const DoubleVector &x, unsigned int row, unsigned int col) const = 0;
};

class MINIMUMSHARED_EXPORT IGradient
{
public:
    virtual void gradient(const DoubleVector &x, DoubleVector &g) const = 0;

public:
    static void Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g);
    static void Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g, unsigned int start, unsigned int end);
    static void Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g, unsigned int *inx, unsigned int size);
};

//struct MINIMUMSHARED_EXPORT OrdDifEquation
//{
//    virtual double fx(double x, double y) const = 0;
//};

//typedef std::vector<RnFunction> RnFunctionList;
//typedef std::vector<RnFunction*> PRnFunctionList;

#endif // FUNCTION_H
