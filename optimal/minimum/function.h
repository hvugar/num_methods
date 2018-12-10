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
    virtual ~R1Function();
    virtual double fx(double x) const = 0;
};

struct MINIMUMSHARED_EXPORT R2Function
{
    virtual ~R2Function();
    virtual double fx(double x, double y) const = 0;
};

struct MINIMUMSHARED_EXPORT R3Function
{
    virtual ~R3Function();
    virtual double fx(double x, double y, double z) const = 0;
};

struct MINIMUMSHARED_EXPORT RnFunction
{
    virtual ~RnFunction();
    virtual double fx(const DoubleVector &x) const = 0;
};

struct MINIMUMSHARED_EXPORT VectorFunction
{
    virtual ~VectorFunction();
    virtual double fx(double x, unsigned int num) const = 0;
};

struct MINIMUMSHARED_EXPORT VectorRnFunction
{
    virtual ~VectorRnFunction();
    virtual double fx(const DoubleVector &x, unsigned int num = 0) const = 0;
};

struct MINIMUMSHARED_EXPORT MatrixR1Function
{
    virtual ~MatrixR1Function();
    virtual double fx(double x, unsigned int row, unsigned int col) const = 0;
};

struct MINIMUMSHARED_EXPORT MatrixRnFunction
{
    virtual ~MatrixRnFunction();
    virtual double fx(const DoubleVector &x, unsigned int row, unsigned int col) const = 0;
};

struct MINIMUMSHARED_EXPORT IGradient
{
    virtual ~IGradient();
    virtual void gradient(const DoubleVector &x, DoubleVector &g) const = 0;

    static void Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g);
    static void Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g, unsigned int start, unsigned int end);
    static void Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g, unsigned int *inx, unsigned int size);
};

#endif // FUNCTION_H
