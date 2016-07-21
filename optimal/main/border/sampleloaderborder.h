#ifndef SAMPLELOADERBORDER_H
#define SAMPLELOADERBORDER_H

#include <function.h>
#include <doublevector.h>
#include <rungekutta.h>
#include <printer.h>
#include "cmethods.h"

double A11(double t);
double A12(double t);
double A21(double t);
double A22(double t);

double B11(double t);
double B12(double t);
double B21(double t);
double B22(double t);

double C1(double t);
double C2(double t);

struct RnFunctionR0 : public RnFunction { virtual double fx(const DoubleVector &x); };
struct RnFunctionS0 : public RnFunction { virtual double fx(const DoubleVector &x); };

struct RnFunctionAlpha1_11 : public RnFunction { virtual double fx(const DoubleVector &x); };
struct RnFunctionAlpha1_12 : public RnFunction { virtual double fx(const DoubleVector &x); };
//struct RnFunctionAlpha1_21 : public RnFunction { virtual double fx(const DoubleVector &x); };
//struct RnFunctionAlpha1_22 : public RnFunction { virtual double fx(const DoubleVector &x); };

//struct RnFunctionAlpha2_11 : public RnFunction { virtual double fx(const DoubleVector &x); };
//struct RnFunctionAlpha2_12 : public RnFunction { virtual double fx(const DoubleVector &x); };
//struct RnFunctionAlpha2_21 : public RnFunction { virtual double fx(const DoubleVector &x); };
//struct RnFunctionAlpha2_22 : public RnFunction { virtual double fx(const DoubleVector &x); };

//struct RnFunctionAlpha3_11 : public RnFunction { virtual double fx(const DoubleVector &x); };
//struct RnFunctionAlpha3_12 : public RnFunction { virtual double fx(const DoubleVector &x); };
//struct RnFunctionAlpha3_21 : public RnFunction { virtual double fx(const DoubleVector &x); };
//struct RnFunctionAlpha3_22 : public RnFunction { virtual double fx(const DoubleVector &x); };

struct RnFunctionBetta11 : public RnFunction { virtual double fx(const DoubleVector &x); };
struct RnFunctionBetta12 : public RnFunction { virtual double fx(const DoubleVector &x); };
//struct RnFunctionBetta21 : public RnFunction { virtual double fx(const DoubleVector &x); };
//struct RnFunctionBetta22 : public RnFunction { virtual double fx(const DoubleVector &x); };

struct RnFunctionQamma1 : public RnFunction { virtual double fx(const DoubleVector &x); };
struct RnFunctionQamma2 : public RnFunction { virtual double fx(const DoubleVector &x); };

struct RnFunctionM : public RnFunction { virtual double fx(const DoubleVector &x); };

class SampleLoaderBorder
{
public:
    static void main();
};

#endif // SAMPLELOADERBORDER_H
