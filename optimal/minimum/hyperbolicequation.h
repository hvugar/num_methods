#ifndef HYPERBOLICEQUATION_H
#define HYPERBOLICEQUATION_H

#include "global.h"
#include "doublevector.h"
#include <cmethods.h>
#include "tomasmethod.h"
#include "printer.h"

struct MINIMUMSHARED_EXPORT IHyperbolicEquation
{
public:
    virtual double fi1(unsigned int i) const = 0;
    virtual double fi2(unsigned int i) const = 0;
    virtual double m1(unsigned int j) const = 0;
    virtual double m2(unsigned int j) const = 0;
    virtual double f(unsigned int i, unsigned int j) const = 0;

    virtual void calculateU(DoubleMatrix &u, double hx, double ht, unsigned int M, unsigned int N, double a=1.0, double lamda=0.25) const;
    virtual void calculateU(DoubleVector &u, double hx, double ht, unsigned int M, unsigned int N, double a=1.0, double lamda=0.25) const;
};

struct MINIMUMSHARED_EXPORT IBackwardHyperbolicEquation
{
public:
    virtual double bfi1(unsigned int i) const = 0;
    virtual double bfi2(unsigned int i) const = 0;
    virtual double bm1(unsigned int j) const = 0;
    virtual double bm2(unsigned int j) const = 0;
    virtual double bf(unsigned int i, unsigned int j) const = 0;

    virtual void calculateU(DoubleMatrix &p, double hx, double ht, unsigned int M, unsigned int N, double a=1.0, double lamda=0.25) const;
};


#endif // HYPERBOLICEQUATION_H
