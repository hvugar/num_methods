#ifndef HYPERBOLICEQUATION_H
#define HYPERBOLICEQUATION_H

#include "global.h"
#include "doublevector.h"
#include <cmethods.h>
#include "tomasmethod.h"
#include "printer.h"

class MINIMUMSHARED_EXPORT IHyperbolicEquation
{
public:
    virtual double fi1(unsigned int i) const = 0;
    virtual double fi2(unsigned int i) const = 0;
    virtual double m1(unsigned int j) const = 0;
    virtual double m2(unsigned int j) const = 0;
    virtual double f(unsigned int i, unsigned int j) const = 0;

    virtual void calculateU(DoubleVector &u, double hx, double ht, unsigned int M, unsigned int N, double a=1.0, double lamda=0.25) const;
    virtual void calculateU(DoubleMatrix &u, double hx, double ht, unsigned int M, unsigned int N, double a=1.0, double lamda=0.25) const;
};

class MINIMUMSHARED_EXPORT IBackwardHyperbolicEquation
{
public:
    virtual double bfi1(unsigned int i) const = 0;
    virtual double bfi2(unsigned int i) const = 0;
    virtual double bm1(unsigned int j) const = 0;
    virtual double bm2(unsigned int j) const = 0;
    virtual double bf(unsigned int i, unsigned int j) const = 0;

    virtual void calculateU(DoubleVector &p, double hx, double ht, unsigned int M, unsigned int N, double a=1.0, double lamda=0.25) const;
    virtual void calculateU(DoubleMatrix &p, double hx, double ht, unsigned int M, unsigned int N, double a=1.0, double lamda=0.25) const;
};

class MINIMUMSHARED_EXPORT IHyperbolicEquation2D
{
public:
    virtual double fi1(unsigned int i, unsigned int j) const = 0;
    virtual double fi2(unsigned int i, unsigned int j) const = 0;
    virtual double m1(unsigned int j, unsigned int k) const = 0;
    virtual double m2(unsigned int j, unsigned int k) const = 0;
    virtual double m3(unsigned int i, unsigned int k) const = 0;
    virtual double m4(unsigned int i, unsigned int k) const = 0;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const = 0;

    virtual void calculateMVD(DoubleMatrix &u, double h1, double h2, double ht, double N1, double N2, double M, double a1=1.0, double a2=1.0) const;
    virtual void calculateMVD(DoubleCube &u, double h1, double h2, double ht, double N1, double N2, double M, double a1=1.0, double a2=1.0) const;
    //virtual void calculateMFS(DoubleCube &u, double h1, double h2, double ht, double N1, double N2, double M, double a1=1.0, double a2=1.0) const;
    //virtual void calculateMFS(DoubleMatrix &u, double h1, double h2, double ht, double N1, double N2, double M, double a1=1.0, double a2=1.0) const;
    virtual void calculateU1(DoubleCube &u, double h1, double h2, double ht, double N1, double N2, double M, double a1=1.0, double a2=1.0, double qamma=1.0) const;
    virtual void calculateU1(DoubleMatrix &u, double h1, double h2, double ht, double N1, double N2, double M, double a1=1.0, double a2=1.0, double qamma=1.0) const;
};

class MINIMUMSHARED_EXPORT IBackwardHyperbolicEquation2D
{
public:
    virtual double bfi1(unsigned int i, unsigned int j) const = 0;
    virtual double bfi2(unsigned int i, unsigned int j) const = 0;
    virtual double bm1(unsigned int j, unsigned int k) const = 0;
    virtual double bm2(unsigned int j, unsigned int k) const = 0;
    virtual double bm3(unsigned int i, unsigned int k) const = 0;
    virtual double bm4(unsigned int i, unsigned int k) const = 0;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const = 0;

    virtual void calculateU(DoubleMatrix &u, double h1, double h2, double ht, double N1, double N2, double M, double a1=1.0, double a2=1.0) const;
    virtual void calculateU(DoubleCube &p, double h1, double h2, double ht, double N1, double N2, double M, double a1=1.0, double a2=1.0) const;
    virtual void calculateU1(DoubleCube &p, double h1, double h2, double ht, double N1, double N2, double M, double a1=1.0, double a2=1.0, double qamma=1.0) const;
};


#endif // HYPERBOLICEQUATION_H
