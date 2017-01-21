#ifndef HYPERBOLICEQUATION_H
#define HYPERBOLICEQUATION_H

#include "iibvp.h"

class MINIMUMSHARED_EXPORT IHyperbolicEquation : public IPDEInitialBoundaryValueProblem
{
public:
    virtual double f(unsigned int i, unsigned int j) const = 0;

    virtual double initial1(unsigned int i) const = 0;
    virtual double initial2(unsigned int i) const = 0;
    virtual double boundary(Boundary type, unsigned int j) const = 0;

    virtual void calculateU(DoubleVector &u, double hx, double ht, unsigned int M, unsigned int N, double a=1.0, double lamda=0.25) const;
    virtual void calculateU(DoubleMatrix &u, double hx, double ht, unsigned int M, unsigned int N, double a=1.0, double lamda=0.25) const;
};

class MINIMUMSHARED_EXPORT IBackwardHyperbolicEquation : public IPDEInitialBoundaryValueProblem
{
public:
    virtual double binitial1(unsigned int i) const = 0;
    virtual double binitial2(unsigned int i) const = 0;
    virtual double bboundary(Boundary type, unsigned int j) const = 0;
    virtual double bf(unsigned int i, unsigned int j) const = 0;
    virtual void calculateU(DoubleVector &p, double hx, double ht, unsigned int M, unsigned int N, double a=1.0, double lamda=0.25) const;
    virtual void calculateU(DoubleMatrix &p, double hx, double ht, unsigned int M, unsigned int N, double a=1.0, double lamda=0.25) const;
};

class MINIMUMSHARED_EXPORT IHyperbolicEquation2D  : public IPDEInitialBoundaryValueProblem
{
public:
    virtual double initial1(unsigned int i, unsigned int j) const = 0;
    virtual double initial2(unsigned int i, unsigned int j) const = 0;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const = 0;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const = 0;
    virtual void calculateMVD(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;
    virtual void calculateMVD(DoubleCube &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;
    virtual void calculateU1(DoubleCube &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0, double qamma=1.0) const;
    virtual void calculateU1(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0, double qamma=1.0) const;
};

class MINIMUMSHARED_EXPORT IBackwardHyperbolicEquation2D  : public IPDEInitialBoundaryValueProblem
{
public:
    virtual double binitial1(unsigned int i, unsigned int j) const = 0;
    virtual double binitial2(unsigned int i, unsigned int j) const = 0;
    virtual double bboundary(unsigned int i, unsigned int j, unsigned int k) const = 0;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const = 0;
    virtual void calculateU(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;
    virtual void calculateU(DoubleCube &p, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;
    virtual void calculateU1(DoubleCube &p, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0, double qamma=1.0) const;
};


#endif // HYPERBOLICEQUATION_H
