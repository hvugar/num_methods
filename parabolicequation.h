#ifndef PARABOLICEQUATION_H
#define PARABOLICEQUATION_H

#include "doublevector.h"

enum Boundary
{
    Left = 0,
    Right = 1
};

void tridioganalMatrixAlgorithm(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int n);

class IParabolicEquation
{
public:
    virtual double initial(unsigned int i) const = 0;
    virtual double boundary(Boundary type, unsigned int j) const = 0;
    virtual double f(unsigned int i, unsigned int j) const = 0;
    virtual void calculateU(DoubleVector &u, double hx, double ht, unsigned int N, unsigned int M, double a) const;
	virtual void calculateU(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;
};

class IParabolicEquation2D
{
public:
    virtual double initial(unsigned int i, unsigned int j) const = 0;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const = 0;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const = 0;
    void caluclateMVD(DoubleMatrix &u, double hx1, double hx2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const;
};

class IBackwardParabolicEquation
{
public:
    virtual double binitial(unsigned int i) const = 0;
    virtual double bboundary(Boundary type, unsigned int j) const = 0;
    virtual double bf(unsigned int i, unsigned int j) const = 0;
    virtual void calculateU(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const;
};

class IBackwardParabolicEquation2D
{
public:
    virtual double binitial(unsigned int i, unsigned int j) const = 0;
    virtual double bboundary(unsigned int i, unsigned int j, unsigned int k) const = 0;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const = 0;
    void caluclateMVD(DoubleCube &psi, double hx1, double hx2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const;
};

#endif // PARABOLICEQUATION_H
