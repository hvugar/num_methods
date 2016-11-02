#ifndef PARABOLICQUATION_H
#define PARABOLICQUATION_H

#include "global.h"
#include "vector2d.h"
#include "matrix2d.h"
#include "matrix3d.h"
#include "printer.h"

class MINIMUMSHARED_EXPORT IParabolicEquation
{
public:
    virtual double initial(unsigned int i) const = 0;
    virtual double boundary(Boundary type, unsigned int j) const = 0;
    virtual double f(unsigned int i, unsigned int j) const = 0;
    virtual void calculateU(DoubleVector &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;
    virtual void calculateU(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;
    virtual void calculateN(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;

    virtual void calculateL(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;
    virtual void calculateN1(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;

    virtual void calculateU1(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;
};

class MINIMUMSHARED_EXPORT IBackwardParabolicEquation
{
public:
    virtual double binitial(unsigned int i) const = 0;
    virtual double bboundary(Boundary type, unsigned int j) const = 0;
    virtual double bf(unsigned int i, unsigned int j) const = 0;
    virtual void calculateU(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;
};

class MINIMUMSHARED_EXPORT IParabolicEquation2D
{
public:
    virtual double initial(unsigned int i, unsigned int j) const = 0;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const = 0;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const = 0;

    // Метод дробных шагов
    void caluclateMVD(DoubleMatrix &u, double hx1, double hx2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;
    void caluclateMVD(DoubleCube &u, double hx1, double hx2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;
    void caluclateMVD1(DoubleMatrix &u, double hx1, double hx2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;
    void caluclateMVD1(DoubleCube &u, double hx1, double hx2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;

    // Метод дробных шагов
    void caluclateMFS(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;
    void caluclateMFS(DoubleCube &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;
};

class MINIMUMSHARED_EXPORT IBackwardParabolicEquation2D
{
public:
    virtual double binitial(unsigned int i, unsigned int j) const = 0;
    virtual double bboundary(unsigned int i, unsigned int j, unsigned int k) const = 0;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const = 0;
    void caluclateMVD(DoubleCube &psi, double hx1, double hx2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1=1.0, double a2=1.0) const;
};

#endif // PARABOLICQUATION_H
