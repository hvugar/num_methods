#ifndef HYPERBOLICEQUATION_H
#define HYPERBOLICEQUATION_H

#include "global.h"
#include "doublevector.h"

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
//protected:
//    virtual void calculateU(DoubleVector& u, unsigned int M = 1000, unsigned int N = 1000, double t0 = 0.0, double t1 = 1.0, double x0 = 0.0, double x1 = 1.0, double a = 1.0, double lamda=0.25) const;
//    virtual void calculateU(DoubleMatrix& u, unsigned int M = 1000, unsigned int N = 1000, double t0 = 0.0, double t1 = 1.0, double x0 = 0.0, double x1 = 1.0, double a = 1.0, double lamda=0.25) const;
//    virtual void calculateP(DoubleVector& u, unsigned int M = 1000, unsigned int N = 1000, double t0 = 0.0, double t1 = 1.0, double x0 = 0.0, double x1 = 1.0, double a = 1.0, double lamda=0.25) const;
};

struct MINIMUMSHARED_EXPORT HyperbolicEquation
{
public:
    HyperbolicEquation(unsigned int M = 100, unsigned int N = 100, double t0 = 0.0, double t1 = 1.0, double x0 = 0.0, double x1 = 1.0, double a = 1.0);
    virtual ~HyperbolicEquation();

    void setTimeInterval(double t0, double t1);
    void setLengthInterval(double x0, double x1);
    void setPartNumber(unsigned int M, unsigned int N);

    virtual double fi1(unsigned int i) const = 0;
    virtual double fi2(unsigned int i) const = 0;
    virtual double m1(unsigned int j) const = 0;
    virtual double m2(unsigned int j) const = 0;
    virtual double f(unsigned int i, unsigned int j) const = 0;

    virtual void calculate(DoubleVector& u);

protected:
    /* number of parts of time */
    unsigned int M;
    /* number of parts of length */
    unsigned int N;
    double t0;
    double t1;
    double x0;
    double x1;
    double a;
    double ht;
    double hx;
    double lamda;
};

struct MINIMUMSHARED_EXPORT HyperbolicEquation2D
{
    HyperbolicEquation2D(double t0 = 0.0, double t1 = 1.0, double x10 = 0.0, double x11 = 1.0, double x20 = 0.0, double x21 = 1.0, double a1 = 1.0, double a2 = 1.0, unsigned int M = 100, unsigned int N1 = 100, unsigned int N2 = 100);
    virtual ~HyperbolicEquation2D();

    virtual double fi1(unsigned int i, unsigned int j) const = 0;
    virtual double fi2(unsigned int i, unsigned int j) const = 0;
    virtual double m1(unsigned int j, unsigned int k) const = 0;
    virtual double m2(unsigned int j, unsigned int k) const = 0;
    virtual double m3(unsigned int i, unsigned int k) const = 0;
    virtual double m4(unsigned int i, unsigned int k) const = 0;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const = 0;

    void calculateImplicitly(DoubleMatrix& m);
    void calculateExplicitly(DoubleMatrix& m);

protected:
    double t0;
    double t1;
    double x10;
    double x11;
    double x20;
    double x21;
    double a1;
    double a2;
    unsigned int M;
    unsigned int N1;
    unsigned int N2;
    double ht;
    double h1;
    double h2;
};

#endif // HYPERBOLICEQUATION_H
