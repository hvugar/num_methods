#ifndef HYPERBOLICEQUATION_H
#define HYPERBOLICEQUATION_H

#include "global.h"
#include "doublevector.h"

class MINIMUMSHARED_EXPORT HyperbolicEquation
{
public:
    HyperbolicEquation(double t0 = 0.0, double t1 = 1.0, double x0 = 0.0, double x1 = 1.0, double a = 1.0, unsigned int M = 100, unsigned int N = 100);
    virtual ~HyperbolicEquation();

    void setTimeInterval(double t0, double t1);
    void setLengthInterval(double x0, double x1);
    void setPartNumber(unsigned int M, unsigned int N);

    virtual double fi1(double x) const = 0;
    virtual double fi2(double x) const = 0;
    virtual double m1(double t) const = 0;
    virtual double m2(double t) const = 0;
    virtual double f(double x, double t) const = 0;

    virtual void calculate(DoubleVector& u);

protected:
    double t0;
    double t1;
    double x0;
    double x1;
    double a;
    /* number of parts of time */
    unsigned int M;
    /* number of parts of length */
    unsigned int N;
    double ht;
    double hx;
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
