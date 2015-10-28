#ifndef HEATEQUATION_H
#define HEATEQUATION_H

#include "global.h"
#include "tomasmethod.h"
#include "printer.h"

struct MINIMUMSHARED_EXPORT ParabolicEquation
{
    ParabolicEquation(double t0 = 0.0, double t1 = 1.0, double x0 = 0.0, double x1 = 1.0, double a = 1.0, unsigned int M = 100, unsigned int N = 100);
    virtual ~ParabolicEquation();

    void setTimeInterval(double t0, double t1);
    void setLengthInterval(double x0, double x1);
    void setPartNumber(unsigned int M, unsigned int N);

    virtual double fi(unsigned int i, double x) const = 0;
    virtual double m1(unsigned int j, double t) const = 0;
    virtual double m2(unsigned int j, double t) const = 0;
    virtual double f(unsigned int i, double x, unsigned int j, double t) const = 0;

    virtual void calculateU(DoubleVector& u);

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

struct MINIMUMSHARED_EXPORT ConjuctionParabolicEquation
{
    ConjuctionParabolicEquation(double t0 = 0.0, double t1 = 1.0, double x0 = 0.0, double x1 = 1.0, double a = -1.0, unsigned int M = 100, unsigned int N = 100);
    virtual ~ConjuctionParabolicEquation();

    void setTimeInterval(double t0, double t1);
    void setLengthInterval(double x0, double x1);
    void setPartNumber(unsigned int M, unsigned int N);

    virtual double fi(double x) const = 0;
    virtual double m1(double t) const = 0;
    virtual double m2(double t) const = 0;
    virtual double f(double x, double t) const = 0;

    virtual void calculateU(DoubleVector& u);

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

struct MINIMUMSHARED_EXPORT ParabolicEquation2D
{
    ParabolicEquation2D(double t0 = 0.0, double t1 = 1.0, double x10 = 0.0, double x11 = 1.0, double x20 = 0.0, double x21 = 1.0, double a1 = 1.0, double a2 = 1.0, unsigned int M = 100, unsigned int N1 = 100, unsigned int N2 = 100);
    virtual double fi(double x1, double x2) const = 0;
    virtual double m1(double x2, double t) const = 0;
    virtual double m2(double x2, double t) const = 0;
    virtual double m3(double x1, double t) const = 0;
    virtual double m4(double x1, double t) const = 0;
    virtual double f(double x1, double x2, double t) const = 0;

    void setBorders(double t0, double t1, double x10, double x11, double x20, double x21);
    void setPartNumbers(unsigned int N1, unsigned int N2, unsigned M);

    void calculate(DoubleMatrix& u);
    void calculateBack(DoubleMatrix& u);

private:
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

#endif // HEATEQUATION_H
