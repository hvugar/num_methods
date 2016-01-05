#ifndef HEATEQUATION_H
#define HEATEQUATION_H

#include "global.h"
#include "tomasmethod.h"
#include "printer.h"

struct MINIMUMSHARED_EXPORT IParabolicEquation
{
    virtual double fi(unsigned int i) const = 0;
    virtual double m1(unsigned int j) const = 0;
    virtual double m2(unsigned int j) const = 0;
    virtual double f(unsigned int i, unsigned int j) const = 0;

    virtual void calculateU(DoubleVector &u, double hx, double ht, unsigned int N, unsigned int M, double a) const;
    virtual void calculateU(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const;
};

struct MINIMUMSHARED_EXPORT ParabolicEquation
{
    ParabolicEquation(double t0 = 0.0, double t1 = 1.0, double x0 = 0.0, double x1 = 1.0, unsigned int M = 100, unsigned int N = 100, double a = 1.0);
    virtual ~ParabolicEquation();

    void setTimeInterval(double t0, double t1);
    void setLengthInterval(double x0, double x1);
    void setPartNumber(unsigned int M, unsigned int N);

    virtual double fi(unsigned int i) const = 0;
    virtual double m1(unsigned int j) const = 0;
    virtual double m2(unsigned int j) const = 0;
    virtual double f(unsigned int i, unsigned int j) const = 0;

    virtual void calculateU(DoubleVector& u);

protected:
    /* initial time */
    double t0;
    /* end time */
    double t1;
    /* start point of length */
    double x0;
    /* end point of length */
    double x1;
    /* number of parts of time */
    unsigned int M;
    /* number of parts of length */
    unsigned int N;
    double ht;
    double hx;
    double a;
};

struct MINIMUMSHARED_EXPORT ParabolicEquation2D
{
    ParabolicEquation2D(double t0 = 0.0, double t1 = 1.0, double x10 = 0.0, double x11 = 1.0, double x20 = 0.0, double x21 = 1.0, unsigned int M = 100, unsigned int N1 = 100, unsigned int N2 = 100, double a1 = 1.0, double a2 = 1.0);
    virtual double fi(unsigned int i, unsigned int j) const = 0;
    virtual double m1(unsigned int j, unsigned int k) const = 0;
    virtual double m2(unsigned int j, unsigned int k) const = 0;
    virtual double m3(unsigned int i, unsigned int k) const = 0;
    virtual double m4(unsigned int i, unsigned int k) const = 0;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const = 0;

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
    unsigned int M;
    unsigned int N1;
    unsigned int N2;
    double a1;
    double a2;
    double ht;
    double h1;
    double h2;
};

#endif // HEATEQUATION_H
