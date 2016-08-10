#ifndef CAUCHYPROBLEMTEST_H
#define CAUCHYPROBLEMTEST_H

#include "headers.h"

struct A : public CauchyProblem
{
    double f(double x, const DoubleVector &y) const
    {
        double y1 = y[0];
        double y2 = y[1];
        double y3 = y[2];
        return x*y1 + y2 + 2.0*y3 + (2.0*x - 2.0*x*x*x - 2.0*sin(x));
    }
};

struct B : public CauchyProblem
{
    double f(double x, const DoubleVector &y) const
    {
        double y1 = y[0];
        double y2 = y[1];
        double y3 = y[2];
        return y1 + 2.0*x*x;
    }
};

struct C : public CauchyProblem
{
    double f(double x, const DoubleVector &y) const
    {
        double y1 = y[0];
        double y2 = y[1];
        double y3 = y[2];
        return y1 + 2.0*y2 + y3 + (cos(x) - x*x - 2.0*x*x*x - sin(x));
    }
};

class P : public IParabolicEquation
{
public:
    double initial(unsigned int i) const
    {
        double x = i*hx;
        return x*x*x;
    }
    double boundary(Boundary type, unsigned int j) const
    {
        double t = j*ht;
        if (type == Left) return t*t;
        if (type == Right) return t*t+1.0;
        return 0.0;
    }
    double f(unsigned int i, unsigned int j) const
    {
        double x = i*hx;
        double t = j*ht;
        return 2.0*t - 6.0*x*a*a;
    }
    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
    double a;
};

class P1 : public IParabolicEquation
{
public:
    double initial(unsigned int i) const
    {
        double x = i*hx;
        return x*x*x*x;
    }
    double boundary(Boundary type, unsigned int j) const
    {
        double t = j*ht;
        if (type == Left) return 0.0;
        if (type == Right) return 4.0;
        return 0.0;
    }
    double f(unsigned int i, unsigned int j) const
    {
        double x = i*hx;
        double t = j*ht;
        return 4.0*t*t*t - 12.0*a*a*x*x;
    }
    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
    double a;
};

class CauchyProblemTest
{
public:
    static void main();
};

#endif // CAUCHYPROBLEMTEST_H
