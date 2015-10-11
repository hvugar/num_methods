#ifndef HEADCONTROL2D1_H
#define HEADCONTROL2D1_H

#include <function.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <printer.h>

/**
 * @brief The HeadControl2D class
 * du/dt = a1 * d^2u/dx1^2 + a2 * d^2u/dx2^2 + SUM(fi(t)*delta(x1-E1i)*delta(x2-E2i)
 */

class HeadControl2D1 : public RnFunction
{
public:
    HeadControl2D1();

    virtual double fx(const DoubleVector& x);
    virtual void gradient(double step, const DoubleVector& x, DoubleVector& g);

    void calculateU(const DoubleVector& E, DoubleMatrix& u);
    void calculateP(const DoubleVector& E, DoubleVector& g);

    void initialize();

    static void main();

    unsigned int C;

protected:
    double u(double x1, double x2, double t);

    double fi(double x1, double x2);
    double m1(double x2, double t);
    double m2(double x2, double t);
    double m3(double x1, double t);
    double m4(double x1, double t);

    double f(double x1, double x2, double t);

    double psi_fi(int i, int j);
    double psi_m1(double x2, double t);
    double psi_m2(double x2, double t);
    double psi_m3(double x1, double t);
    double psi_m4(double x1, double t);

private:
    double t0;
    double t1;
    double x10;
    double x11;
    double x20;
    double x21;

    double a1;
    double a2;

    unsigned int N1;
    unsigned int N2;
    unsigned int M;

    double h1;
    double h2;
    double ht;

    DoubleMatrix U;
    DoubleMatrix mu;
};

struct HeadControl2D1Printer : public Printer
{
    virtual void print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const;
};

#endif // HEADCONTROL2D1_H
