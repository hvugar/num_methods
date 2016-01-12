#ifndef CFUNCTION3_H
#define CFUNCTION3_H

#include <math.h>
#include <function.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

/**
 * @brief Optimal control problem
 * @param x1' = 3.0 * x2^2;
 * @param x2' = x1 + x2 - 2.0 * u - t^3 + 1.0;
 * @param x1(0) = 0.0;
 * @param x2(0) = 0.0;
 * @param T(x(T)) = (x2[T] - 1.0)^2;
 * @param J = integral( (x1-t^3)^2 + (x2-t)^2 + (2*u-t)^2 ) + T(x(T));
 *
 * u*(t) = t/2;
 * x1*(t) = t^3;
 * x2*(t) = t;
 */
struct CFunction3 : public RnFunction, public IGradient, public Printer
{
public:
    CFunction3(double t0, double t1, double h);
    virtual ~CFunction3();

    virtual double fx(const DoubleVector& u);

    virtual void gradient(const DoubleVector& u, DoubleVector &g);

    virtual void print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const;

private:
    double fx0(double t, const DoubleVector& x, const DoubleVector& u) const;
    double T(double t, const DoubleVector& x) const;
    double fx1(double t, const DoubleVector& x, const DoubleVector& u) const;
    double fx2(double t, const DoubleVector& x, const DoubleVector& u) const;
    double fp1(double t, const DoubleVector& x, const DoubleVector& psi, const DoubleVector& u);
    double fp2(double t, const DoubleVector& x, const DoubleVector& psi, const DoubleVector& u);
    double H(double t, const DoubleVector& x, const DoubleVector& u, const DoubleVector& psi);
    void calculate_x(const DoubleVector& u);
    void calculate_psi(const DoubleVector& u);

public:
    double t0;
    double t1;
    double h;
    int n;

    DoubleVector t;
    DoubleVector x1;
    DoubleVector x2;
    DoubleVector psi1;
    DoubleVector psi2;

public:
    static void main();
};

#endif // CFUNCTION3_H
