#ifndef CFUNCTION2_H
#define CFUNCTION2_H

#include <math.h>
#include <function.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

/**
 * @brief Optimal control problem
 * @param x1' = x2
 * @param x2' = -6*x1 + x2 + u - t + 1;
 * @param x1(0) = 0.0;
 * @param x2(0) = 0.0;
 * @param F(x(T)) = 0.0;
 * @param J = integral( (x1-t^2/2)^2 + (x2-t)^2 ) + F(x(T));
 *
 * u*(t) = t^2/2;
 * x1*(t) = t^2/2;
 * x2*(t) = t;
 */
struct CFunction2 : public RnFunction, public IGradient, public IPrinter
{
    CFunction2(double t0, double t1, double h);
    virtual ~CFunction2();

    virtual double fx(const DoubleVector& u);
    virtual void gradient(const DoubleVector& u, DoubleVector &g);

    virtual void print(unsigned int i, const DoubleVector& u, const DoubleVector &g, double alpha, RnFunction* f) const;

    double fx0(double t, const DoubleVector& x, double u) const;
    double F(double t, const DoubleVector& x, double u) const;
    double fx1(double t, const DoubleVector& x, double u) const;
    double fx2(double t, const DoubleVector& x, double u) const;
    double fp1(double t, const DoubleVector& x, const DoubleVector& psi, double u) const;
    double fp2(double t, const DoubleVector& x, const DoubleVector& psi, double u) const;
    double H(double t, const DoubleVector& x, double u, const DoubleVector& psi) const;
    void calculate_x(const DoubleVector& u);
    void calculate_psi(const DoubleVector& u);

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

#endif // CFUNCTION2_H
