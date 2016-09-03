#ifndef CFUNCTION4_H
#define CFUNCTION4_H


/**
 * @brief The ControlFunction4 class
 * dx/dt = A(t)x(t) + B(t)u(t) + C(t), x(0) = x0;
 * J[u] = ||x(T)-Y|| --> min;
 * u(t) = Ks*x(ti);
 */

#include <math.h>
#include <function.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

class ControlFunction4 : public RnFunction, public IGradient, public IPrinter
{
public:
    ControlFunction4();
    virtual ~ControlFunction4();

    virtual double fx(const DoubleVector &u);
    virtual void gradient(const DoubleVector &u, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const;

    void calculate_X(const DoubleVector &u1, const DoubleVector &u2, DoubleVector &x1, DoubleVector &x2);
    void calculate_P(const DoubleVector &u1, const DoubleVector &u2, DoubleVector &x1, DoubleVector &x2);

    double fx1(double t, double x1, double x2, double u1, double u2);
    double fx2(double t, double x1, double x2, double u1, double u2);

    double px1(double t, double x1, double x2, double u1, double u2);
    double px2(double t, double x1, double x2, double u1, double u2);

    double H(double t, double x1, double x2, double u1, double u2, double p1, double p2);

    static void main(int argc, char* argv[]);

private:
    double t0;
    double t1;
    unsigned int N;
    double h;
};

#endif // CFUNCTION4_H
