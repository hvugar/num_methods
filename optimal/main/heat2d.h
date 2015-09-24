#ifndef HEAT2D_H
#define HEAT2D_H

#include <function.h>
#include <doublevector.h>

struct Heat2DControl : public RnFunction
{
public:
    Heat2DControl();
    ~Heat2DControl();

    virtual double fx(const DoubleVector& x);
    virtual void gradient(double step, const DoubleVector& x, DoubleVector& g);

    double u(double x1, double x2, double t);
    double U(double x1, double x2);

    double f(double x1, double x2, double t);
    double fi(double x1, double x2);
    double m1(double x2, double t);
    double m2(double x2, double t);
    double m3(double x1, double t);
    double m4(double x1, double t);

    double psi_fi(double x1, double x2);
    double psi_m1(double x2, double t);
    double psi_m2(double x2, double t);
    double psi_m3(double x1, double t);
    double psi_m4(double x1, double t);

    void calculateX(const DoubleMatrix& u, const DoubleMatrix &f, DoubleMatrix& x1, DoubleMatrix& x2);
    void calculatePsi(const DoubleMatrix& u, DoubleMatrix& x1, DoubleMatrix& x2);

private:
    double t0;
    double t1;
    double x11;
    double x12;
    double x21;
    double x22;

    double dt;
    double h1;
    double h2;

    unsigned int N1;
    unsigned int N2;
    unsigned int M;

    double a1;
    double a2;

    DoubleMatrix mu;
    DoubleMatrix u0;
    DoubleMatrix mf;
    DoubleMatrix mp;
    DoubleMatrix mx1;
    DoubleMatrix mx2;
};

#endif // HEAT2D_H
