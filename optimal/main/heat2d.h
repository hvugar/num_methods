#ifndef HEAT2D_H
#define HEAT2D_H

#include <function.h>
#include <printer.h>
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

    double psi_fi(int i, int j);
    double psi_m1(double x2, double t);
    double psi_m2(double x2, double t);
    double psi_m3(double x1, double t);
    double psi_m4(double x1, double t);

    void calculateU(const DoubleVector& f);
    void calculateP(const DoubleVector& f, DoubleVector& g);

    static void main();

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

    int N1;
    int N2;
    int M;
    int C;

    double a1;
    double a2;

    DoubleMatrix mu;
};

struct Heat2DControlPrinter : public Printer
{
    virtual void print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const;
};

#endif // HEAT2D_H
