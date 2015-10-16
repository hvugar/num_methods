#ifndef HEATCONTROL2DELTA_H
#define HEATCONTROL2DELTA_H

#include <function.h>
#include <doublevector.h>
#include <printer.h>
#include <projection.h>

class HeatControl2Delta : public RnFunction
{
public:
    HeatControl2Delta(unsigned int M, unsigned int N2, unsigned int N1);
    ~HeatControl2Delta();

    double fx(const DoubleVector& e);
    void gradient(const DoubleVector& e, DoubleVector& g, double gradient_step);

    void calculateU(const DoubleVector& e, DoubleMatrix& u);
    void calculateP(const DoubleVector& e, DoubleVector& g);
    void calculateG(const DoubleVector& e, const DoubleMatrix& psi, DoubleVector& g, unsigned int k);

    double u(double x1, double x2, double t) { return x1*x1 + x2*x2 + t*t; }

    DoubleMatrix U;
    DoubleMatrix uT;

    static void main();

private:
    double fi(double x1, double x2) { return u(x1, x2, t0); }
    double m1(double x2, double t) { return u(x10, x2, t); }
    double m2(double x2, double t) { return u(x11, x2, t); }
    double m3(double x1, double t) { return u(x1, x20, t); }
    double m4(double x1, double t) { return u(x1, x21, t); }

    double fxt(double x1, double x2, double t, const DoubleVector& e, unsigned int k);

    double pm1(double x2, double t) { return 0.0; }
    double pm2(double x2, double t) { return 0.0; }
    double pm3(double x1, double t) { return 0.0; }
    double pm4(double x1, double t) { return 0.0; }

    double f1(double t) { return t; }
    double f2(double t) { return t*t; }
    double f3(double t) { return t*t*t; }

    double t0;
    double t1;

    double x10;
    double x11;
    double x20;
    double x21;

    unsigned int M;
    unsigned int N1;
    unsigned int N2;
    unsigned int C;
    unsigned int L;

    double a1;
    double a2;

    double h1;
    double h2;
    double ht;

    void initialize();
};

struct HeatControl2DeltaPrinter : public Printer
{
    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &s, double a, RnFunction* f) const;
};

struct HeatControl2DeltaProjection : public Projection
{
    virtual void project(DoubleVector &x, int index);
};

#endif // HEATCONTROL2DELTA_H
