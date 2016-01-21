#ifndef HEATCONTROL2DELTAF_H
#define HEATCONTROL2DELTAF_H

#include <function.h>
#include <doublevector.h>
#include <printer.h>
#include <projection.h>
#include <stdlib.h>

class HeatControl2DeltaF : public RnFunction, IPrinter, Projection
{
public:
    HeatControl2DeltaF(unsigned int M, unsigned int N2, unsigned int N1);
    virtual ~HeatControl2DeltaF();

    double fx(const DoubleVector& e);
    void gradient(const DoubleVector& e, DoubleVector& g, double gradient_step);

    void calculateU(const DoubleVector& f, DoubleMatrix& u);
    void calculateP(const DoubleVector& f, DoubleVector& g);
    void calculateG(const DoubleVector &f, DoubleVector &g, const std::vector<DoubleMatrix>& psi);

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

    double fxt(unsigned int i, unsigned int j, unsigned k, const DoubleVector &f);

    double pm1(double x2, double t) { return 0.0; }
    double pm2(double x2, double t) { return 0.0; }
    double pm3(double x1, double t) { return 0.0; }
    double pm4(double x1, double t) { return 0.0; }

    double f1(double t) { return t; }
    double f2(double t) { return t*t; }
    double f3(double t) { return t*t*t; }

    void psiDerivative(double &psiX1, double& psiX2, double e1, double e2, const DoubleMatrix& psi);

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

    double alpha;

    DoubleVector E;

    void initialize();
    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &s, double a, RnFunction* f) const;
    virtual void project(DoubleVector &x, int index);
};

#endif // HEATCONTROL2DELTAF_H
