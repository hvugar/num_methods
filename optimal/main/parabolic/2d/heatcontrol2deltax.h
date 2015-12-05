#ifndef HEATCONTROL2DELTAX_H
#define HEATCONTROL2DELTAX_H

#include <function.h>
#include <printer.h>
#include <projection.h>

class HeatControl2DeltaX : public RnFunction, Printer, Projection
{
public:
    HeatControl2DeltaX(unsigned int M, unsigned int N2, unsigned int N1);
    virtual ~HeatControl2DeltaX();

public:
    //RnFunction
    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step=0.000001);
    //Printer
    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;
    //Projection
    virtual void project(DoubleVector &x, int index);

    //
    inline void calculateU(const DoubleVector& e, DoubleMatrix& u);
    inline void calculateP(const DoubleVector& e, DoubleVector& g);
    inline void calculateGX(const DoubleVector& e, const DoubleMatrix& psi, DoubleVector& g, unsigned int k);
    inline void calculateGF(const DoubleVector &e, const DoubleMatrix& psi, DoubleVector& g, unsigned int k);
    inline void calculateG2(const DoubleVector& e, DoubleVector& g);
    inline void psiDerivative(double &psiX1, double& psiX2, double e1, double e2, const DoubleMatrix& psi);
    void initialize();
    void test();

    double calculateIntegral(const DoubleVector& x);
    double calculareNorm(const DoubleVector& x);

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

    DoubleMatrix U;
    DoubleMatrix uT;

    static void main();

private:
    double u(double x1, double x2, double t) { return x1*x1 + x2*x2 + t*t; }

    double fi(double x1, double x2) { return u(x1, x2, t0); }
    double m1(double x2, double t) { return u(x10, x2, t); }
    double m2(double x2, double t) { return u(x11, x2, t); }
    double m3(double x1, double t) { return u(x1, x20, t); }
    double m4(double x1, double t) { return u(x1, x21, t); }
    double fxt(unsigned int i, unsigned int j, unsigned k, const DoubleVector& f);

    double pm1(double x2, double t) { return 0.0; }
    double pm2(double x2, double t) { return 0.0; }
    double pm3(double x1, double t) { return 0.0; }
    double pm4(double x1, double t) { return 0.0; }

    double g1(double t) { return t; }
    double g2(double t) { return t*t; }
    double g3(double t) { return t*t*t; }

    void write(const char* fileName, const DoubleMatrix& m);
};

#endif // HEATCONTROL2DELTAX_H
