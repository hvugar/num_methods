#ifndef HEATCONTROL2DELTA_H
#define HEATCONTROL2DELTA_H

#include <function.h>
#include <printer.h>
#include <projection.h>

class HeatControl2Delta : public RnFunction, Printer, Projection
{
public:
    HeatControl2Delta(unsigned int M, unsigned int N2, unsigned int N1);
    virtual ~HeatControl2Delta() {}

    //RnFunction
    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step);
    //Printer
    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;
    //Projection
    virtual void project(DoubleVector &x, int index);
    //
    void calculateU(const DoubleVector &x, DoubleMatrix &u);
    void calculateP(const DoubleVector &x, DoubleVector &g, const DoubleMatrix &u);
    void calculateGX(const DoubleVector &x, const DoubleMatrix& psi, DoubleVector& g, unsigned int k);
    void calculateGF(const DoubleVector &x, const DoubleMatrix& psi, DoubleVector& g, unsigned int k);
    void calculateG2(const DoubleVector &x, DoubleVector& g);

    double norm(const DoubleVector& x) const;

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

    void initialize();

    DoubleMatrix U;

    static void main();

    double gause1;
    double gause2;
    double gause3;
    //FILE* file;

private:
    double u(double x1, double x2, double t) { return x1*x1 + x2*x2 + t*t; }

    inline double fi(double x1, double x2) { return u(x1, x2, t0); }
    inline double m1(double x2, double t) { return u(x10, x2, t); }
    inline double m2(double x2, double t) { return u(x11, x2, t); }
    inline double m3(double x1, double t) { return u(x1, x20, t); }
    inline double m4(double x1, double t) { return u(x1, x21, t); }
    inline double f(unsigned int i, unsigned int j, unsigned int k);

    inline double pm1(double x2, double t) { return 0.0; }
    inline double pm2(double x2, double t) { return 0.0; }
    inline double pm3(double x1, double t) { return 0.0; }
    inline double pm4(double x1, double t) { return 0.0; }

    inline double g1(double t) const { return 100*t; }
    inline double g2(double t) const { return 200*t; }
    inline double g3(double t) const { return 300*t; }

    inline void psiDerivative(double &psiX1, double &psiX2, double e1, double e2, const DoubleMatrix &psi);
    void write(const char* fileName, const DoubleMatrix& m);

    const DoubleVector *px;
    DoubleVector optimal;
    DoubleVector initial;
};

#endif // HEATCONTROL2DELTA_H
