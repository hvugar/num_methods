#ifndef HEATCONTROL2DELTA_H
#define HEATCONTROL2DELTA_H

#include <function.h>
#include <printer.h>
#include <projection.h>
#include <parabolicequation.h>

class HeatControl2Delta : public RnFunction, public IGradient, public IParabolicEquation2D, public IBackwardParabolicEquation2D, public IPrinter, public Projection
{
public:
    HeatControl2Delta(unsigned int M, unsigned int N2, unsigned int N1);
    virtual ~HeatControl2Delta() {}

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double alpha, RnFunction* fn) const;
    virtual void project(DoubleVector &x, int index);

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
    //unsigned int C;
    unsigned int L;

    double a1;
    double a2;

    double h1;
    double h2;
    double ht;

    double alpha;
    double gause_a;
    double gause_b;

    void initialize();

    DoubleMatrix U;
    static void main();

private:
    double u(double x1, double x2, double t) const { return x1*x1 + x2*x2 + t*t; }

    virtual double fi(unsigned int i, unsigned int j) const;
    virtual double m1(unsigned int j, unsigned int k) const;
    virtual double m2(unsigned int j, unsigned int k) const;
    virtual double m3(unsigned int i, unsigned int k) const;
    virtual double m4(unsigned int i, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    virtual double bfi(unsigned int i, unsigned int j) const;
    virtual double bm1(unsigned int j, unsigned int k) const;
    virtual double bm2(unsigned int j, unsigned int k) const;
    virtual double bm3(unsigned int i, unsigned int k) const;
    virtual double bm4(unsigned int i, unsigned int k) const;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const;

    inline double g1(double t) const { return 10*t; }
    inline double g2(double t) const { return 20*t; }
    inline double g3(double t) const { return 30*t; }

    inline void psiDerivative(double &psiX1, double &psiX2, double e1, double e2, const DoubleMatrix &psi);
    void write(const char* fileName, const DoubleMatrix& m);

    const DoubleVector *px;
    const DoubleMatrix *pu;
    DoubleVector O;
    DoubleVector I;
};

#endif // HEATCONTROL2DELTA_H
