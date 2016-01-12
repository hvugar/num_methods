#ifndef HEATCONTROLDELTAX_H
#define HEATCONTROLDELTAX_H

#include <function.h>
#include <printer.h>
#include <projection.h>

class HeatControlDeltaX : public RnFunction, public IGradient, public Printer, public Projection
{
public:
    HeatControlDeltaX(unsigned int M, unsigned int N, double a1);
    virtual ~HeatControlDeltaX() {}

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& f, DoubleVector& g);

    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &s, double a, RnFunction* f) const;
    virtual void project(DoubleVector &x, int index);

    inline double u(double x, double t) { return x*x + t*t; }
    inline double fxt(unsigned int i, unsigned int j, const DoubleVector &f);

    inline double fi(double x) { return u(x, 0.0); }
    inline double m1(double t) { return u(x0, t); }
    inline double m2(double t) { return u(x1, t); }

    inline double pfi(double x) { return 0.0; }
    inline double pm1(double t) { return 0.0; }
    inline double pm2(double t) { return 0.0; }

    double f1(double t) { return t*t; }
    double f2(double t) { return 0.5*t*t; }

    static void main();

    unsigned int N;
    unsigned int M;
    unsigned int C;

protected:
    void calculateU(const DoubleVector& f, DoubleVector& u);
    void calculateP(const DoubleVector& f, DoubleVector& g);
    void calculateG(const DoubleVector& f, const std::vector<DoubleVector>& psi, DoubleVector& g, unsigned int j);
    void initializeU();

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double a1;
    double hx;
    double ht;

    unsigned int L;

    DoubleVector U;
    DoubleVector uT;
};

#endif // HEATCONTROLDELTAX_H
