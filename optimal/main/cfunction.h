#ifndef CFUNCTION_H
#define CFUNCTION_H

#include <math.h>
#include <function.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

struct CFunction
{
    virtual double fx(double t, const DoubleVector &x, double u) { return 0.0; }
    virtual double fx(double t, const DoubleVector &x, const DoubleVector &psi, double u)  { return 0.0; }
    virtual double fx(double t, const DoubleVector &x)  { return 0.0; }
};

struct ControlFunction : public RnFunction, public Printer
{
    ControlFunction(double t0, double t1, double h);
    virtual ~ControlFunction();

protected:
    virtual double fx(const DoubleVector& u);
    virtual void gradient(const DoubleVector& u, DoubleVector &g, double gradient_step);

    virtual void print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const;
    void print(const char* s, const DoubleVector& x) const;

public:
    CFunction *fx0;
    CFunction *T;
    CFunction *fx1;
    CFunction *fx2;
//    CFunction *fp1;
//    CFunction *fp2;

    CFunction **dx;

private:
    double H(double t, const DoubleVector &x, double u, const DoubleVector &psi);
    double Integral(const DoubleVector &u);
    void calculate_x(const DoubleVector& u, DoubleVector& x1, DoubleVector& x2);
    void calculate_psi(const DoubleVector& u, DoubleVector& psi1, DoubleVector& psi2, DoubleVector& x1, DoubleVector& x2);

    double t0;
    double t1;
    double h;
    int n;

    DoubleVector t;
//    DoubleVector x1;
//    DoubleVector x2;
//    DoubleVector psi1;
//    DoubleVector psi2;

public:
    static void main();
};

#endif // CFUNCTION_H
