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

struct ControlFunction : public RnFunction, public IGradient, public Printer
{
    ControlFunction(double t0, double t1, double h);
    virtual ~ControlFunction();

protected:
    virtual double fx(const DoubleVector& u);
    virtual void gradient(const DoubleVector& u, DoubleVector &g);

    virtual void print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const;
    void print(const char* s, const DoubleVector& x) const;

public:
//    struct {
//        virtual double fx(double t, const DoubleVector &x, const DoubleVector &u) const {
//            //double t = x[0];
//            double x1 = x[1];
//            double x2 = x[2];
//            double u1 = x[3];
//            return (x1 - t*t*t)*(x1 - t*t*t)+(x2-t)*(x2-t)+(2*u1-t)*(2*u1-t);
//        }
//    } f0;

//    struct {
//        virtual double fx(double t, const DoubleVector &x, const DoubleVector &u) {
//            double x2 = x[2];
//            return 3.0 * x2 * x2;
//        }
//    } f1;

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
