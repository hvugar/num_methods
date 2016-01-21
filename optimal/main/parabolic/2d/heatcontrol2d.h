#ifndef HEATCONTROL2D_H
#define HEATCONTROL2D_H

#include <function.h>
#include <doublevector.h>
#include <printer.h>
#include <parabolicequation.h>

class HeatControl2D :public RnFunction, public IGradient, public IPrinter, public IParabolicEquation2D
{
public:
    HeatControl2D(unsigned int M, unsigned int N2, unsigned int N1);
    virtual ~HeatControl2D();

    unsigned int N1;
    unsigned int N2;
    unsigned int M;
    unsigned int C;
    double h1;
    double h2;
    double ht;

    double t0;
    double t1;
    double x10;
    double x11;
    double x20;
    double x21;

    double a1;
    double a2;

    DoubleMatrix U;
    //DoubleMatrix uT;

    //void calculateU(const DoubleVector& f);
    //void calculateU1(const DoubleVector &f);
    void calculateP(const DoubleVector &f, const DoubleMatrix &u, DoubleVector& g);

    static void main();

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g);
    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &s, double a, RnFunction* f) const;

private:
    double u1(double x1, double x2, double t) const;

    double fi(unsigned int i, unsigned int j) const;
    double m1(unsigned int j, double k) const;
    double m2(unsigned int j, double k) const;
    double m3(unsigned int i, double k) const;
    double m4(unsigned int i, double k) const;
    double f(unsigned int i, unsigned j, double k) const;

    double pm1(double x2, double t) { return 0.0; }
    double pm2(double x2, double t) { return 0.0; }
    double pm3(double x1, double t) { return 0.0; }
    double pm4(double x1, double t) { return 0.0; }

    double fxt(double x1, double x2, double t);


    const DoubleVector *pf;
};

#endif // HEATCONTROL2D_H
