#ifndef HEATCONTROL2D_H
#define HEATCONTROL2D_H

#include <function.h>
#include <doublevector.h>
#include <printer.h>

class HeatControl2D :public RnFunction
{
public:
    HeatControl2D();
    ~HeatControl2D();

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
    DoubleCube f;
    DoubleCube g;

    void calculateU(const DoubleVector& f);
    void calculateU1(const DoubleVector& f);
    void calculateP(const DoubleVector& f, DoubleCube& g);
    void calculateP1(const DoubleVector& f);
    void calculateP2(const DoubleVector& f);

    static void main();

protected:
    virtual double fx(const DoubleVector& x);
    virtual void gradient(double step, const DoubleVector& x, DoubleVector& g);

private:
    double u(double x1, double x2, double t);

    double fi(double x1, double x2);
    double m1(double x2, double t);
    double m2(double x2, double t);
    double m3(double x1, double t);
    double m4(double x1, double t);

    double fxt(double x1, double x2, double t);
};

#endif // HEATCONTROL2D_H
