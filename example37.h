#ifndef PARABOLICCONTROL2D333_H
#define PARABOLICCONTROL2D333_H

#include "parabolicequation.h"
#include "printer.h"
#include "gradient_cjt.h"
#include <math.h>

class Parabolic2DControl37 : public RnFunction, public IGradient, public IParabolicEquation2D, public IBackwardParabolicEquation2D, public IPrinter, public IProjection
{
public:
    Parabolic2DControl37(unsigned int m, unsigned int n2, unsigned int n1);
    virtual ~Parabolic2DControl37() {}

    virtual double fx(const DoubleVector &v);
    virtual void gradient(const DoubleVector &v, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &v, const DoubleVector &g, double alpha, RnFunction* fn) const;
    virtual void project(DoubleVector &v, int index);

    virtual double initial(unsigned int i, unsigned int j) const;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    virtual double binitial(unsigned int i, unsigned int j) const;
    virtual double bboundary(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double bf(unsigned int i, unsigned int j, unsigned int k) const;

    void calculateGF(const DoubleVector &v, const DoubleMatrix& psi, DoubleVector& g, unsigned int k);

    static void main(int argc, char *argv[]);
private:
    double u(double x1, double x2, double t) const { return x1*x1 + x2*x2 + t*t; }
    double norm(const DoubleVector& v) const;

    inline double v1(double t) const { return 22*t; }
    inline double v2(double t) const { return 18*t; }
    inline double v3(double t) const { return 25*t; }

    double t0;
    double t1;

    double x10;
    double x11;
    double x20;
    double x21;

    unsigned int M;
    unsigned int N1;
    unsigned int N2;
    unsigned int L;

    double a1;
    double a2;

    double h1;
    double h2;
    double ht;

    double alpha;
    double gause_a;
    double gause_b;

    DoubleMatrix U;
    const DoubleVector *pv;
    const DoubleMatrix *pu;
    DoubleVector E;
};

#endif // PARABOLICCONTROL2D333_H
