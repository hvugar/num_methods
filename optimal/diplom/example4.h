#ifndef BORDERPARABOLIC2D2_H
#define BORDERPARABOLIC2D2_H

#include <parabolicequation.h>
#include <doublevector.h>
#include <math.h>

class BorderParabolic2D2 : public IParabolicEquation2D
{
public:
    BorderParabolic2D2();
    virtual ~BorderParabolic2D2() {}

    virtual double fi(unsigned int i, unsigned int j) const;
    virtual double m1(unsigned int j, unsigned int k) const;
    virtual double m2(unsigned int j, unsigned int k) const;
    virtual double m3(unsigned int i, unsigned int k) const;
    virtual double m4(unsigned int i, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    static void main();

private:
    double x10;
    double x11;
    double x20;
    double x21;
    double t0;
    double t1;
    double h1;
    double h2;
    double ht;
    unsigned int N1;
    unsigned int N2;
    unsigned int M;
    unsigned int S;

    DoubleVector e;

    double v1(double t) const;
    double v2(double t) const;
    double v3(double t) const;
    double v4(double t) const;
};

#endif // BORDERPARABOLIC2D2_H
