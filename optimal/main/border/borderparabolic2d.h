#ifndef BORDERPARABOLIC2D_H
#define BORDERPARABOLIC2D_H

#include <parabolicequation.h>
#include <math.h>

class BorderParabolic2D : public IParabolicEquation2D
{
public:
    BorderParabolic2D();
    virtual ~BorderParabolic2D();

    virtual double fi(unsigned int i, unsigned int j) const;
    virtual double m1(unsigned int j, unsigned int k) const;
    virtual double m2(unsigned int j, unsigned int k) const;
    virtual double m3(unsigned int i, unsigned int k) const;
    virtual double m4(unsigned int i, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;

    double u(unsigned int i, unsigned int j, unsigned int k) const;

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
    double a1;
    double a2;
};

#endif // BORDERPARABOLIC2D_H
