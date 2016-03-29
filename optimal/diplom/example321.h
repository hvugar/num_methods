#ifndef BORDERPARABOLIC2D321_H
#define BORDERPARABOLIC2D321_H

#include <parabolicequation.h>
#include <doublevector.h>
#include <math.h>

class BorderParabolic2D321 : public IParabolicEquation2D
{
public:
    BorderParabolic2D321();
    virtual ~BorderParabolic2D321() {}

    virtual double initial(unsigned int i, unsigned int j) const;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const;
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
};

#endif // BORDERPARABOLIC2D321_H
