#ifndef BORDERHYPERBOLIC2DS1_H
#define BORDERHYPERBOLIC2DS1_H

#include <hyperbolicequation.h>
#include <math.h>

class BorderHyperbolic2DS1 : public IHyperbolicEquation2D
{
public:
    BorderHyperbolic2DS1();
    virtual ~BorderHyperbolic2DS1() {}

    virtual double fi1(unsigned int i, unsigned int j) const;
    virtual double fi2(unsigned int i, unsigned int j) const;
    virtual double m1(unsigned int j, unsigned int k) const;
    virtual double m2(unsigned int j, unsigned int k) const;
    virtual double m3(unsigned int i, unsigned int k) const;
    virtual double m4(unsigned int i, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const;
    double u(unsigned int i, unsigned int j, unsigned int k) const;

private:
    double h1;
    double h2;
    double ht;
    unsigned int N1;
    unsigned int N2;
    unsigned int M;

public:
    static void main();
};

#endif // BORDERHYPERBOLIC2DS1_H
