#ifndef BORDERHYPERBOLIC2DS1_H
#define BORDERHYPERBOLIC2DS1_H

#include <hyperbolicequation.h>

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
};

#endif // BORDERHYPERBOLIC2DS1_H
