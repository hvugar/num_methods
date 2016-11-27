#ifndef EXAMPLE3_H
#define EXAMPLE3_H

#include <vector2d.h>
#include <printer.h>
#include <math.h>

class Example3
{
public:
    Example3();

    double a(unsigned int i, unsigned int j, unsigned int k) const;
    double b(unsigned int i, unsigned int k) const;

    double X1(unsigned int i) const;
    double X2(unsigned int i) const;

    void calculate();
    void calculate1();

    double h = 0.01;
    double N = 100;
};

#endif // EXAMPLE3_H
