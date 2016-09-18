#ifndef EXAMPLE1_H
#define EXAMPLE1_H

#include <global.h>
#include <math.h>
#include <matrix2d.h>
#include <vector2d.h>
#include <rungekutta.h>
#include <printer.h>

class Example1
{
public:
    Example1();
    virtual ~Example1() {}

    void calculate();

    double A(unsigned int i, unsigned int j, double t) const;
    double B(unsigned int i, double t) const;
    double S0(const DoubleVector &a, double c, double t) const;

    unsigned int N = 4;
};

#endif // EXAMPLE1_H
