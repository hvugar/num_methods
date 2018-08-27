#ifndef TEMPLATEA_H
#define TEMPLATEA_H

#include <vector2d.h>
#include <cmethods.h>
#include <printer.h>

class TemplateA
{
public:
    static void Main(int argc, char** argv);

    double r(unsigned int i, double x) const;
    double p(unsigned int i, double x) const;
    double q(unsigned int i, double x) const;
    double f(unsigned int i, double x) const;

    void calculate(double h, unsigned int N, DoubleVector &y, double a, double b) const;
    double fx(double x) const;
};

#endif // TEMPLATEA_H
