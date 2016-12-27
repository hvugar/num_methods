#ifndef EXAMPLE5_H
#define EXAMPLE5_H

#include <vector2d.h>
#include <printer.h>
#include <math.h>

//#define SAMPLE_1
#define SAMPLE_2

class Example5
{
public:
    Example5();
    void calculate_n2();
    void calculate_n4();
    void calculate_n6();

    double a(unsigned int) const;
    double b(unsigned int) const;
    double c(unsigned int) const;
    double x(unsigned int) const;

    double h = 0.01;
    unsigned int N = 100;

    void static Main(int argc, char **argv);
};

#endif // EXAMPLE5_H
