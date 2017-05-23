#ifndef EXAMPLE7_H
#define EXAMPLE7_H

#include <vector2d.h>
#include <matrix2d.h>
#include <math.h>
#include <printer.h>
#include <cmethods.h>

class Example7
{
public:
    void static Main(int argc, char *argv[]);

    Example7();

    unsigned int N;
    double h;
    unsigned int n;
    unsigned int K;
    unsigned int w;
    unsigned int p;

    void calculateK4_L_2_R();
    void calculateK4_R_2_L();

    double a(unsigned int k) const;
    double b(unsigned int k) const;
    double f(unsigned int k) const;
};

#endif // EXAMPLE7_H
