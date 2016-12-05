#ifndef EXAMPLE4_H
#define EXAMPLE4_H

#include <vector2d.h>
#include <matrix2d.h>
#include <matrix3d.h>
#include <printer.h>
#include <math.h>
#include <vector>

class Example4
{
public:
    Example4();

    double h;
    unsigned int N;
    unsigned int K;
    unsigned int n;

    void static Main(int argc, char *argv[]);

    double X1(unsigned int) const;
    double X2(unsigned int) const;
    double X3(unsigned int) const;

    void calculate1();
    void calculate2();
    void calculate3();
    void calculate3M();

    double a(unsigned int i, unsigned int j, unsigned int k) const;
    double b(unsigned int i, unsigned int k) const;
};

#endif // EXAMPLE4_H
