#ifndef EXAMPLE6_H
#define EXAMPLE6_H

#include <vector2d.h>
#include <matrix2d.h>
#include <matrix3d.h>
#include <printer.h>

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class Example6
{
public:
    static void Main(int argc, char* arg[]);

public:
    Example6();

    double a(unsigned int i, unsigned int j, unsigned int k) const;
    double b(unsigned int i, unsigned int k) const;
    double fx(unsigned int i, unsigned int k) const;

private:
    double h = 0.01;
    unsigned int N = 100;
    unsigned int K = 4;
    unsigned int F = 10;
};

#endif // EXAMPLE6_H
