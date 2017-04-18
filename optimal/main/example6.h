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

#define SAMPLE_2

class Example6
{
public:
    static void Main(int argc, char* arg[]);

public:
    Example6();

    double a(unsigned int i, unsigned int j, unsigned int k) const;
    double b(unsigned int i, unsigned int k) const;
    double fx(unsigned int i, unsigned int k) const;

    void method1K2();

private:
    double h = 0.001;
    unsigned int N = 1000;
    unsigned int K = 4;
    unsigned int F = 10;
    unsigned int n = 3;
};

#endif // EXAMPLE6_H
