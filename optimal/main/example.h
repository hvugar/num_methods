#ifndef EXAMPLE_H
#define EXAMPLE_H

#include <vector2d.h>
#include <matrix2d.h>
#include <math.h>
#include <printer.h>
#include <cmethods.h>

class Example
{
public:
    void static Main(int argc, char *argv[]);

    Example();

    virtual double f(unsigned int k) const;
    virtual double a(unsigned int k) const;
    virtual double b(unsigned int k) const;

    unsigned int N;
    double h;
    unsigned int w;
    unsigned int p;

    void calculate2R2LV1();
    void calculate2R2LV2();

    void calculate4R2LV1();
    void calculate4R2LV2();
    void calculate6R2LV2();
};

#endif // EXAMPLE_H
