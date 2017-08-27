#ifndef SINGLE_DIF_EQUATION_H
#define SINGLE_DIF_EQUATION_H

#include <vector2d.h>
#include <matrix2d.h>
#include <math.h>
#include <printer.h>
#include <cmethods.h>
#include <grid/ibvp.h>

#define SAMPLE_1

class SingleDifEquation
{
public:
    void static Main(int argc, char *argv[]);

    SingleDifEquation();

    virtual double f(unsigned int k) const;
    virtual double a(unsigned int k) const;
    virtual double b(unsigned int k) const;

    unsigned int N;
    double h;
    unsigned int w;
    unsigned int p;

    void calculate2R2LV1(const DoubleVector &rx);
    void calculate2R2LV11(const DoubleVector &rx);
    void calculate4R2LV1(const DoubleVector &rx);
    void calculate4R2LV11(const DoubleVector &rx);
    void calculate6R2LV1(const DoubleVector &rx);
    void calculate6R2LV11(const DoubleVector &rx);

    void calculateRX(DoubleVector &rx);
    double calculateEta();
};

#endif // SINGLE_DIF_EQUATION_H
