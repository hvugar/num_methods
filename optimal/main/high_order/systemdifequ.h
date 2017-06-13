#ifndef EXAMPLEM_H
#define EXAMPLEM_H

#include <vector2d.h>
#include <matrix2d.h>
#include <math.h>
#include <printer.h>
#include <cmethods.h>

#define SAMPLE_1

class SystemDifEquation
{
public:
    unsigned int N;
    double h;
    unsigned int w;
    unsigned int p;

    void static Main(int argc, char *argv[]);

    SystemDifEquation();

    virtual double a(unsigned int i, unsigned int j, unsigned int k) const;
    virtual double b(unsigned int i, unsigned int k) const;
    virtual double f(unsigned int i, unsigned int k) const;

    void calculate2R2LV1(const DoubleMatrix &rx);
    void calculate4R2LV1(const DoubleMatrix &rx);
    void calculate6R2LV1(const DoubleMatrix &rx);

    void calculateRX(DoubleMatrix &rx);
    //double calculateEta();
};

#endif // EXAMPLEM_H
