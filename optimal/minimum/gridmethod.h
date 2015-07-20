#ifndef GRIDMETHOD_H
#define GRIDMETHOD_H

#include "function.h"
#include "doublevector.h"
#include <stdio.h>

class MINIMUMSHARED_EXPORT GridMethod
{
public:
    GridMethod();

    void implicitDifferenceScheme();
    void explicitDifferenceScheme();
    void tomas_algorithm(std::vector<DoubleVector> &a, const DoubleVector& b, DoubleVector& x);

    void setLengthInterval(double x0, double x1);
    void setTimeInterval(double t0, double t1);
    void setLengthTimeStep(double dx, double dt);
    void setLengthTimeStepCount(unsigned int n, unsigned int m);

    std::vector<DoubleVector> u;

    void setF(R2Function* f);
    void setM1(R1Function* m1);
    void setM2(R1Function* m2);
    void setFi(R1Function* fi);

private:
    double t0;
    double t1;

    double x0;
    double x1;

    unsigned int m; // time layer count
    unsigned int n; // length part count

    double dt;
    double dx;

    double alpha;

    R2Function* f;
    R1Function* m1;
    R1Function* m2;
    R1Function* fi;
};

#endif // GRIDMETHOD_H
