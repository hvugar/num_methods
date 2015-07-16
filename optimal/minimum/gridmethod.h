#ifndef GRIDMETHOD_H
#define GRIDMETHOD_H

#include "function.h"
#include "doublevector.h"

class MINIMUMSHARED_EXPORT GridMethod
{
public:
    GridMethod();

    void implicitDifferenceScheme();
    void explicitDifferenceScheme();
    void tomas_algorithm(std::vector<DoubleVector> &a, const DoubleVector& b, DoubleVector& x);

    std::vector<DoubleVector> u;

    void setF(R2Function *f);

    R2Function* f;
    R1Function* m1;
    R1Function* m2;
    R1Function* fi;

    double t0;
    double t1;

    double x0;
    double x1;

    int m; // time layer count
    int n; // length part count

    double dt;
    double dx;

    double alpha;

};

#endif // GRIDMETHOD_H
