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

//    static double _u(double x, double y, double t);
//    static double _f(double x, double y, double t);
//    static double _fi(double x, double y);
//    static double _m1(double y, double t);
//    static double _m2(double y, double t);
//    static double _m3(double x, double t);
//    static double _m4(double x, double t);
    static void VariableDirectionsMethod(R2Function *fi, R2Function *m1, R2Function *m2, R2Function *m3, R2Function *m4, R3Function *f);
    static void TomasAlgorithm(const DoubleVector& a, const DoubleVector& b, const DoubleVector& c, const DoubleVector& d, DoubleVector& x);
    static void printResult(unsigned int k, unsigned int N, unsigned int M);

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
