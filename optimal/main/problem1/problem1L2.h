#ifndef PROBLEM1L2_H
#define PROBLEM1L2_H

#include <function.h>
#include <gradient.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <cmethods.h>
#include <gradient/igradient.h>

#include "problem1newton.h"
#include "iproblem1.h"

class Problem1L2 : public IProblem1, public IPrinter, public IProjection
{
public:
    static void Main(int argc, char* argv[]);

    Problem1L2();
    virtual ~Problem1L2() {}

    void initialize();
    void startOptimize();
    void optimize(DoubleVector &x0) const;
    virtual void print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double fx, GradientMethod::MethodResult result) const;
    virtual void project(DoubleVector &x, int index);
};

#endif // PROBLEM1L2_H
