#ifndef PROBLEM1L3_H
#define PROBLEM1L3_H

#include <function.h>
#include <gradient.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>
#include <cmethods.h>

#include "iproblem1.h"

class Problem1L3 : public IProblem1, public IProjection
{
public:
    Problem1L3();
    virtual void project(DoubleVector &x, unsigned int index);
    static void Main(int argc, char* argv[]);
};

#endif // PROBLEM1L3_H
