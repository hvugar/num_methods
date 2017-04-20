#ifndef ART_PROBLEM1L3_H
#define ART_PROBLEM1L3_H

#include "../iproblem1.h"
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>

class ArticleProblem1L3 : protected IProblem1, public IPrinter, public IProjection
{
public:
    static void Main(int argc, char *argv[]);

public:
    ArticleProblem1L3();
    virtual ~ArticleProblem1L3() {}

    void optimize(DoubleVector &y0) const;

    virtual void print(unsigned int i, const DoubleVector &y, const DoubleVector &g, double fx, GradientMethod::MethodResult result) const;
    virtual void project(DoubleVector &x, int index);
};

#endif // ART_PROBLEM1L3_H
