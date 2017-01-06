#ifndef BOUNDARYVALUEPROBLEM1_H
#define BOUNDARYVALUEPROBLEM1_H

#include <bvp.h>
#include <printer.h>
#include <math.h>

#define SAMPLE_1

class BoundaryValueProblem1 : public BoundaryValueProblem
{
public:
    virtual double r(unsigned int i) const;
    virtual double p(unsigned int i) const;
    virtual double q(unsigned int i) const;
    virtual double f(unsigned int i) const;
    virtual double boundary(Boundary bound) const;

    double fx(unsigned int i) const;
    double h;
    unsigned int N;

    void static Main(int argc, char** argv);
};

#endif // BOUNDARYVALUEPROBLEM1_H
