#ifndef BOUNDARYVALUEPROBLEM1_H
#define BOUNDARYVALUEPROBLEM1_H

#include <grid/lbvpode.h>

#define SAMPLE_6

class BoundaryValueProblem1 : public LinearBoundaryValueProblemODE
{
public:
    virtual double r(unsigned int i) const;
    virtual double p(unsigned int i) const;
    virtual double q(unsigned int i) const;
    virtual double f(unsigned int i) const;
    virtual double boundary(BoundaryType bound) const;

    double fx(unsigned int i) const;
    double h;
    unsigned int N;

    void static Main(int argc, char** argv);
};

#endif // BOUNDARYVALUEPROBLEM1_H
