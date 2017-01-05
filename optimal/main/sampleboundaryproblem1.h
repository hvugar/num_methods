#ifndef SAMPLEBOUNDARYPROBLEM1_H
#define SAMPLEBOUNDARYPROBLEM1_H

#include <boundaryproblem.h>
#include <printer.h>
#include <math.h>

#define SAMPLE_1

class SampleBoundaryProblem1 : public BoundaryProblem
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

#endif // SAMPLEBOUNDARYPROBLEM1_H
