#ifndef SAMPLEBORDERHYPERBOLIC_H
#define SAMPLEBORDERHYPERBOLIC_H

#include <hyperbolicequation.h>

class SampleBorderHyperBolic : public IHyperbolicEquation
{
public:
    SampleBorderHyperBolic();
    virtual ~SampleBorderHyperBolic() {}

    virtual double f(unsigned int i, unsigned int j) const;
    virtual double initial1(unsigned int i) const;
    virtual double initial2(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;

    static void main(int argc, char **argv);

private:
    double x0;
    double x1;
    double t0;
    double t1;
    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
    double a;
};

#endif // SAMPLEBORDERHYPERBOLIC_H
