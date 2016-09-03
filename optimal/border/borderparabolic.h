#ifndef BORDERPARABOLIC_H
#define BORDERPARABOLIC_H

#include <parabolicequation.h>

class BorderParabolic : public IParabolicEquation
{
public:
    BorderParabolic();
    virtual ~BorderParabolic();

    virtual double initial(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    static void main(int argc, char** argv);

    double t0;
    double t1;
    double x0;
    double x1;
    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
};

#endif // BORDERPARABOLIC_H
