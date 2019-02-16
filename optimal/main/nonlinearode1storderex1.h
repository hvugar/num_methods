#ifndef NONLINEARODE1STORDEREX2_H
#define NONLINEARODE1STORDEREX2_H

#include <ode/nlode1o.h>
#include <utils/random.h>
#include <float.h>
#include <math.h>

class NonLinearODE1stOrderEx2 : public NonLinearODE1stOrder
{
public:
    static void Main(int argc, char** argv);

    NonLinearODE1stOrderEx2();
    virtual ~NonLinearODE1stOrderEx2();

protected:
    virtual unsigned int count() const;
    virtual double x(const PointNodeODE &node, unsigned int row = 1) const;
    virtual double f(const PointNodeODE &node, const DoubleVector &x, unsigned int r) const;
};

#endif // NONLINEARODE1STORDEREX2_H
