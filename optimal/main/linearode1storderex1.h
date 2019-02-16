#ifndef LINEARODE1STORDEREX1_H
#define LINEARODE1STORDEREX1_H

#include <ode/lode1o.h>
#include <utils/random.h>
#include <float.h>
#include <math.h>

class LinearODE1stOrderEx1 : public LinearODE1stOrder
{
public:
    LinearODE1stOrderEx1();
    virtual ~LinearODE1stOrderEx1();

    static void Main(int argc, char** argv);

protected:
    virtual double A(const PointNodeODE &node, unsigned int row = 0, unsigned int col = 0) const;
    virtual double B(const PointNodeODE &node, unsigned int row = 0) const;
    virtual double x(const PointNodeODE &node, unsigned int row = 0) const;
};

#endif // LINEARODE1STORDEREX1_H
