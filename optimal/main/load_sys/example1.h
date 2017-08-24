#ifndef EXAMPLE1_H
#define EXAMPLE1_H

#include <ode/lode1o.h>

class Example1 : public LinearODE1stOrder
{
public:
    void static Main(int argc, char *argv[]);

protected:
    virtual double A(double x, unsigned int i, unsigned int row = 0, unsigned int col = 0) const;
    virtual double B(double x, unsigned int i, unsigned int row = 0) const;
    virtual unsigned int equationsNumber() const;
};

#endif // EXAMPLE1_H
