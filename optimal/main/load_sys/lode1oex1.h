#ifndef LOADLINEARODE1STORDEREX1_H
#define LOADLINEARODE1STORDEREX1_H

#define EXAMPLE_3

#include <ode/lode1o.h>

class LinearODE1stOrderEx1 : public LinearODE1stOrder
{
public:
    static void Main(int agrc, char *argv[]);

    void example1();
    void example2();
    void example3();

    double X(double t, int i) const;
public:
    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row = 0, unsigned int col = 0) const;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row = 0) const;
    virtual unsigned int equationsNumber() const;
};

#endif // LOADLINEARODE1STORDEREX1_H
