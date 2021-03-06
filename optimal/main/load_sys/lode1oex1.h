#ifndef LOADLINEARODE1STORDEREX1_H
#define LOADLINEARODE1STORDEREX1_H

#define EXAMPLE_2

#include <ode/lode1o.h>
#include <utils/random.h>

class LinearODE1stOrderEx1 : public LinearODE1stOrder
{
public:
    static void Main(int agrc, char *argv[]);

    void example1();
    void example2();
    void example3();
    void example4();
    void example5();
    void example6();
    void example7();
    void example8();

    double X(double t, int i=0) const;
public:
    virtual double A(const PointNodeODE &node, unsigned int row = 0, unsigned int col = 0) const;
    virtual double B(const PointNodeODE &node, unsigned int row = 0) const;
    virtual unsigned int count() const;
};

#endif // LOADLINEARODE1STORDEREX1_H
