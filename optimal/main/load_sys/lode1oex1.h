#ifndef LOADLINEARODE1STORDEREX1_H
#define LOADLINEARODE1STORDEREX1_H

#define EXAMPLE_4

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

    double X(double t, int i=0) const;
public:
    virtual double A(const GridNodeODE &node, unsigned int row = 0, unsigned int col = 0) const;
    virtual double B(const GridNodeODE &node, unsigned int row = 0) const;
    virtual unsigned int equationsNumber() const;
};

#endif // LOADLINEARODE1STORDEREX1_H
