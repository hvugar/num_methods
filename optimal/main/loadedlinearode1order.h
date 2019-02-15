#ifndef LOADEDLINEARODE1ORDER_H
#define LOADEDLINEARODE1ORDER_H

#include <ode/lode1o.h>

class LoadedLinearODE1Order : public LinearODE1stOrder
{
public:
    static void Main(int argc, char *argv[]);

    void solveProblem(std::vector<Condition> &a, std::vector<Condition> &b);
protected:
    virtual double A(const PointNodeODE &node, unsigned int row = 0, unsigned int col = 0) const;
    virtual double B(const PointNodeODE &node, unsigned int row = 0) const;
    virtual double C(const PointNodeODE &node, unsigned int s, unsigned int row = 0, unsigned int col = 0) const;

public:
    double x(const PointNodeODE &node, unsigned int num) const;
};

#endif // LOADEDLINEARODE1ORDER_H
