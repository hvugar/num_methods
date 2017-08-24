#ifndef LOADLINEARODE1STORDEREX1_H
#define LOADLINEARODE1STORDEREX1_H

#include <ode/lode1o.h>

class LoadLinearODE1stOrderEx1 : public LinearODE1stOrder
{
public:
    static void Main(int agrc, char *argv[]);

    void initialize(std::vector<Condition> &nscs, DoubleVector &betta);
    double X(double t, int i) const;

public:
    virtual double A(double t UNUSED_PARAM, unsigned int k, unsigned int row = 0, unsigned int col = 0) const;
    virtual double B(double t UNUSED_PARAM, unsigned int k, unsigned int row = 0) const;
};

#endif // LOADLINEARODE1STORDEREX1_H
