#ifndef NON_LINEAR_EQUATION_EX1_H
#define NON_LINEAR_EQUATION_EX1_H

#include <nonlinearequation.h>

class NonLinearEquationEx1 : public NonLinearEquation
{
public:
    static void Main(int agrc, char *argv[]);

    virtual double fx(const DoubleVector &x, unsigned int num = 0) const;
};

#endif // NON_LINEAR_EQUATION_EX1_H
