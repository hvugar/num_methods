#ifndef NONLINEARFUNCTIONSYSTEM_H
#define NONLINEARFUNCTIONSYSTEM_H

#include <function.h>
#include <printer.h>

class NonLinearFunctionSystem : public VectorRnFunction
{
public:
    static void Main(int agrc, char *argv[]);

    virtual double fx(const DoubleVector &x, unsigned int num) const = 0;

    void calculate(const DoubleVector &x0, DoubleVector &x, double epsilon);
};

#endif // NONLINEARFUNCTIONSYSTEM_H
