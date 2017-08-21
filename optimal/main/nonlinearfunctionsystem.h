#ifndef NONLINEARFUNCTIONSYSTEM_H
#define NONLINEARFUNCTIONSYSTEM_H

#include <function.h>
#include <printer.h>

class INonLinearFunctionSystem : public VectorRnFunction
{
protected:
    virtual double fx(const DoubleVector &x, unsigned int num) const = 0;

public:
    void calculateSimpleIdetartion(const DoubleVector &x0, DoubleVector &x, double epsilon);
    void calculateNewtonMethod(const DoubleVector &x0, DoubleVector &rx, double diffEspilon, double espilon);
};

class NonLinearFunctionSystem : public INonLinearFunctionSystem
{
public:
    static void Main(int agrc, char *argv[]);
    virtual double fx(const DoubleVector &x, unsigned int num) const;
};

#endif // NONLINEARFUNCTIONSYSTEM_H
