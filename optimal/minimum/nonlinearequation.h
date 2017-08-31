#ifndef NONLINEAREQUATION_H
#define NONLINEAREQUATION_H

#include <function.h>
#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>

class MINIMUMSHARED_EXPORT NonLinearEquation
{
protected:
    virtual double fx(const DoubleVector &x, unsigned int num = 0) const = 0;

public:
    void calculateSimpleIdetartion(const DoubleVector &x0, DoubleVector &x, double epsilon);
    void calculateNewtonMethod(const DoubleVector &x0, DoubleVector &rx, double diffEspilon, double espilon);
    void calculateNewtonMethodMod(const DoubleVector &x0, DoubleVector &rx, double diffEspilon, double espilon);
    void calculateNewtonMethodMod2(const DoubleVector &x0, DoubleVector &rx, double diffEspilon, double espilon);

private:
    double minimize(const DoubleMatrix &W, const DoubleMatrix &WI, const DoubleVector& xk, unsigned int n);
};

#endif // NONLINEAREQUATION_H
