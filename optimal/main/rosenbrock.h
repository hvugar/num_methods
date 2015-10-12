#ifndef ROSENBROCK_H
#define ROSENBROCK_H

#include <function.h>
#include <printer.h>

struct Rosenbrock : public RnFunction
{
public:
    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step);

    static void main();

private:
    double grad_step;
};

struct RosenbrockPrinter : public Printer
{
    void print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const;
};

#endif // ROSENBROCK_H
