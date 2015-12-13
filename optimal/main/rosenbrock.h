#ifndef ROSENBROCK_H
#define ROSENBROCK_H

#include <function.h>
#include <printer.h>

struct Rosenbrock : public RnFunction, public Printer
{
public:
    Rosenbrock() : count(0) {}
    virtual ~Rosenbrock() {}

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step);

    void print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const;

    static void main();
private:
    double grad_step;
    unsigned int count;
};

#endif // ROSENBROCK_H
