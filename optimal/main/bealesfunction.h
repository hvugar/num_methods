#ifndef BEALESFUNCTION_H
#define BEALESFUNCTION_H

#include <function.h>
#include <printer.h>

/**
 * @brief The Beales Function. Range -4.5 <= x,y <= 4.5
 */
struct BealesFunction : public RnFunction
{
public:
    virtual double fx(const DoubleVector& x);
    virtual void gradient(double gradient_step, const DoubleVector& x, DoubleVector& g);

    static void main();

private:
    double grad_step;
};

struct BealesPrinter : public Printer
{
    void print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const;
};

#endif // BEALESFUNCTION_H
