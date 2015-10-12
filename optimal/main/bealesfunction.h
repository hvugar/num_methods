#ifndef BEALESFUNCTION_H
#define BEALESFUNCTION_H

#include <function.h>
#include <printer.h>

/**
 * @brief The Beales Function. Range -4.5 <= x,y <= 4.5. Optimal f(3, 0.5)=0;
 */
struct BealesFunction : public RnFunction
{
public:
    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step);

    static void main();

private:
    double grad_step;
};

struct BealesPrinter : public Printer
{
    void print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const;
};

#endif // BEALESFUNCTION_H
