#ifndef BOOTHFUNCTION_H
#define BOOTHFUNCTION_H

#include <function.h>
#include <printer.h>

/**
 * @brief The Booth's function. Range -10.0 <= x,y <= +10.0. Optimal f(1, 3)=0;
 */
struct BoothFunction : public RnFunction
{
public:
    virtual double fx(const DoubleVector& x);
    virtual void gradient(double gradient_step, const DoubleVector& x, DoubleVector& g);

    static void main();

private:
    double grad_step;
};

struct BoothPrinter : public Printer
{
    void print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const;
};

#endif // BOOTHFUNCTION_H
