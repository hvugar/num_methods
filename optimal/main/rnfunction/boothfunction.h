#ifndef BOOTHFUNCTION_H
#define BOOTHFUNCTION_H

#include <function.h>
#include <printer.h>
#include <projection.h>

/**
 * @brief The Booth's function. Range -10.0 <= x,y <= +10.0. Optimal f(1, 3)=0;
 */
struct BoothFunction : public RnFunction, public IGradient, public Printer, public Projection
{
public:
    virtual ~BoothFunction() {}
    //RnFunction
    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g);
    //virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step);
    //Printer
    void print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const;
    //Projection
    virtual void project(DoubleVector &x, int index);

    static void main();

    double a;
    double b;

private:
    double grad_step;
};

#endif // BOOTHFUNCTION_H
