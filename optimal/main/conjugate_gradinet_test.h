#ifndef CONJUGATE_GRADINET_TEST_H
#define CONJUGATE_GRADINET_TEST_H

#include <gradient_cjt.h>
#include <function.h>
#include <printer.h>
#include <projection.h>

class ConjugateGradinetTest : public RnFunction, public IGradient, public IPrinter
{
public:
    static void Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM);

    virtual double fx(const DoubleVector &x) const;
    virtual void gradient(const DoubleVector &x, DoubleVector &g) const;
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const;
};

#endif // CONJUGATE_GRADINET_TEST_H
