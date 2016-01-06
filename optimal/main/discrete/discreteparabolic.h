#ifndef DISCRETEPARABOLIC_H
#define DISCRETEPARABOLIC_H

#include <function.h>
#include <doublevector.h>

class DiscreteParabolic : public RnFunction, public R1Function
{
public:
    DiscreteParabolic();
    ~DiscreteParabolic();

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step=0.000001);
    virtual double fx(double x);

    double fi(unsigned int i);
    double m1(unsigned int j);
    double m2(unsigned int j);
    double f(unsigned int i, unsigned int j);

    void calculateU(const DoubleVector& f, DoubleVector& u);

    static void main();

private:
    DoubleVector *pf;

    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
    double x0;
    double x1;
    double t0;
    double t1;
    double a;
    double lamda;
};

#endif // DISCRETEPARABOLIC_H
