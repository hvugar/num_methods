#ifndef BORDERTEST_H
#define BORDERTEST_H

#include <parabolicequation.h>

class BorderTest : public IParabolicEquation
{
public:
    BorderTest();
    virtual ~BorderTest() {}

    virtual double initial(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    double hx = 0.001;
    double ht = 0.01;
    unsigned int N = 1000;
    unsigned int M = 100;
    double a = 1.0;

    double U(unsigned int i, unsigned int j) const;

    virtual void calculateUError2(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;
    virtual void calculateU3Forward(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;
    virtual void calculateU3Backward(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a=1.0) const;

    static void Main(int argc, char* argv[]);
};

#endif // BORDERTEST_H
