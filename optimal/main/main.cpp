#include <stdio.h>
#include <function.h>
#include <r1minimize.h>
#include <gradient.h>
#include <methods.h>

class R1Function1 : public R1Function
{
public:
    virtual double fx(double x) { return x*x; }
};

class RnFunction1 : public RnFunction
{
public:
    virtual double fx(double *x, int n) { return x[0]*x[0] + 2*x[1]*x[1]; }
    virtual double fx(std::vector<double> x) { return 0.0; }
};

class R1Minimize1 : public R1Minimize
{
public:
    virtual double fx(double x) { return (x+1.0)*(x+1.0) + 1.0; }
};

class SampleGradient : public Gradient1
{
protected:
    virtual double fx(std::vector<double> x)
    {
        double x1 = x[0];
        double x2 = x[1];
        return ((1-x1)*(1-x1)) + 100*(x2-x1*x1)*(x2-x1*x1);
    }
};

int main()
{
    std::vector<double> x;
    x.push_back(-1.0);
    x.push_back(+1.2);

    SampleGradient gm;
    gm.setEpsilon(0.00001);
    gm.setPoint(x);
    gm.fastProximalGradientMethod();
}

