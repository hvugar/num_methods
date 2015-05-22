#include <stdio.h>
#include <function.h>
#include <r1minimize.h>
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

int main()
{
    R1Minimize1 r1;
    r1.setX0(0.0);
    r1.setStep(0.02);
    r1.setEpsilon(0.0);

    double c = 0.0;
    r1.straightLineSearch();
    printf("[%.6f, %.6f]\n", r1.a(), r1.b());
    r1.setEpsilon(0.00001);
    c = r1.goldenSectionSearch();
    printf("%.6f [%.6f, %.6f]\n", c, r1.a(), r1.b());
    return 0;
}

