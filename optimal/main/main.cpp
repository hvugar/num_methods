#include <stdio.h>
#include <function.h>
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

int main()
{
    R1Function1 r1;
    double a, b;
    a=b=0.0;
    straight_line_search_metod(&r1, 10.0, 0.1, a, b);
    double min = golden_section_search_min(&r1, a, b, 0.00001);
    printf("%f %f %f\n", a, b, min);

    RnFunction1 r2;
    double x[] = { 2.0, 2.0 };
    double g[] = { 0.0, 0.0 };
    gradient(&r2, x, 2, 0.00001, g);
    printf("%f %f\n", g[0], g[1]);

    return 0;
}

