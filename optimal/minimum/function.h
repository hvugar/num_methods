#ifndef R1FUNCTION_H
#define R1FUNCTION_H

#include <vector>

using namespace std;

class R1Function
{
public:
    virtual double fx(double x) = 0;
};

class R2Function
{
public:
    virtual double fx(double x, double y) = 0;
};

class RnFunction
{
public:
    virtual double fx(double *x, int n) = 0;
    virtual double fx(std::vector<double> x) = 0;
};

#endif // R1FUNCTION_H
