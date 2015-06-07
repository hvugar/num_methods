#ifndef R1FUNCTION_H
#define R1FUNCTION_H

#include <vector>

using namespace std;

struct R1Function
{
    virtual double fx(double x) = 0;
};

struct R2Function
{
    virtual double fx(double x, double y) = 0;
};

struct RnFunction
{
    virtual double fx(const std::vector<double>& x) = 0;
};

#endif // R1FUNCTION_H
