#ifndef GRADIENT_H
#define GRADIENT_H

#include "global.h"
#include <vector>

using namespace std;

class MINIMUMSHARED_EXPORT Gradient
{
public:
    Gradient();

    void setPoint(const std::vector<double>& x);
    const std::vector<double>& x() const;
    void setEpsilon(double);
    double epsilon() const;

    virtual void gradient();
    virtual double minimize();
    void fastProximalGradientMethod();

protected:
    virtual double fx(std::vector<double> x) = 0;
    virtual void iterationInfo();

private:
    std::vector<double> mx;
    std::vector<double> mgrads;
    double eps;
};

#endif // GRADIENT_H
