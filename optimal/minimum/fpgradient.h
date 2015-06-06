#ifndef FPGRADIENT_H
#define FPGRADIENT_H

#include "gradient.h"

class MINIMUMSHARED_EXPORT FastProximalGradient : public Gradient
{
public:
    FastProximalGradient();
    virtual ~FastProximalGradient();

    virtual double minimize();
    virtual void calculate();
    void print();

private:
    struct ArgMin : public R1Function
    {
        ArgMin(std::vector<double> &x, std::vector<double> &g, Gradient *gradient);
        std::vector<double> &x;
        std::vector<double> &g;
        Gradient *gradient;
        virtual double fx(double alpha);
    };
};

#endif // FPGRADIENT_H
