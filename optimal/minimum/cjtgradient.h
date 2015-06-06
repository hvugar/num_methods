#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "gradient.h"

class MINIMUMSHARED_EXPORT ConjugateGradient : public Gradient
{
public:
    ConjugateGradient();
    virtual ~ConjugateGradient();

    virtual void calculate();
    virtual double minimize();
    void print();

private:
    std::vector<double> s;
    struct ArgMin : public R1Function
    {
        ArgMin(std::vector<double> &x, std::vector<double> &g, Gradient *gradient);
        std::vector<double> &x;
        std::vector<double> &g;
        Gradient *gradient;
        virtual double fx(double alpha);
    };
};

#endif // CONJUGATEGRADIENT_H
