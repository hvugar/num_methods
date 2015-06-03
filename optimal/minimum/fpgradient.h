#ifndef FPGRADIENT_H
#define FPGRADIENT_H

#include "gradient.h"

class MINIMUMSHARED_EXPORT FastProximalGradient : public Gradient
{
public:
    FastProximalGradient();
    virtual ~FastProximalGradient();

    virtual void setF(RnFunction* f);
    virtual RnFunction *f() const;

    virtual void gradient();
    virtual double argmin(double alpha);
    virtual double minimize();
    void calculate();

private:
    RnFunction *mf;
    std::vector<double> mgrads;
};

#endif // FPGRADIENT_H
