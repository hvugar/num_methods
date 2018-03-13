#ifndef JFUNCTIONAL_H
#define JFUNCTIONAL_H

#include "ifunctional.h"

class JFunctional : public IFunctional
{
public:
    JFunctional();

    virtual double fx(const DoubleVector &prms) const;
    virtual void gradient(const DoubleVector &prms, DoubleVector &g) const;

    void setInitTemperatures(const DoubleVector &fis, const DoubleVector &p_fis);
    void setEnvrTemperatures(const DoubleVector &thetas, const DoubleVector &p_thetas);

    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, GradientMethod::MethodResult result) const;

private:
    DoubleVector fis;
    DoubleVector p_fis;
    DoubleVector thetas;
    DoubleVector p_thetas;
};

#endif // JFUNCTIONAL_H
