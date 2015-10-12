#ifndef HEATCONTROL2DELTA_H
#define HEATCONTROL2DELTA_H

#include <function.h>

class HeatControl2Delta : public RnFunction
{
public:
    HeatControl2Delta();
    ~HeatControl2Delta();

    double fx(const DoubleVector& x);
    void gradient(double step, const DoubleVector& x, DoubleVector& g);

private:

};

#endif // HEATCONTROL2DELTA_H
