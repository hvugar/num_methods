#ifndef HYPERBOLICCONTROL1DT_H
#define HYPERBOLICCONTROL1DT_H

#include <function.h>

class HyperbolicControl1DT : public R1Function
{
public:
    HyperbolicControl1DT();
    virtual ~HyperbolicControl1DT() {}

    virtual double fx(double x);

    void static main();

protected:
};

#endif // HYPERBOLICCONTROL1DT_H
