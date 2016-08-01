#ifndef PENALTYMETHOD_H
#define PENALTYMETHOD_H

#include "global.h"
#include "function.h"
#include <vector>

class MINIMUMSHARED_EXPORT PenaltyMethod
{
public:
    PenaltyMethod();
    virtual ~PenaltyMethod();

    const std::vector<RnFunction*>& h() const;
    const std::vector<RnFunction*>& g() const;

    virtual void calculate();

    virtual double P(double x);

private:
    RnFunction* m_f;
    std::vector<RnFunction*> m_h;
    std::vector<RnFunction*> m_g;
    double m_R;
    double m_r;
};

#endif // PENALTYMETHOD_H
