#ifndef PENALTYMETHOD_H
#define PENALTYMETHOD_H

#include "global.h"
#include "function.h"

class MINIMUMSHARED_EXPORT PenaltyMethod
{
public:
    PenaltyMethod();
    virtual ~PenaltyMethod();

    const PRnFunctionList& h() const;
    const PRnFunctionList& g() const;

    virtual void calculate();

private:
    RnFunction* m_f;
    PRnFunctionList m_h;
    PRnFunctionList m_g;
    double m_R;
    double m_r;
};

#endif // PENALTYMETHOD_H
