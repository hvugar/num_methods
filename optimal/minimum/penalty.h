#ifndef PENALTYMETHOD_H
#define PENALTYMETHOD_H

#include "global.h"
#include "function.h"

class MINIMUMSHARED_EXPORT PenaltyMethod
{
public:
    PenaltyMethod();
    virtual ~PenaltyMethod();

    std::vector<RnFunction*>& h() const;
    std::vector<RnFunction*>& g() const;

    virtual void calculate();

private:
    RnFunction* m_f;
    std::vector<RnFunction*> m_h;
    std::vector<RnFunction*> m_g;
    double m_R;
    double m_r;
};

#endif // PENALTYMETHOD_H
