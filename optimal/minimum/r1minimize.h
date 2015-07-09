#ifndef R1MINIMIZE_H
#define R1MINIMIZE_H

#include "function.h"
#include "r1minimize.h"
#include "global.h"

class MINIMUMSHARED_EXPORT R1Minimize
{
public:
    R1Minimize();
    virtual ~R1Minimize();

    void setF(R1Function *f);
    R1Function* f() const;

    double x0() const;
    void setX0(double);

    double step() const;
    void setStep(double);

    double epsilon() const;
    void setEpsilon(double);

    double a() const;
    double b() const;

    double straightLineSearch();

    double goldenSectionSearch();
    double halphIntervalMethod();

private:
    R1Function *m_f;
    double m_x0;
    double m_step;
    double m_epsilon;
    double m_a;
    double m_b;
};

#endif // R1MINIMIZE_H
