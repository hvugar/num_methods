#ifndef R1MINIMIZE_H
#define R1MINIMIZE_H

#include "global.h"

class MINIMUMSHARED_EXPORT R1Minimize
{
public:
    R1Minimize();
    virtual ~R1Minimize();

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

protected:
    virtual double fx(double x) = 0;

private:
    double m_x0;
    double m_step;
    double m_epsilon;
    double m_a;
    double m_b;
};

#endif // R1MINIMIZE_H
