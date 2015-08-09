#ifndef R1MINIMIZE_H
#define R1MINIMIZE_H

#include "global.h"
#include "function.h"

class MINIMUMSHARED_EXPORT R1Minimize
{
public:
    R1Minimize();
    virtual ~R1Minimize();

    void setFunction(R1Function *f);
    R1Function* function() const;

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

    static void StranghLineSearch(double x, double step, double &a, double &b, R1Function *f);
    static void Swann(double x, double step, double &a, double &b, R1Function *f);

    static double GoldenSectionSearch(double a, double b, double &x, R1Function *f, double epsilon);
    static double HalphIntervalMethod(double a, double b, double &x, R1Function *f, double epsilon);

private:
    R1Function *m_f;
    double m_x0;
    double m_step;
    double m_epsilon;
    double m_a;
    double m_b;
};

#endif // R1MINIMIZE_H
