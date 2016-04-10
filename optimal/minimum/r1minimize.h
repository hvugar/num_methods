#ifndef R1MINIMIZE_H
#define R1MINIMIZE_H

#include "global.h"
#include "function.h"

#ifdef __cplusplus
extern "C" {
#endif
void MINIMUMSHARED_EXPORT stranghLineSearch(double x, double step, double &a, double &b, R1Function *fn);
double MINIMUMSHARED_EXPORT goldenSectionSearch(double &a, double &b, double &x, R1Function *f, double epsilon);
double MINIMUMSHARED_EXPORT goldenSectionSearch1(double &a, double &b, double &x, R1Function *f, double epsilon);
#ifdef __cplusplus
}
#endif

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

    static double GoldenSectionSearch(double &a, double &b, double &x, R1Function *f, double epsilon);
    static void FibonachiMethod(double &a, double &b, double &c, double step, double epsilon, R1Function *f);
    static double HalphIntervalMethod(double &a, double &b, double &x, R1Function *f, double epsilon);
    static void UniformLineSearch(double &a, double &b, double &c, unsigned int n, R1Function *f);
    static void BruteForceLineSearch(double &a, double &b, double &c, unsigned int n, R1Function *f);
    static void DichotomyMethod(double &a, double &b, double &c, R1Function *f, double step, double epsilon);

private:
    R1Function *m_f;
    double m_x0;
    double m_step;
    double m_epsilon;
    double m_a;
    double m_b;
};

#endif // R1MINIMIZE_H
