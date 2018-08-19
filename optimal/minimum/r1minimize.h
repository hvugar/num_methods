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
double MINIMUMSHARED_EXPORT goldenSectionSearch2(double &a, double &b, double &x, R1Function *f, double epsilon);
#ifdef __cplusplus
}
#endif

class MINIMUMSHARED_EXPORT R1Minimize1
{
public:
    R1Minimize1();
    virtual ~R1Minimize1();

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
    static void UniformLineSearch(double &a, double &b, double &c, unsigned int n, R1Function *f);

    static double GoldenSectionSearch(double &a, double &b, double &x, R1Function *f, double epsilon);
    static void FibonachiMethod(double &a, double &b, double &c, double step, double epsilon, R1Function *f);
    static double HalphIntervalMethod(double &a, double &b, double &x, R1Function *f, double epsilon);
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


class MINIMUMSHARED_EXPORT R1FxMinimizer
{
public:
    struct Callback
    {
        friend class R1FxMinimizer;

        virtual void straightLineSearchCallback(unsigned int iteration, double x, double a, double b, double fxa, double fxb, unsigned int fx_count) const;
        virtual void goldenSectionSearchCallback(unsigned int iteration, double x, double a, double b, double fxx, double fxa, double fxb, unsigned int fx_count) const;

        R1Function* function() const;

    private:
        R1Function* mfunction;
    };

    void setFunction(R1Function *function);
    R1Function* function() const;

    /**
     * @brief Этап устонавления границ интервала
     */
    void straightLineSearch(double x, double step, double &a, double &b, double &fxa, double &fxb) const;
    void swann(double x, double step, double &a, double &b, double &fx, double &fxb) const;

    /**
     * @brief Этап уменьшения интервала
     */
    void goldenSectionSearch(double &x, double &a, double &b, double epsilon) const;
    void halphIntervalMethod(double &x, double &a, double &b, double epsilon) const;
    void uniformLineSearch(double &x, double &a, double &b, unsigned int n) const;

    void setCallback(R1FxMinimizer::Callback *callback);
    R1FxMinimizer::Callback* callback() const;

private:
    R1Function *mfunction = NULL;
    R1FxMinimizer::Callback *mcallback = NULL;
};

#endif // R1MINIMIZE_H
