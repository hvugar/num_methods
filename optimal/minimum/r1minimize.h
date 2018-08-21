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

class MINIMUMSHARED_EXPORT R1FxMinimizer
{
public:
    struct MINIMUMSHARED_EXPORT Callback
    {
        friend class R1FxMinimizer;
        R1Function* function() const;

        virtual void straightLineSearchCallback(unsigned int iteration, double x, double a, double b, double fxa, double fxb, unsigned int fx_count) const;
        virtual void swannCallback(unsigned int iteration, double x, double a, double b, double fxa, double fxb, unsigned int fx_count) const;

        virtual void goldenSectionSearchCallback(unsigned int iteration, double x, double a, double b, double fxx, double fxa, double fxb, unsigned int fx_count) const;
        virtual void halphIntervalMethodCallback(unsigned int iteration, double x, double a, double b, double fxx, double fxa, double fxb, unsigned int fx_count) const;
        virtual void dichotomyCallback(unsigned int iteration, double x, double a, double b, double fxx, double fxa, double fxb, unsigned int fx_count) const;
        virtual void uniformLineSearchCallback(unsigned int iteration, double x, double a, double b, double fxx, double fxa, double fxb, unsigned int fx_count) const;
        virtual void fibonachiMethodCallback(unsigned int iteration, double x, double a, double b, double fxx, double fxa, double fxb, unsigned int fx_count) const;

    private:
        R1Function* mfunction;
    };

    enum UnimodalResultResult
    {
        UnimodalFunction = 0,
        NonUnimodalFunction = 1
    };

    void setFunction(R1Function *function);
    R1Function* function() const;

    /**
     * @brief Этап устонавления границ интервала
     */
    void straightLineSearch(double x, double step, double &a, double &b, double &fxa, double &fxb, bool &unimodal) const;
    void swann(double x, double step, double &a, double &b, double &fx, double &fxb, bool &unimodal) const;

    /**
     * @brief Этап уменьшения интервала
     */
    void goldenSectionSearch(double &x, double &a, double &b, double epsilon) const;
    void halphIntervalMethod(double &x, double &a, double &b, double epsilon) const;
    void dichotomyMethod(double &x, double &a, double &b, double epsilon, double step) const;
    void uniformLineSearch(double &x, double &a, double &b, unsigned int n) const;
    void fibonachiMethod(double &x, double &a, double &b, double epsilon, double step) const;

    void setCallback(R1FxMinimizer::Callback *callback);
    R1FxMinimizer::Callback* callback() const;

private:
    R1Function *mfunction = NULL;
    R1FxMinimizer::Callback *mcallback = NULL;
};

#endif // R1MINIMIZE_H
