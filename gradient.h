#ifndef GRADIENT_H
#define GRADIENT_H

#include "r1minimize.h"
#include "doublevector.h"

class RnFunction;
class IGradient;
class IPrinter;
class IProjection;

class GradientMethod
{
public:
    GradientMethod();
    virtual ~GradientMethod();

    virtual void calculate(DoubleVector &x) = 0;

    virtual RnFunction* function() const;
    virtual void setFunction(RnFunction *function);

    virtual IGradient* gradient() const;
    virtual void setGradient(IGradient *gradient);

    double epsilon1() const;
    void setEpsilon1(double epsilon);

    double epsilon2() const;
    void setEpsilon2(double epsilon);

    double epsilon3() const;
    void setEpsilon3(double epsilon);

    void setR1MinimizeEpsilon(double step, double epsilon);

    int count() const;

    void setPrinter(IPrinter *printer);
    void setProjection(IProjection *projection);
    void setNormalize(bool normalize);
    void showEndMessage(bool showEndMessage);

protected:
    virtual double minimize(const DoubleVector &x, const DoubleVector &g) = 0;

    RnFunction *m_fn;
    IGradient *m_gr;
    double m_epsilon1;
    double m_epsilon2;
    double m_epsilon3;
    double min_epsilon;
    double min_step;
    int iterationCount;
    bool m_normalize;
    IPrinter* m_printer;
    IProjection *m_projection;
    bool mshowEndMessage;
};

#endif // GRADIENT_H
