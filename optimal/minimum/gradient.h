#ifndef GRADIENT_H
#define GRADIENT_H

#include "global.h"

#include <vector>
#include <math.h>

#include "function.h"
#include "r1minimize.h"
#include "doublevector.h"
#include "printer.h"
#include "projection.h"
#include "exceptions.h"

/**
 * @brief The Abstract Gradient Method class
 */
class MINIMUMSHARED_EXPORT GradientMethod
{
public:
    GradientMethod();
    virtual ~GradientMethod();

    virtual void calculate(DoubleVector& x) = 0;

    virtual RnFunction* function() const;
    virtual void setFunction(RnFunction* function);

    virtual IGradient* gradient() const;
    virtual void setGradient(IGradient *gradient);

    /**
     * @brief Epsilon for gradient norm
     * @return
     */
    double epsilon1() const;
    void setEpsilon1(double epsilon);

    /**
     * @brief Epsilon for points distance
     * @return
     */
    double epsilon2() const;
    void setEpsilon2(double epsilon);

    /**
     * @brief Epsilon for function result
     * @return
     */
    double epsilon3() const;
    void setEpsilon3(double epsilon);

    void setR1MinimizeEpsilon(double step, double epsilon1);
    int count() const;

    void setPrinter(Printer* printer);
    void setProjection(Projection* projection);
    void setNormalize(bool normalize);
    void showEndMessage(bool showEndMessage);

protected:
    virtual double minimize(const DoubleVector &x, const DoubleVector &g) = 0;

    RnFunction *m_fn;
    IGradient *m_gr;
    double m_epsilon1;
    double m_epsilon2;
    double m_epsilon3;
    //double grad_step;
    double min_epsilon;
    double min_step;
    int iterationCount;
    bool m_normalize;
    Printer* m_printer;
    Projection *m_projection;

    bool mshowEndMessage;
};

#endif // GRADIENT_H
