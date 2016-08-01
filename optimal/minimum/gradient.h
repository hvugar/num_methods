#ifndef GRADIENT_H
#define GRADIENT_H

#include "global.h"

#include "r1minimize.h"
#include "exceptions.h"

class RnFunction;
class IGradient;
class IPrinter;
class IProjection;

struct GradientIterationInfo
{
    unsigned int number;
    DoubleVector x;
    DoubleVector g;
    DoubleVector s;
    double alpha;
    double fxResult;
    RnFunction *fx;
    double gradientNorm;
    double distance;
    double difference;
};

/**
 * @brief The Abstract Gradient Method class
 */
class MINIMUMSHARED_EXPORT GradientMethod
{
public:
    GradientMethod();
    virtual ~GradientMethod();

    virtual void calculate(DoubleVector &x) = 0;

    virtual RnFunction* function() const;
    virtual void setFunction(RnFunction *function);

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

    /**
     * @brief setR1MinimizeEpsilon
     * @param step
     * @param epsilon1
     */
    void setR1MinimizeEpsilon(double step, double epsilon);

    /**
     * @brief count
     * @return
     */
    int count() const;

    /**
     * @brief setPrinter
     * @param printer
     */
    void setPrinter(IPrinter *printer);

    /**
     * @brief setProjection
     * @param projection
     */
    void setProjection(IProjection *projection);

    /**
     * @brief setNormalize
     * @param normalize
     */
    void setNormalize(bool normalize);

    /**
     * @brief showEndMessage
     * @param showEndMessage
     */
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
    IPrinter* m_printer;
    IProjection *m_projection;

    bool mshowEndMessage;
};

#endif // GRADIENT_H
