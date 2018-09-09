#ifndef GRADIENT_H
#define GRADIENT_H

#include "global.h"

#include "r1minimize.h"
#include "exceptions.h"

class RnFunction;
class IGradient;
class IPrinter;
class IProjection;
class IGradientPrinter;

/**
 * @brief The Abstract Gradient Method class
 */
class MINIMUMSHARED_EXPORT GradientMethod
{
public:
    enum MethodResult
    {
        NEXT_ITERATION = 0,
        FIRST_ITERATION = 1,
        BREAK_FIRST_ITERATION = 2,
        BREAK_GRADIENT_NORM_LESS = 3,
        BREAK_DISTANCE_LESS = 4
    };

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
    virtual double minimize(const DoubleVector &x, const DoubleVector &g) const = 0;

    /**
     * @brief Этап установления границ интервала
     */
    //void straightLineSearch(double x, double step, double &a, double &b, double &fxa, double &fxb, bool &unimodal) const;
    //void swann(double x, double step, double &a, double &b, double &fx, double &fxb, bool &unimodal) const;

    RnFunction *m_fn;
    IGradient *m_gr;
    IPrinter* m_printer;
    IGradientPrinter *m_printer_gr;
    IProjection *m_projection;
    double m_epsilon1;
    double m_epsilon2;
    double m_epsilon3;
    double min_step;
    double min_epsilon;
    int m_iteration_count;
    bool m_normalize;
    bool m_show_end_message;
};

struct MINIMUMSHARED_EXPORT IGradientPrinter
{
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double fxResult, GradientMethod::MethodResult result, double alpha);
};

#endif // GRADIENT_H
