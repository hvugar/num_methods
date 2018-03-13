#ifndef IGRADIENT_H
#define IGRADIENT_H

#include "global.h"

#include "../r1minimize.h"
#include "../function.h"

class RnFunction;
struct R1Function;
class IGradient;
class IPrinter;
class IProjection;

/**
 * @brief The Gradient Method Interface
 */
class MINIMUMSHARED_EXPORT IGradientMethod : protected RnFunction, protected IGradient, protected R1Function
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
    IGradientMethod();
    virtual ~IGradientMethod();

    virtual double fx(double x) const = 0;
    virtual double fx(const DoubleVector &x) const = 0;
    virtual void gradient(const DoubleVector &x, DoubleVector &g) const = 0;
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, MethodResult result) const = 0;
    virtual void project(DoubleVector &x, int index) const = 0;


    virtual void calculate(DoubleVector &x) = 0;

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

    double m_epsilon1;
    double m_epsilon2;
    double m_epsilon3;
    double min_epsilon;
    double min_step;
    int iterationCount;
    bool m_normalize;
    bool mshowEndMessage;
};

class MINIMUMSHARED_EXPORT IConjugateGradient : protected IGradientMethod
{
public:
    IConjugateGradient();
    virtual ~IConjugateGradient();

    virtual double fx(double alpha) const = 0;
    virtual double fx(const DoubleVector &x) const = 0;
    virtual void gradient(const DoubleVector &x, DoubleVector &g) const = 0;
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, MethodResult result) const;
    virtual void project(DoubleVector &x, int index) const;

    virtual void calculate(DoubleVector &x);
protected:
    virtual double minimize(const DoubleVector &x, const DoubleVector &g);

    DoubleVector *mx;
    DoubleVector *ms;
};

#endif // IGRADIENT_H
