#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "gradient.h"

/**
 * @brief
 * En: Method of Conjugate Gradient.
 * Ru: Метод сопряженных градиентов. (Метод Флетчера-Ривса)
 * Az: Qoşma qradiyentlər üsulu.
 */
class MINIMUMSHARED_EXPORT ConjugateGradient : public GradientBasedMethod
{
public:
    enum class Algorithm
    {
        FLETCHER_REEVES = 0,
        POLAK_RIBIERE = 1,
        HESTENES_STIEFEL = 2,
        DAI_YUAN = 3
    };

public:
    ConjugateGradient();
    virtual ~ConjugateGradient();

    virtual void calculate(DoubleVector &x);

    void setAlgorithm(Algorithm algorithm);
    void setResetIteration(bool reset);

protected:
    virtual double minimize(const DoubleVector &x, const DoubleVector &g) const;
    virtual double fx(double alpha) const;

protected:
    DoubleVector *mx;
    DoubleVector *ms;
    Algorithm m_algoritm;
    bool m_reset_iteration;
};


/**
 * @brief Метод сопряжённых градиентов (для решения СЛАУ)
 * Метод сопряженных градиентов — численный метод решения систем линейных алгебраических уравнений,
 * является итерационным методом Крыловского типа.
 */
class MINIMUMSHARED_EXPORT ConjugateGradientSLE  : public ConjugateGradient
{
public:
    ConjugateGradientSLE();
    virtual ~ConjugateGradientSLE();

    virtual void calculate(DoubleVector &x);

protected:
    virtual double minimize(const DoubleVector &x, const DoubleVector &g) const;

    DoubleVector *mg;
};

#endif // CONJUGATE_GRADIENT_H
