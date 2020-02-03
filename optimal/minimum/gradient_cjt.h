#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "gradient.h"

/**
 * @brief
 * En: Method of Conjugate Gradient.
 * Ru: Метод сопряженных градиентов. (Метод Флетчера-Ривса)
 * Az: Qoşma qradiyentlər üsulu.
 */
class MINIMUMSHARED_EXPORT ConjugateGradient : public GradientMethod
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
    virtual void calculate1(DoubleVector &x);

    void setAlgorithm(Algorithm algorithm);
    void setResetIteration(bool reset);

protected:
    virtual double minimize(const DoubleVector &x, const DoubleVector &g) const;
    virtual double fx(double alpha) const;

private:
    DoubleVector *mx;
    DoubleVector *ms;
    Algorithm m_algoritm;
    bool m_reset_iteration;
};

#endif // CONJUGATE_GRADIENT_H
