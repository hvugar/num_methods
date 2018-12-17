#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

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
    enum Algorithm
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

    virtual void callBeforeFindingLine() const {}
    virtual void callAftreFindingLine() const {}
    virtual void callBeforeGoldenSearch() const {}
    virtual void callAfterGoldenSearch() const {}

    DoubleVector *mx;
    DoubleVector *ms;

private:
    Algorithm malgoritm;
    bool mResetIteration;
};

#endif // CONJUGATEGRADIENT_H
