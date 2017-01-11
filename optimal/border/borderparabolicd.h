#ifndef BORDERPARABOLICD_H
#define BORDERPARABOLICD_H

#include <global.h>
#include <parabolicequation.h>

/**
 * @brief The BorderParabolicD class
 * u_t = a^2 u_xx + f
 * u(x,t) = x^2 + t^2
 */

#define SAMPLE_4

class MINIMUMSHARED_EXPORT BorderParabolicD : public IParabolicEquation
{
public:
    virtual double initial(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    double U(unsigned int i, unsigned int j) const;

    void calculateN4L2RM(DoubleMatrix &u);
    void calculateN4R2LM(DoubleMatrix &u);

    void calculateN6L2RM(DoubleMatrix &u);
    void calculateN6R2LM(DoubleMatrix &u);

    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
    double a = 1.0;

    static void Main(int argc, char* argv[]);
};

#endif // BORDERPARABOLICD_H
