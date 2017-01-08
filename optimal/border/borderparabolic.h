#ifndef BORDERPARABOLIC_H
#define BORDERPARABOLIC_H

#include <global.h>
#include <parabolicequation.h>

/**
 * @brief The BorderParabolic class
 * u_t = a^2 u_xx + f
 * u(x,t) = x^2 + t^2
 */

#define SAMPLE_3

class MINIMUMSHARED_EXPORT BorderParabolic : public IParabolicEquation
{
public:
    virtual double initial(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    double u(unsigned int i, unsigned int j) const;

    void calculateN4L2RM(DoubleMatrix &u);
    void calculateN4R2LM(DoubleMatrix &u);

    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
    double a = 1.0;

    static void Main(int argc, char* argv[]);
};

#endif // BORDERPARABOLIC_H
