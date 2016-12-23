#ifndef BORDERPARABOLIC_H
#define BORDERPARABOLIC_H

#include <global.h>
#include <parabolicequation.h>

/**
 * @brief The BorderParabolic class
 * u_t = a^2 u_xx + f
 * u(x,t) = x^2 + t^2
 */

class MINIMUMSHARED_EXPORT BorderParabolic : public IParabolicEquation
{
public:
    virtual double initial(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    double hx = 0.01;
    double ht = 0.01;
    unsigned int N = 100;
    unsigned int M = 100;
    double a = 1.0;

    static void Main(int argc, char* argv[]);

    void calculateN41(DoubleMatrix &u);
    void calculateN42(DoubleMatrix &u);
    void calculateN6(DoubleMatrix &u);
};

#endif // BORDERPARABOLIC_H
