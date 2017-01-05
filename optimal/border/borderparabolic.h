#ifndef BORDERPARABOLIC_H
#define BORDERPARABOLIC_H

#include <global.h>
#include <parabolicequation.h>

/**
 * @brief The BorderParabolic class
 * u_t = a^2 u_xx + f
 * u(x,t) = x^2 + t^2
 */

#define SAMPLE_1

class MINIMUMSHARED_EXPORT BorderParabolic : public IParabolicEquation
{
public:
    virtual double initial(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    double u(unsigned int i, unsigned int j) const;
    void calculateN4(DoubleMatrix &u);

    void calculateN4L2R(DoubleMatrix &u);
    void calculateN4R2L(DoubleMatrix &u);

    double hx = 0.0001;
    double ht = 0.01;
    unsigned int N = 10000;
    unsigned int M = 100;
    double a = 1.0;

    static void Main(int argc, char* argv[]);

//    void calculateN41(DoubleMatrix &u);
//    void calculateN42(DoubleMatrix &u);
//    void calculateN43(DoubleMatrix &u);
//    void calculateN44(DoubleMatrix &u);
//    void calculateN45(DoubleMatrix &u);
};

#endif // BORDERPARABOLIC_H
