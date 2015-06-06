#include "fpgradient.h"

FastProximalGradient::FastProximalGradient() : Gradient()
{
    mf1 = new FastProximalGradient::ArgMin(mx, mg, this);
}

FastProximalGradient::~FastProximalGradient()
{}

double FastProximalGradient::minimize()
{
    R1Minimize r1;
    r1.setF(mf1);
    r1.setX0(0.0);
    r1.setStep(min_step);
    r1.setEpsilon(min_epsilon);
    r1.straightLineSearch();
    double alpha = r1.goldenSectionSearch();
    return alpha;
}

void FastProximalGradient::calculate()
{
    mcount = 0;
    do
    {
        /* calculating function gradient at current point */
        calcGradient();

        /* if norm of gradinet at current point is less than epsilon break. no minimize */
        if (gradientNorm() < epsilon())
            break;

        mcount++;

        /* R1 minimization in direct of antigradient */
        malpha = minimize();

        print();

        for (unsigned int i=0; i<mx.size(); i++)
        {
            mx[i] = mx[i] - malpha * mg[i];
        }

        /* calculating distance previous and new point */
    } while (distance() > epsilon());
}

void FastProximalGradient::print()
{
    printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", mx[0], mx[1], mg[0], mg[1], malpha, gradientNorm(), f()->fx(mx));
}

FastProximalGradient::ArgMin::ArgMin(std::vector<double> &x, std::vector<double> &g, Gradient *gradient) :
    x(x), g(g), gradient(gradient)
{}

double FastProximalGradient::ArgMin::fx(double alpha)
{
    std::vector<double> x1 = x;
    for (unsigned int i=0; i < x.size(); i++)
    {
        x1[i] = x[i] - alpha * g[i];
    }
    return gradient->f()->fx(x1);
}


