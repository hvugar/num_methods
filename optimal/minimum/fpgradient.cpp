#include "fpgradient.h"

FastProximalGradient::FastProximalGradient() : Gradient()
{
    min_step = 0.1;
    min_epsilon = 0.000001;
    grad_step = 0.000001;
    mepsilon = 0.000001;
}

FastProximalGradient::~FastProximalGradient()
{}

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


