#include "cjtgradient.h"

ConjugateGradient::ConjugateGradient() : Gradient()
{
    mf1 = new ConjugateGradient::ArgMin(mx, s, this);
}

ConjugateGradient::~ConjugateGradient()
{}

double ConjugateGradient::minimize()
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

void ConjugateGradient::calculate()
{
    unsigned int k = 0;
    double gr1_mod = 0.0;
    double gr2_mod = 0.0;
    s = mx;

    do
    {
        // First iteration
        if (k == 0)
        {
            // Gradient of objectiv function in current point
            calcGradient();

            // First direction is antigradient
            for (unsigned int i=0; i<mg.size(); i++) s[i] = -mg[i];

            // Module of gradient
            gr1_mod = 0.0;
            for (unsigned int i=0; i<mg.size(); i++) gr1_mod = gr1_mod + mg[i]*mg[i];
        }
        else
        {
            // Gradient of objectiv function in next point
            calcGradient();

            // Module of next gradient
            gr2_mod = 0.0;
            for (unsigned int i=0; i<mg.size(); i++) gr2_mod = gr2_mod + mg[i]*mg[i];

            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;

            // Direction in next (k+1) iteration
            for (unsigned int i=0; i<mg.size(); i++) s[i] = -mg[i] + s[i] * w;
        }

        if (gradientNorm() < epsilon())
            break;

        mcount++;

        malpha = minimize();

        print();

        for (unsigned int i=0; i<mx.size(); i++)
        {
            mx[i] = mx[i] + malpha * s[i];
        }

        if ( k == mx.size() ) { k = 0; } else { k++; }

        /* calculating distance previous and new point */
    } while (distance() > epsilon());
}

void ConjugateGradient::print()
{
    printf("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", mx[0], mx[1], mg[0], mg[1], s[0], s[1],  malpha, gradientNorm(), f()->fx(mx));
}

ConjugateGradient::ArgMin::ArgMin(std::vector<double> &x, std::vector<double> &g, Gradient *gradient) :
    x(x), g(g), gradient(gradient)
{}

double ConjugateGradient::ArgMin::fx(double alpha)
{
    std::vector<double> x1 = x;
    for (unsigned int i=0; i < x.size(); i++)
    {
        x1[i] = x[i] + alpha * g[i];
    }
    return gradient->f()->fx(x1);
}
