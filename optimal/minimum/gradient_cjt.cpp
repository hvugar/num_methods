#include "gradient_cjt.h"
#include "printer.h"
#include "projection.h"
#include "function.h"
#include <math.h>

ConjugateGradient::ConjugateGradient() : GradientMethod()
{
    setNormalize(true);
}

ConjugateGradient::~ConjugateGradient()
{}

void ConjugateGradient::calculate(DoubleVector& x)
{
    unsigned int n = x.size();

    unsigned int k = 0;
    double gr0_mod = 0.0;
    double gr1_mod = 0.0;
    double gr2_mod = 0.0;
    double distance = 0.0;
    double alpha = 0.0;
    double f1 = 0.0;
    double f2 = 0.0;
    bool firstIteration = true;

    DoubleVector g(n);
    DoubleVector s(n);

    mx = &x;
    ms = &s;
    iterationCount = 0;

    do
    {
        // Gradient of objectiv function in current point
        m_gr->gradient(x, g);

        double gradient_norm = g.L2Norm();
        if (gradient_norm < epsilon1())
        {
            if (mshowEndMessage) puts("Optimisation ends, because L2 norm of gradient is less than epsilon...");
            break;
        }

        iterationCount++;

        // Module of gradient
        gr0_mod = 0.0;
        for (unsigned int i=0; i<n; i++) gr0_mod = gr0_mod + (g[i]*g[i]);
        //gr0_mod = sqrt(gr0_mod);

        if (k == 0)
        {
            gr1_mod = gr0_mod;
            // First direction is antigradient
            for (unsigned int i=0; i<n; i++) s[i] = -g[i];
        }
        else
        {
            gr2_mod = gr0_mod;
            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;
            // Direction in next iteration (k != 0)
            for (unsigned int i=0; i<n; i++) s[i] = -g[i] + s[i] * w;
        }

        /* Normalize vector */
        if (m_normalize) s.L2Normalize();

        /* R1 minimization in direct of antigradient */
        alpha = minimize(x, s);

        f1 = f2;
        if (firstIteration)
        {
            f1 = m_fn->fx(x);
            firstIteration = false;
        }

        distance = 0.0;
        for (unsigned int i=0; i<n; i++)
        {
            double cx = x[i];
            x[i] = x[i] + alpha * s[i];

            if (m_projection != NULL) m_projection->project(x, i);

            distance += (x[i]-cx)*(x[i]-cx);
        }
        distance = sqrt(distance);
        f2 = m_fn->fx(x);

        if ( k == x.size() ) { k = 0; } else { k++; }

        if (m_printer != NULL) m_printer->print(iterationCount, x, g, alpha, m_fn);

        /* calculating distance previous and new point */
        if (distance < epsilon2() && fabs(f2 - f1) < epsilon3())
        {
            if (mshowEndMessage) puts("Optimisation ends, because distance between last and current point less than epsilon...");
            break;
        }

    } while (true);
}

double ConjugateGradient::minimize(const DoubleVector &x, const DoubleVector &s)
{
    C_UNUSED(x);
    C_UNUSED(s);

    double alpha0 = 1.0;
    double a,b,alpha;

    stranghLineSearch(alpha0, min_step, a, b, this);
    goldenSectionSearch(a, b, alpha, this, min_epsilon);

    if (this->fx(alpha) > this->fx(alpha0)) alpha = alpha0;
    return alpha;
}

double ConjugateGradient::fx(double alpha)
{
    DoubleVector &x = *mx;
    DoubleVector &s = *ms;
    unsigned int n = x.size();

    DoubleVector cx(n);
    for (unsigned int i=0; i<n; i++)
    {
        cx[i] = x[i] + alpha * s[i];
        if (m_projection != NULL) m_projection->project(cx, i);
    }

    return m_fn->fx(cx);
}
