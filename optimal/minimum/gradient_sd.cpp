#include "gradient_sd.h"
#include "printer.h"
#include "projection.h"
#include "function.h"
#include <math.h>

SteepestDescentGradient::SteepestDescentGradient() : GradientMethod()
{
    setNormalize(true);
    r1m.setFunction(this);
}

SteepestDescentGradient::~SteepestDescentGradient()
{}

void SteepestDescentGradient::calculate(DoubleVector &x)
{
    double distance = 0.0;
    double alpha = 0.0;
    double f1 = 0.0;
    double f2 = 0.0;

    unsigned int n = x.length();
    DoubleVector g(n);

    mx = &x;
    mg = &g;
    iterationCount = 0;

    // Gradient of objectiv function in current point
    m_gr->gradient(x, g);
    f1 = m_fn->fx(x);

    /* checking gradient norm */
    /* if gradient norm at current point is less than epsilon then return. no minimize */
    double gradient_norm = g.L2Norm();
    if (gradient_norm < epsilon1())
    {
        if (m_printer != NULL) m_printer->print(iterationCount, x, g, f1, BREAK_FIRST_ITERATION);
        if (mshowEndMessage) puts("Optimisation ends, because norm of gradient is less than epsilon...");
        return;
    }
    if (m_printer != NULL) m_printer->print(iterationCount, x, g, f1, FIRST_ITERATION);

    do
    {

        iterationCount++;

        /* Normalize vector */
        if (m_normalize) g.L2Normalize();

        /* R1 minimization in direct of antigradient */
        alpha = minimize(x, g);
        printf("%f %d\n", alpha, n);

        distance = 0.0;
        double f1 = m_fn->fx(x);
        for (unsigned int i=0; i < n; i++)
        {
            double cx = x[i];
            x[i] = x[i] - alpha * g[i];

            if (m_projection != NULL) m_projection->project(x, i);

            // calculating distance
            distance += (x[i]-cx)*(x[i]-cx);
        }
        distance = sqrt(distance);
        f2 = m_fn->fx(x);

        // Gradient of objectiv function in current point
        m_gr->gradient(x, g);

        /* checking gradient norm */
        /* if gradient norm at current point is less than epsilon then return. no minimize */
        double gradient_norm = g.L2Norm();
        if (gradient_norm < epsilon1())
        {
            if (m_printer != NULL) m_printer->print(iterationCount, x, g, f2, GradientMethod::BREAK_GRADIENT_NORM_LESS);
            if (mshowEndMessage) puts("Optimisation ends, because norm of gradient is less than epsilon...");
            break;
        }

        /* calculating distance previous and new point */
        if (distance < epsilon2() && fabs(f2 - f1) < epsilon3())
        {
            if (m_printer != NULL) m_printer->print(iterationCount, x, g, f2, GradientMethod::BREAK_DISTANCE_LESS);
            if (mshowEndMessage) puts("Optimisation ends, because distance between last and current point less than epsilon...");
            break;
        }

        if (m_printer != NULL) m_printer->print(iterationCount, x, g, f2, GradientMethod::NEXT_ITERATION);

        f1 = f2;

    } while (true);

    g.clear();
}

double SteepestDescentGradient::minimize(const DoubleVector &x, const DoubleVector &g) const
{
    C_UNUSED(x);
    C_UNUSED(g);

    double alpha0 = 0.0;
    double a,b,alpha;

    //stranghLineSearch(alpha0, min_step, a, b, this);
    //goldenSectionSearch(a, b, alpha, this, min_epsilon);
    //if (fx(alpha) > fx(alpha0)) alpha = alpha0;

    double fxa, fxb;
    r1m.straightLineSearch(alpha0, min_step, a, b, fxa, fxb);
    r1m.goldenSectionSearch(alpha, a, b, min_epsilon);
    if (fx(alpha) > fx(alpha0)) alpha = alpha0;

    return alpha;
}

double SteepestDescentGradient::fx(double alpha) const
{
    DoubleVector &x = *mx;
    DoubleVector &g = *mg;
    unsigned int n = x.length();

    DoubleVector cx(n);
    for (unsigned int i=0; i<n; i++)
    {
        cx[i] = x[i] - alpha * g[i];
        if (m_projection != NULL) m_projection->project(cx, i);
    }

    return m_fn->fx(cx);
}

R1FxMinimizer &SteepestDescentGradient::r1Minimizer()
{
    return r1m;
}

const R1FxMinimizer& SteepestDescentGradient::r1Minimizer() const
{
    return r1m;
}
