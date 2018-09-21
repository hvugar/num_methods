#include "gradient_sd.h"
#include "printer.h"
#include "projection.h"
#include "function.h"
#include <math.h>

SteepestDescentGradient::SteepestDescentGradient() : GradientMethod()
{
    setNormalize(true);
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
    m_iteration_count = 0;

    /**************************************************************************************
     * Gradient of objective functionin on initial point.
     **************************************************************************************/
    m_gr->gradient(x, g);

    /**************************************************************************************
     * Checking for gradient vector norm that is less of epsilon.
     * If gradient vector norm in current point is less than epsilon then break the iteration.
     * Finish minimization.
     **************************************************************************************/
    double gradient_norm = g.L2Norm();
    if (gradient_norm < epsilon1())
    {
        if (m_printer != NULL) m_printer->print(m_iteration_count, x, g, m_fn->fx(x), alpha, BREAK_FIRST_ITERATION);
        if (m_show_end_message) puts("Optimisation ends, because norm of gradient is less than epsilon...");
        return;
    }

    /**************************************************************************************
     * Value of objective functionin on initial point.
     **************************************************************************************/
    f1 = m_fn->fx(x);

    if (m_printer != NULL) m_printer->print(m_iteration_count, x, g, f1, alpha, FIRST_ITERATION);

    do
    {
        m_iteration_count++;

        /**************************************************************************************
         * Normalization of a vector
         **************************************************************************************/
        if (m_normalize) g.L2Normalize();

        /**************************************************************************************
         * One-dimensional minimization along the direction of a anti-gradient
         **************************************************************************************/
        alpha = minimize(x, g);

        /**************************************************************************************
         * Calculation next point.
         **************************************************************************************/
        distance = 0.0;
        for (unsigned int i=0; i < n; i++)
        {
            double cx = x[i];
            x[i] = x[i] - alpha * g[i];

            if (m_projection != NULL) m_projection->project(x, i);

            /**************************************************************************************
             * Calculating distance.
             **************************************************************************************/
            distance += (x[i]-cx)*(x[i]-cx);
        }
        if (m_projection != NULL) m_projection->project(x);
        distance = sqrt(distance);
        f2 = m_fn->fx(x);

        /**************************************************************************************
         * Gradient of objectiv function in next point
         **************************************************************************************/
        m_gr->gradient(x, g);

        /**************************************************************************************
         * Checking for gradient vector norm that is less of epsilon.
         * If gradient vector norm in current point is less than epsilon then break the iteration.
         * Finish minimization.
         **************************************************************************************/
        double gradient_norm = g.L2Norm();
        if (gradient_norm < epsilon1())
        {
            if (m_printer != NULL) m_printer->print(m_iteration_count, x, g, f2, alpha, GradientMethod::BREAK_GRADIENT_NORM_LESS);
            if (m_show_end_message) puts("Optimisation ends, because norm of gradient is less than epsilon...");
            break;
        }

        /**************************************************************************************
         * Calculating distance between the previous and current points.
         * Calculating difference values of functions in previous and current points.
         * If distance and difference is less than epsilon then break the iteration.
         * Finish minimization.
         **************************************************************************************/
        if (distance < epsilon2() && fabs(f2 - f1) < epsilon3())
        {
            if (m_printer != NULL) m_printer->print(m_iteration_count, x, g, f2, alpha, GradientMethod::BREAK_DISTANCE_LESS);
            if (m_show_end_message) puts("Optimisation ends, because distance between last and current point less than epsilon...");
            break;
        }

        /**************************************************************************************
         * Printing iteration information.
         **************************************************************************************/
        if (m_printer != NULL) m_printer->print(m_iteration_count, x, g, f2, alpha, GradientMethod::NEXT_ITERATION);

        f1 = f2;

    } while (true);

    g.clear();
}

double SteepestDescentGradient::minimize(const DoubleVector &x, const DoubleVector &g) const
{
    C_UNUSED(x);
    C_UNUSED(g);

    double alpha0 = min_step;
    double a,b,alpha;

    double fxa, fxb;
    bool unimodal;

    straightLineSearch(alpha0, min_step, a, b, fxa, fxb, unimodal);
    //swann(alpha0, min_step, a, b, fxa, fxb, unimodal);
    printf("unimodal %d\n", unimodal);

    if (unimodal)
    {
        goldenSectionSearch(alpha, a, b, min_epsilon);
    }
    else
    {
        fxa < fxb ? alpha = a : alpha = b;
    }

    if (fx(alpha) > fx(alpha0)) alpha = alpha0;

    return alpha;
}

double SteepestDescentGradient::fx(double alpha) const
{
    DoubleVector &x = *mx;
    DoubleVector &g = *mg;
    unsigned int n = x.length();

    DoubleVector cx = x;
    for (unsigned int i=0; i<n; i++)
    {
        cx[i] = x[i] - alpha * g[i];
        if (m_projection != NULL) m_projection->project(cx, i);
    }

    if (m_projection != NULL) m_projection->project(cx);

    return m_fn->fx(cx);
}
