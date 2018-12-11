#include "gradient_cjt.h"
#include "printer.h"
#include "projection.h"
#include "function.h"
#include "vectornormalizer.h"
#include <math.h>

ConjugateGradient::ConjugateGradient() : GradientMethod()
{
    setNormalize(true);
    setAlgorithm(FLETCHER_REEVES);
    setResetIteration(true);
}

ConjugateGradient::~ConjugateGradient()
{}

void ConjugateGradient::calculate(DoubleVector& x)
{
    //double distance = 0.0;
    double alpha = 0.0;
    double f1 = 0.0;
    double f2 = 0.0;
    unsigned int k = 0;
    double grad_mod_cur = 0.0;
    double grad_mod_prv = 0.0;

    unsigned int n = x.length();
    DoubleVector g(n);
    DoubleVector s(n);

    mx = &x;
    ms = &s;
    m_iterationCount = 0;
    //unsigned int m_funcEvaluationCount = 0;

    /**************************************************************************************
     * Gradient of objective functionin on initial point.
     **************************************************************************************/
    m_gr->gradient(x, g);

    /**************************************************************************************
     * Checking for gradient vector norm that is less of optimality tolerance.
     * If gradient vector norm in current point is less than optimality tolerance then break
     * the iteration.
     * Finish minimization.
     **************************************************************************************/
    double gradient_norm = 0.0;
    if (m_normalizer != nullptr) gradient_norm = m_normalizer->norm(g);
    if (gradient_norm < optimalityTolerance())
    {
        if (m_printer) m_printer->print(m_iterationCount, x, g, m_fn->fx(x), alpha, GradientMethod::BREAK_FIRST_ITERATION);
        if (m_show_end_message) puts("Optimisation ends, because norm of gradient is less than optimality tolerance...");
        return;
    }

    /**************************************************************************************
     * Value of objective functionin on initial point.
     **************************************************************************************/
    f1 = m_fn->fx(x);

    if (m_printer) m_printer->print(m_iterationCount, x, g, f1, alpha, GradientMethod::FIRST_ITERATION);

    do
    {
        m_iterationCount++;

        if (malgoritm == FLETCHER_REEVES)
        {
            // Module of gradient
            grad_mod_cur = 0.0;
            for (unsigned int i=0; i<n; i++) grad_mod_cur += g[i]*g[i];

            if (k == 0)
            {
                // First direction is antigradient
                grad_mod_prv = grad_mod_cur;
                for (unsigned int i=0; i<n; i++) s[i] = -g[i];
            }
            else
            {
                // Direction in next iteration (k != 0)
                double w = grad_mod_cur / grad_mod_prv;
                grad_mod_prv = grad_mod_cur;
                for (unsigned int i=0; i<n; i++) s[i] = -g[i] + s[i] * w;
            }
        }

        if (malgoritm == POLAK_RIBIERE)
        {
            if (k == 0)
            {
                for (unsigned int i=0; i<n; i++) s[i] = -g[i];
            }
            else
            {
                double w = 1.0;
                for (unsigned int i=0; i<n; i++) s[i] = -g[i] + s[i] * w;
            }
        }

        /**************************************************************************************
         * Reset the direction.
         **************************************************************************************/
        if (mResetIteration)
        {
            if ( k == n ) { k = 0; } else { k++; }
        }
        else
        {
            k++;
        }

        /**************************************************************************************
         * Normalization of a gradient vector
         **************************************************************************************/
        if (m_normalize && m_normalizer != nullptr) m_normalizer->normalize(s);

        /**************************************************************************************
         * One-dimensional minimization along the direction of a anti-gradient
         **************************************************************************************/
        alpha = minimize(x, s);

        /**************************************************************************************
         * Calculation next point.
         **************************************************************************************/
        DoubleVector cx = x;
        for (unsigned int i=0; i<n; i++)
        {
            //double cx = x[i];
            x[i] = x[i] + alpha * s[i];

            //if (m_projection) m_projection->project(x, i);

            /**************************************************************************************
             * Calculating distance.
             **************************************************************************************/
            //distance += (x[i]-cx)*(x[i]-cx);
        }
        if (m_projection != nullptr) m_projection->project(x);

        double distance = 0.0;
        for (unsigned int i=0; i<n; i++) { distance += (cx[i]-x[i])*(cx[i]-x[i]); }
        distance = sqrt(distance);
        cx.clear();

        f2 = m_fn->fx(x);

        /**************************************************************************************
         * Gradient of objectiv function in next point
         **************************************************************************************/
        m_gr->gradient(x, g);

        /**************************************************************************************
         * Checking for gradient vector norm that is less of optimality tolerance.
         * If gradient vector norm in current point is less than optimality tolerance then break
         * the iteration.
         * Finish minimization.
         **************************************************************************************/
        double gradient_norm = 0.0;
        if (m_normalizer) gradient_norm = m_normalizer->norm(g);
        if (gradient_norm < optimalityTolerance())
        {
            if (m_printer) m_printer->print(m_iterationCount, x, g, f2, alpha, GradientMethod::BREAK_GRADIENT_NORM_LESS);
            if (m_show_end_message) puts("Optimisation ends, because norm of gradient is less than optimality tolerance...");
            break;
        }

        /**************************************************************************************
         * Calculating distance between the previous and current points.
         * Calculating difference values of functions in previous and current points.
         * If distance and difference is less than step tolerance then break the iteration.
         * Finish minimization.
         **************************************************************************************/
        if (distance < stepTolerance() && fabs(f2 - f1) < functionTolerance())
        {
            if (m_printer != nullptr) m_printer->print(m_iterationCount, x, g, f2, alpha, GradientMethod::BREAK_DISTANCE_LESS);
            if (m_show_end_message) puts("Optimisation ends, because distance between previous and current point less than step tolerance...");
            break;
        }

        /**************************************************************************************
         *
         *
         *
         **************************************************************************************/
        if (m_iterationCount >= maxIterations())
        {
            if (m_printer != nullptr) m_printer->print(m_iterationCount, x, g, f2, alpha, GradientMethod::BREAK_DISTANCE_LESS);
            if (m_show_end_message) puts("Optimisation ends, because iteration count reached max allowed iterations number...");
            break;
        }

        /**************************************************************************************
         *
         *
         *
         **************************************************************************************/
        if (maxFunctionEvaluations() != maxFunctionEvaluations())
        {
            if (m_printer != nullptr) m_printer->print(m_iterationCount, x, g, f2, alpha, GradientMethod::BREAK_DISTANCE_LESS);
            if (m_show_end_message) puts("Optimisation ends, because max function evaluation count reached max allowed number...");
            break;
        }

        /**************************************************************************************
         * Printing iteration information.
         **************************************************************************************/
        if (m_printer != nullptr) m_printer->print(m_iterationCount, x, g, f2, alpha, GradientMethod::NEXT_ITERATION);

        f1 = f2;

    } while (true);

    g.clear();
    s.clear();
}

double ConjugateGradient::minimize(const DoubleVector &x, const DoubleVector &s) const
{
    C_UNUSED(x);
    C_UNUSED(s);

    double alpha0 = min_step;
    double a=0.0, b=0.0, alpha=0.0;

    double fxa = 0.0, fxb = 0.0;
    bool unimodal = false;

    straightLineSearch(alpha0, min_step, a, b, fxa, fxb, unimodal);
    //swann(alpha0, min_step, a, b, fxa, fxb, unimodal);

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

double ConjugateGradient::fx(double alpha) const
{
    const DoubleVector &x = *mx;
    const DoubleVector &s = *ms;
    unsigned int n = x.length();

    DoubleVector cx = x;
    for (unsigned int i=0; i<n; i++)
    {
        cx[i] = x[i] + alpha * s[i];
        //if (m_projection) m_projection->project(cx, i);
    }

    if (m_projection != nullptr) m_projection->project(cx);
    double r_fx = m_fn->fx(cx);
    cx.clear();

    return r_fx;
}

void ConjugateGradient::setAlgorithm(Algorithm algorithm)
{
    malgoritm = algorithm;
}

void ConjugateGradient::setResetIteration(bool reset)
{
    mResetIteration = reset;
}
