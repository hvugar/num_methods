#include "gradient_cs.h"
#include "printer.h"
#include "projection.h"
#include "function.h"
#include "vectornormalizer.h"
#include <math.h>

ConstStepGradient::ConstStepGradient() : GradientBasedMethod(),
    mx(nullptr), mg(nullptr)
{}

ConstStepGradient::~ConstStepGradient()
{}

void ConstStepGradient::calculate(DoubleVector &x)
{
    double alpha = 0.0;
    double f1 = 0.0;
    double f2 = 0.0;

    size_t n = x.length();
    DoubleVector g(n);

    mx = &x;
    mg = &g;
    m_iterationNumber = 0;
    m_functionEvaluationNumber = 0;

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
        if (m_printer) m_printer->print(m_iterationNumber, x, g, m_fn->fx(x), alpha, MethodResult::BREAK_FIRST_ITERATION);
        if (m_show_end_message) puts("Optimisation ends, because norm of gradient is less than optimality tolerance...");
        return;
    }

    /**************************************************************************************
     * Value of objective functionin on initial point.
     **************************************************************************************/
    f1 = m_fn->fx(x);

    if (m_printer != nullptr) m_printer->print(m_iterationNumber, x, g, f1, alpha, MethodResult::FIRST_ITERATION);

    do
    {
        m_iterationNumber++;

        /**************************************************************************************
         * Normalization of a gradient vector
         **************************************************************************************/
        if (m_normalize && m_normalizer != nullptr) m_normalizer->normalize(g);

        /**************************************************************************************
         * One-dimensional minimization along the direction of a anti-gradient
         **************************************************************************************/
        alpha = minimize(x, g);

        /**************************************************************************************
         * Calculation next point.
         **************************************************************************************/
        DoubleVector cx = x;
        for (size_t i=0; i<n; i++)
        {
            //double cx = x[i];
            x[i] = x[i] - alpha * g[i];

            if (m_projection != nullptr) m_projection->project(x, i);

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
        if (m_normalizer != nullptr) gradient_norm = m_normalizer->norm(g);
        if (gradient_norm < optimalityTolerance())
        {
            if (m_printer != nullptr) m_printer->print(m_iterationNumber, x, g, f2, alpha, MethodResult::BREAK_OPTIMALITY_TOLERANCE);
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
            if (m_printer != nullptr) m_printer->print(m_iterationNumber, x, g, f2, alpha, MethodResult::BREAK_STEP_TOLERANCE);
            if (m_show_end_message) puts("Optimisation ends, because distance between previous and current point less than step tolerance...");
            break;
        }
        /**************************************************************************************
         *
         *
         *
         **************************************************************************************/
        if (distance <= stepTolerance())
        {
            if (m_printer != nullptr) m_printer->print(m_iterationNumber, x, g, f2, alpha, MethodResult::BREAK_STEP_TOLERANCE);
            if (m_show_end_message) puts("Optimisation ends, because distance between previous and current point less than step tolerance...");
            break;
        }
        /**************************************************************************************
         *
         *
         *
         **************************************************************************************/
        if (fabs(f2 - f1) <= functionTolerance())
        {
            if (m_printer != nullptr) m_printer->print(m_iterationNumber, x, g, f2, alpha, MethodResult::BREAK_FUNCTION_TOLERANCE);
            if (m_show_end_message) puts("Optimisation ends, because previous and current function values difference less than function tolerance...");
            break;
        }
        /**************************************************************************************
         *
         *
         *
         **************************************************************************************/
        if (m_iterationNumber == maxIterationCount())
        {
            if (m_printer != nullptr) m_printer->print(m_iterationNumber, x, g, f2, alpha, MethodResult::BREAK_ITERATION_NUMBER);
            if (m_show_end_message) puts("Optimisation ends, because iteration count reached max allowed iterations number...");
            break;
        }

        /**************************************************************************************
         *
         *
         *
         **************************************************************************************/
        if (m_functionEvaluationNumber == maxFunctionEvaluationCount())
        {
            if (m_printer != nullptr) m_printer->print(m_iterationNumber, x, g, f2, alpha, MethodResult::BREAK_FUNCTION_EVALUATION_NUMBER);
            if (m_show_end_message) puts("Optimisation ends, because max function evaluation count reached max allowed number...");
            break;
        }

        /**************************************************************************************
         * Printing iteration information.
         **************************************************************************************/
        if (m_printer != nullptr) m_printer->print(m_iterationNumber, x, g, f2, alpha, MethodResult::NEXT_ITERATION);

        f1 = f2;

    } while (true);

    g.clear();
}

double ConstStepGradient::minimize(const DoubleVector &x, const DoubleVector &g) const
{
    size_t n = x.length();

    DoubleVector cx(n);

    double alpha = min_step;
    for (size_t i=0; i < n; i++)
    {
        cx[i] = x[i] - alpha * g[i];
        if (m_projection != nullptr) m_projection->project(cx, i);
    }
    if (m_projection != nullptr) m_projection->project(cx);

    while (m_fn->fx(cx) > m_fn->fx(x))
    {
        alpha = alpha * 0.5;
        for (size_t i=0; i < n; i++)
        {
            cx[i] = x[i] - alpha * g[i];
            if (m_projection != nullptr) m_projection->project(cx, i);
        }
        if (m_projection != nullptr) m_projection->project(cx);
    }
    return alpha;
}

double ConstStepGradient::fx(double alpha) const
{
    DoubleVector &x = *mx;
    DoubleVector &g = *mg;
    size_t n = x.length();

    DoubleVector cx = x;
    for (size_t i=0; i<n; i++)
    {
        cx[i] = x[i] - alpha * g[i];
        if (m_projection != nullptr) m_projection->project(cx, i);
    }

    if (m_projection != nullptr) m_projection->project(cx);
    double r_fx = m_fn->fx(cx);
    cx.clear();

    return r_fx;
}
