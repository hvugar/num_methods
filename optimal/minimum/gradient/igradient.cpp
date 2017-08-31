#include "igradient.h"
#include "../printer.h"
#include "../projection.h"
#include "../function.h"
#include <math.h>

IConjugateGradient::IConjugateGradient() : IGradientMethod()
{
    setNormalize(true);
}

IConjugateGradient::~IConjugateGradient()
{}

void IConjugateGradient::calculate(DoubleVector& x)
{
    unsigned int k = 0;
    double gr0_mod = 0.0;
    double gr1_mod = 0.0;
    double gr2_mod = 0.0;
    double distance = 0.0;
    double alpha = 0.0;
    double f1 = 0.0;
    double f2 = 0.0;

    unsigned int n = x.length();
    DoubleVector g(n);
    DoubleVector s(n);

    mx = &x;
    ms = &s;
    iterationCount = 0;

    /**************************************************************************************
     * Gradient of objective functionin initial point.
     **************************************************************************************/
    gradient(x, g);

    /* checking gradient norm */
    /* if gradient norm at current point is less than epsilon then return. no minimize */
    /**************************************************************************************
     * checking gradient norm.
     * if gradient norm at current point is less than epsilon then return. no minimize
     **************************************************************************************/
    double gradient_norm = g.L2Norm();
    if (gradient_norm < epsilon1())
    {
        print(iterationCount, x, g, fx(x), BREAK_FIRST_ITERATION);
        if (mshowEndMessage) puts("Optimisation ends, because norm of gradient is less than epsilon...");
        return;
    }

    f1 = fx(x);
    print(iterationCount, x, g, f1, FIRST_ITERATION);

    do
    {
        iterationCount++;

        // Module of gradient
        gr0_mod = 0.0;
        for (unsigned int i=0; i<n; i++) gr0_mod = gr0_mod + (g[i]*g[i]);

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

        /**************************************************************************************
         * Normalization of a vector
         **************************************************************************************/
        if (m_normalize) s.L2Normalize();

        /**************************************************************************************
         * One-dimensional minimization along the direction of a anti-gradient
         **************************************************************************************/
        alpha = minimize(x, s);

        /**************************************************************************************
         * Calculation current point.
         **************************************************************************************/
        distance = 0.0;
        for (unsigned int i=0; i<n; i++)
        {
            double cx = x[i];
            x[i] = x[i] + alpha * s[i];

            project(x, i);

            distance += (x[i]-cx)*(x[i]-cx);
        }
        distance = sqrt(distance);
        f2 = fx(x);

        if ( k == x.length() ) { k = 0; } else { k++; }

        /**************************************************************************************
         * Gradient of objectiv function in current point
         **************************************************************************************/
        gradient(x, g);

        /**************************************************************************************
         * checking gradient norm.
         * if gradient norm at current point is less than epsilon then return. no minimize
         *
         **************************************************************************************/
        double gradient_norm = g.L2Norm();
        if (gradient_norm < epsilon1())
        {
            print(iterationCount, x, g, f2, BREAK_GRADIENT_NORM_LESS);
            if (mshowEndMessage) puts("Optimisation ends, because norm of gradient is less than epsilon...");
            break;
        }

        /**************************************************************************************
         * Calculating Euclied distance between the previous and current points.
         * Calculating difference values of functions on the previous and current points.
         **************************************************************************************/
        if (distance < epsilon2() && fabs(f2 - f1) < epsilon3())
        {
            print(iterationCount, x, g, f2, BREAK_DISTANCE_LESS);
            if (mshowEndMessage) puts("Optimisation ends, because distance between last and current point less than epsilon...");
            break;
        }

        f1 = f2;

        /**************************************************************************************
         * Printing iteration information.
         **************************************************************************************/
        print(iterationCount, x, g, f2, NEXT_ITERATION);
    } while (true);

    g.clear();
    s.clear();
}

double IConjugateGradient::minimize(const DoubleVector &x, const DoubleVector &s)
{
    C_UNUSED(x);
    C_UNUSED(s);

    double alpha0 = 0.0;
    double a,b,alpha;

    stranghLineSearch(alpha0, min_step, a, b, this);
    goldenSectionSearch(a, b, alpha, this, min_epsilon);
    if (fx(alpha) > fx(alpha0)) alpha = alpha0;

    return alpha;
}

double IConjugateGradient::fx(double alpha) const
{
    DoubleVector &x = *mx;
    DoubleVector &s = *ms;
    unsigned int n = x.length();

    DoubleVector cx(n);
    for (unsigned int i=0; i<n; i++)
    {
        cx[i] = x[i] + alpha * s[i];
        project(cx, i);
    }

    return fx(cx);
}

void IConjugateGradient::project(DoubleVector &, int) const {}

void IConjugateGradient::print(unsigned int, const DoubleVector &, const DoubleVector &, double, MethodResult) const {}

IGradientMethod::IGradientMethod() {}

IGradientMethod::~IGradientMethod() {}

/**
 * @brief Epsilon for gradient norm
 * @return
 */
double IGradientMethod::epsilon1() const
{
    return m_epsilon1;
}

/**
 * @brief Epsilon for gradient norm
 * @param epsilon
 */
void IGradientMethod::setEpsilon1(double epsilon)
{
    m_epsilon1 = epsilon;
}

/**
 * @brief Epsilon for distance between points
 * @return
 */
double IGradientMethod::epsilon2() const
{
    return m_epsilon2;
}

/**
 * @brief Epsilon for distance between points
 * @param epsilon
 */
void IGradientMethod::setEpsilon2(double epsilon)
{
    m_epsilon2 = epsilon;
}

/**
 * @brief Epsilon for distance between points
 * @param epsilon
 */
void IGradientMethod::setEpsilon3(double epsilon)
{
    m_epsilon3 = epsilon;
}

/**
 * @brief Epsilon for distance between points
 * @return
 */
double IGradientMethod::epsilon3() const
{
    return m_epsilon3;
}

void IGradientMethod::setR1MinimizeEpsilon(double step, double epsilon)
{
    min_step = step;
    min_epsilon = epsilon;
}

int IGradientMethod::count() const
{
    return iterationCount;
}

void IGradientMethod::setNormalize(bool normalize)
{
    this->m_normalize = normalize;
}

void IGradientMethod::showEndMessage(bool showEndMessage)
{
    mshowEndMessage = showEndMessage;
}
