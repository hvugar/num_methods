#include "gradient_cjt.h"
#include "printer.h"
#include "projection.h"
#include "function.h"
#include <math.h>

ConjugateGradient::ConjugateGradient() : GradientMethod()
{
    setNormalize(true);
    r1m.setFunction(this);
    setAlgorithm(FLETCHER_REEVES);
    setResetIteration(true);
}

ConjugateGradient::~ConjugateGradient()
{}

void ConjugateGradient::calculate(DoubleVector& x)
{
    unsigned int k = 0;
    double grad_mod_cur = 0.0;
    double grad_mod_prv = 0.0;
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
    m_gr->gradient(x, g);

    /**************************************************************************************
     * Checking for gradient vector norm that is less of epsilon.
     * If gradient vector norm in current point is less than epsilon then break the iteration.
     * Finish minimization.
     **************************************************************************************/
    double gradient_norm = g.L2Norm();
    if (gradient_norm < epsilon1())
    {
        if (m_printer != NULL) m_printer->print(iterationCount, x, g, m_fn->fx(x), GradientMethod::BREAK_FIRST_ITERATION);
        if (mshowEndMessage) puts("Optimisation ends, because norm of gradient is less than epsilon...");
        return;
    }
    f1 = m_fn->fx(x);
    if (m_printer != NULL) m_printer->print(iterationCount, x, g, f1, GradientMethod::FIRST_ITERATION);

    do
    {
        iterationCount++;

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
         * Normalization of a vector
         **************************************************************************************/
        if (m_normalize) s.L2Normalize();

        /**************************************************************************************
         * One-dimensional minimization along the direction of a anti-gradient
         **************************************************************************************/
        alpha = minimize(x, s);
        printf("After minimization: %.8f\n", alpha);

        /**************************************************************************************
         * Calculation current point.
         **************************************************************************************/
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

        /**************************************************************************************
         * Gradient of objectiv function in current point
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
            if (m_printer != NULL) m_printer->print(iterationCount, x, g, f2, GradientMethod::BREAK_GRADIENT_NORM_LESS);
            if (mshowEndMessage) puts("Optimisation ends, because norm of gradient is less than epsilon...");
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
            if (m_printer != NULL) m_printer->print(iterationCount, x, g, f2, GradientMethod::BREAK_DISTANCE_LESS);
            if (mshowEndMessage) puts("Optimisation ends, because distance between last and current point less than epsilon...");
            break;
        }

        /**************************************************************************************
         * Printing iteration information.
         **************************************************************************************/
        if (m_printer != NULL) m_printer->print(iterationCount, x, g, f2, GradientMethod::NEXT_ITERATION);

        f1 = f2;

    } while (true);

    g.clear();
    s.clear();
}

double ConjugateGradient::minimize(const DoubleVector &x, const DoubleVector &s) const
{
    C_UNUSED(x);
    C_UNUSED(s);

//    for (int i=-100; i<=+100; i++)
//    {
//        double a_lpha = i*0.01;
//        printf("%f\n", fx(a_lpha));
//    }

    double alpha0 = +min_step;
    double a,b,alpha;

    //stranghLineSearch(alpha0, min_step, a, b, this);
    //goldenSectionSearch(a, b, alpha, this, min_epsilon);
    //if (fx(alpha) > fx(alpha0)) alpha = alpha0;

    double fxa, fxb;
    bool unimodal;
    r1m.straightLineSearch(alpha0, min_step, a, b, fxa, fxb, unimodal);
    //r1m.swann(alpha0, min_step, a, b, fxa, fxb, unimodal);
    puts("---");
    if (unimodal)
    {
        //printf("%s\n", "unimodal");
        r1m.goldenSectionSearch(alpha, a, b, min_epsilon);
        //printf("alpha %.8f\n", alpha);
    }
    else
    {
        //double h = (b-a)/100.0;
        //for (unsigned int i=0; i<=100; i++)
        //{
        //    double x = a + i*h;
        //    double y = this->fx(x);
        //    printf("%f %f\n", x, y);
        //}
        //double x = -0.3;
        //double y = this->fx(x);
        //printf("%f %f\n", x, y);
        //exit(-1);

        fxa < fxb ? alpha = a : alpha = b;
    }
    //printf("alpha %.8f fx %.8f alpha0 %.8f fx %.8f\n", alpha, fx(alpha), alpha0, fx(alpha0));
    if (fx(alpha) > fx(alpha0)) alpha = alpha0;

    return alpha;
}

double ConjugateGradient::fx(double alpha) const
{
    const DoubleVector &x = *mx;
    const DoubleVector &s = *ms;
    unsigned int n = x.length();

    //printf("--- o: %f %f %f %f c: %f %f %f %f\n", x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
    DoubleVector cx = x;
    for (unsigned int i=0; i<n; i++)
    {
        cx[i] = x[i] + alpha * s[i];
        if (m_projection != NULL) m_projection->project(cx, i);
    }

    //printf("--- o: %f %f %f %f c: %f %f %f %f\n", cx[8], cx[9], cx[10], cx[11], cx[12], cx[13], cx[14], cx[15]);
    double f = m_fn->fx(cx);
    //printf("alpha: %12.8f fx: %18.8f\n", alpha, f);

    return f;
}

void ConjugateGradient::setAlgorithm(Algorithm algorithm)
{
    malgoritm = algorithm;
}

void ConjugateGradient::setResetIteration(bool reset)
{
    mResetIteration = reset;
}

R1FxMinimizer &ConjugateGradient::R1Minimizer()
{
    return r1m;
}

const R1FxMinimizer& ConjugateGradient::R1Minimizer() const
{
    return r1m;
}

