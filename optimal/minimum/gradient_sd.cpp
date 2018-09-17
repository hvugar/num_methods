#include "gradient_sd.h"
#include "printer.h"
#include "projection.h"
#include "function.h"
#include <math.h>

SteepestDescentGradient::SteepestDescentGradient() : GradientMethod()
{
    setNormalize(true);
    //    r1m.setFunction(this);
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

    // Gradient of objectiv function in current point
    m_gr->gradient(x, g);
    f1 = m_fn->fx(x);

    /* checking gradient norm */
    /* if gradient norm at current point is less than epsilon then return. no minimize */
    double gradient_norm = g.L2Norm();
    if (gradient_norm < epsilon1())
    {
        if (m_printer != NULL) m_printer->print(m_iteration_count, x, g, f1, alpha, BREAK_FIRST_ITERATION);
        if (m_show_end_message) puts("Optimisation ends, because norm of gradient is less than epsilon...");
        return;
    }
    if (m_printer != NULL) m_printer->print(m_iteration_count, x, g, f1, alpha, FIRST_ITERATION);

    do
    {
        m_iteration_count++;

        /* Normalize vector */
        if (m_normalize) g.L2Normalize();

        /* R1 minimization in direct of antigradient */
        alpha = minimize(x, g);

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
            if (m_printer != NULL) m_printer->print(m_iteration_count, x, g, f2, alpha, GradientMethod::BREAK_GRADIENT_NORM_LESS);
            if (m_show_end_message) puts("Optimisation ends, because norm of gradient is less than epsilon...");
            break;
        }

        /* calculating distance previous and new point */
        if (distance < epsilon2() && fabs(f2 - f1) < epsilon3())
        {
            if (m_printer != NULL) m_printer->print(m_iteration_count, x, g, f2, alpha, GradientMethod::BREAK_DISTANCE_LESS);
            if (m_show_end_message) puts("Optimisation ends, because distance between last and current point less than epsilon...");
            break;
        }

        if (m_printer != NULL) m_printer->print(m_iteration_count, x, g, f2, alpha, GradientMethod::NEXT_ITERATION);

        f1 = f2;

    } while (true);

    g.clear();
}

double SteepestDescentGradient::minimize(const DoubleVector &x, const DoubleVector &g) const
{
    C_UNUSED(x);
    C_UNUSED(g);

//    for (int i = -100; i <= +100;  i++)
//    {
//        double a = 0.001*i;
//        double f = fx(a);
//        printf("%f %f\n", a, f);
//        //printf("   %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", g[0], g[1], g[2], g[3],
//        //        g[4], g[5], g[6], g[7], g[8], g[9], g[10], g[11], g[12], g[13], g[14], g[15]);
//        //printf("   %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", x[0], x[1], x[2], x[3],
//        //        x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
//        //printf("   %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", x[0]-a*g[0], x[1]-a*g[1], x[2]-a*g[2], x[3]-a*g[3],
//        //                                                            x[4]-a*g[4], x[5]-a*g[5], x[6]-a*g[6], x[7]-a*g[7],
//        //                                                            x[8]-a*g[8], x[9]-a*g[9], x[10]-a*g[10], x[11]-a*g[11],
//        //                                                            x[12]-a*g[12], x[13]-a*g[13], x[14]-a*g[14], x[15]-a*g[15]);
//    }
//    exit(-1);

    double alpha0 = min_step;
    double a,b,alpha;

    double fxa, fxb;
    bool unimodal;

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

    //printf("mg %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", g[0], g[1], g[2], g[3],
    //        g[4], g[5], g[6], g[7], g[8], g[9], g[10], g[11], g[12], g[13], g[14], g[15]);
    //printf("xx %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", x[0], x[1], x[2], x[3],
    //        x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15]);
    //printf("cx %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", cx[0], cx[1], cx[2], cx[3],
    //        cx[4], cx[5], cx[6], cx[7], cx[8], cx[9], cx[10], cx[11], cx[12], cx[13], cx[14], cx[15]);
    return m_fn->fx(cx);
}
