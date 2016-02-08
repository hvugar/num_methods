#include "gradient_sd.h"

SteepestDescentGradient::SteepestDescentGradient() : GradientMethod()
{
    setNormalize(true);
}

SteepestDescentGradient::~SteepestDescentGradient()
{}

void SteepestDescentGradient::calculate(DoubleVector &x)
{
    unsigned int n = x.size();

    DoubleVector g(n);

    mx = &x;
    mg = &g;
    iterationCount = 0;

    double distance = 0.0;
    double alpha = 0.0;
//    double f1 = 0.0;
//    double f2 = 0.0;

    do
    {
        /* calculating function gradient at current point */
        m_gr->gradient(x, g);

        /* if gradinet norm at current point is less than epsilon then break. no minimize */
        double gradient_norm = g.L2Norm();
        if (gradient_norm < epsilon1())
        {
            puts("Optimisation ends, because L2 norm of gradient is less than epsilon...");
            break;
        }


        iterationCount++;

        /* Normalize vector */
        if (m_normalize) g.L2Normalize();

        /* R1 minimization in direct of antigradient */
        alpha = minimize(x, g);

        if (m_printer != NULL)
        {
            m_printer->print(iterationCount, x, g, alpha, function());
        }

        distance = 0.0;
        double f1 = m_fn->fx(x);
        for (unsigned int i=0; i < n; i++)
        {
            double cx = x[i];
            x[i] = x[i] - alpha * g[i];

            // calculating distance
            distance += (x[i]-cx)*(x[i]-cx);
        }
        distance = sqrt(distance);
        double f2 = m_fn->fx(x);

        /* calculating distance previous and new point */
        if (distance < epsilon2() && fabs(f2 - f1) < epsilon3())
        {
            puts("Optimisation ends, because distance beetween last and current point less than epsilon...");
            break;
        }

    } while (true);

    g.clear();
}

double SteepestDescentGradient::fx(double alpha)
{
    DoubleVector &x = *mx;
    DoubleVector &g = *mg;
    unsigned int n = x.size();

    DoubleVector cx(n);
    for (unsigned int i=0; i<n; i++)
    {
        cx[i] = x[i] - alpha * g[i];
        if (m_projection != NULL) m_projection->project(cx, i);
    }

    return m_fn->fx(cx);
}

double SteepestDescentGradient::minimize(const DoubleVector &x, const DoubleVector &g)
{
    C_UNUSED(x);
    C_UNUSED(g);

    double alpha0 = 0.0;
    double a,b,alpha;

    stranghLineSearch(alpha0, min_step, a, b, this);
    goldenSectionSearch(a, b, alpha, this, min_epsilon);

    if (this->fx(alpha) > this->fx(alpha0)) alpha = alpha0;
    return alpha;
}
