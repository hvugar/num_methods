#include "gradient_cs.h"
#include "printer.h"
#include "projection.h"
#include "function.h"
#include <math.h>

ConstStepGradient::ConstStepGradient() : GradientMethod()
{
    setNormalize(true);
}

ConstStepGradient::~ConstStepGradient()
{}

void ConstStepGradient::calculate(DoubleVector &x)
{
    unsigned int n = x.length();

    m_iterationCount = 0;
    double distance = 0.0;
    double alpha = 0.0;

    DoubleVector g(n);

    do
    {
        /* calculating function gradient at current point */
        m_gr->gradient(x, g);

        /* if gradinet norm at current point is less than optimality tolerance then break. no minimize */
        double gradNorm = g.L2Norm();
        if (gradNorm < optimalityTolerance())
        {
            puts("Optimisation ends, because L2 norm of gradient is less than epsilon...");
            break;
        }

        m_iterationCount++;

        /* Normalize vector */
        if (m_normalize) g.L2Normalize();

        /* R1 minimization in direct of antigradient */
        alpha = minimize(x, g);

        //if (m_printer != NULL) m_printer->print(iterationCount, x, g, 0.0);

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
        if (distance < stepTolerance() && fabs(f2 - f1) < functionTolerance())
        {
            puts("Optimisation ends, because distance between last and current point less than epsilon...");
            break;
        }

    } while (true);
}

double ConstStepGradient::minimize(const DoubleVector &x, const DoubleVector &g) const
{
    unsigned int n = x.length();

    DoubleVector cx(n);

    double alpha = 0.01;
    for (unsigned int i=0; i < n; i++)
    {
        cx[i] = x[i] - alpha * g[i];
    }

    while (m_fn->fx(cx) > m_fn->fx(x))
    {
        alpha = alpha / 2.0;
        for (unsigned int i=0; i < n; i++)
        {
            cx[i] = x[i] - alpha * g[i];
        }
    }
    return alpha;
}
