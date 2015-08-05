#include "gradient_sd.h"

SteepestDescentGradient::SteepestDescentGradient() : GradientMethod()
{}

SteepestDescentGradient::~SteepestDescentGradient()
{}

double SteepestDescentGradient::minimize()
{
    double alpha0 = 0.0;
    R1Minimize r1;
    r1.setFunction(this);
    r1.setX0(alpha0);
    r1.setStep(min_step);
    r1.setEpsilon(min_epsilon);
    r1.straightLineSearch();
    double alpha = r1.goldenSectionSearch();
    if (fx(alpha) > fx(alpha0)) alpha = alpha0;
    return alpha;
}

void SteepestDescentGradient::calculate()
{
    iterationCount = 0;
    double distance = 0.0;

    do
    {
        /* calculating function gradient at current point */
        //calculateGradient();
        m_fn->gradient(grad_step, m_x, m_g);

        /* if gradinet norm at current point is less than epsilon then break. no minimize */
        double gradient_norm = m_g.L2Norm();
        if (gradient_norm < epsilon()) break;

        iterationCount++;

        /* Normalize vector */
        m_g.L2Normalize();

        /* R1 minimization in direct of antigradient */
        m_alpha = minimize();

        if (printer != NULL)
        {
            printer->print(iterationCount, m_x, m_g, m_alpha, function());
        }

        distance = 0.0;
        for (unsigned int i=0; i<m_x.size(); i++)
        {
            double x = m_x[i];
            m_x[i] = m_x[i] - m_alpha * m_g[i];

            // calculating distance
            distance += (m_x[i]-x)*(m_x[i]-x);
        }
        distance = sqrt(distance);

        /* calculating distance previous and new point */
    } while (distance > epsilon());
}

void SteepestDescentGradient::calculate(DoubleVector &x)
{
}

double SteepestDescentGradient::fx(double alpha)
{
    DoubleVector x(m_x.size());
    for (unsigned int i=0; i < m_x.size(); i++)
    {
        x[i] = m_x[i] - alpha * m_g[i];
    }
    return m_fn->fx(x);
}
