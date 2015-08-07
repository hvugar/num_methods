#include "gradient_prj.h"

ProjectionGradient::ProjectionGradient() : SteepestDescentGradient(), a(-INFINITY), b(INFINITY)
{
}

ProjectionGradient::~ProjectionGradient()
{

}

void ProjectionGradient::calculate()
{
    iterationCount = 0;

    double distance = 0.0;
    do
    {
        /* calculating function gradient at current point */
        m_fn->gradient(grad_step, m_x, m_g);

        /* if norm of gradinet at current point is less than epsilon break. no minimize */
        double gradient_norm = m_g.L2Norm();
        if (gradient_norm < epsilon()) break;

        iterationCount++;

        /* Normalize vector */
        m_g.L2Normalize();

        /* R1 minimization in direct of antigradient */
        m_alpha = minimize();

        if (printer != NULL) printer->print(iterationCount, m_x, m_g, m_alpha, function());

        distance = 0.0;
        for (unsigned int i=0; i<m_x.size(); i++)
        {
            double x = m_x[i];
            m_x[i] = m_x[i] - m_alpha * m_g[i];

            if (m_x[i] < a) { m_x[i] = a; }
            if (m_x[i] > b) { m_x[i] = b; }

            distance += (m_x[i]-x)*(m_x[i]-x);
        }
        distance = sqrt(distance);

        if (distance < epsilon()) break;

        /* calculating distance previous and new point */
    } while (true);
}

void ProjectionGradient::calculate(DoubleVector& x)
{

}

double ProjectionGradient::fx(double alpha)
{
    DoubleVector x(m_x.size(), 0.0);
    for (unsigned int i=0; i<m_x.size(); i++)
    {
        x[i] = m_x[i] - alpha * m_g[i];
        if ( x[i] < a ) x[i] = a;
        if ( x[i] > b ) x[i] = b;
    }
    return m_fn->fx(x);
}

