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
    do
    {
        /* calculating function gradient at current point */
        calculateGradient();

        /* if norm of gradinet at current point is less than epsilon break. no minimize */
        double gn = gradientNorm();
        if (gn < epsilon())
            break;

        iterationCount++;

        for (unsigned int i=0; i<m_g.size(); i++) m_g[i] = m_g[i] / gn;

        /* R1 minimization in direct of antigradient */
        m_alpha = minimize();

        print();

        for (unsigned int i=0; i<m_x.size(); i++)
        {
            m_x[i] = m_x[i] - m_alpha * m_g[i];

            if (m_x[i] < a) { m_x[i] = a; }
            if (m_x[i] > b) { m_x[i] = b; }
        }

        /* calculating distance previous and new point */
    } while (distance() > epsilon());
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

