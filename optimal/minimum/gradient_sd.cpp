#include "gradient_sd.h"

SteepestDescentGradient::SteepestDescentGradient() : Gradient()
{
}

SteepestDescentGradient::~SteepestDescentGradient()
{}

double SteepestDescentGradient::minimize()
{
    double alpha0 = 0.0;
    R1Minimize r1;
    r1.setF(this);
    r1.setX0(alpha0);
    r1.setStep(min_step);
    r1.setEpsilon(min_epsilon);
    r1.straightLineSearch();
    double alpha = r1.goldenSectionSearch();
    if ( this->fx(alpha) > this->fx(alpha0) ) alpha = alpha0;
    return alpha;
}

void SteepestDescentGradient::calculate()
{
    k = 0;
    do
    {
        /* calculating function gradient at current point */
        calculateGradient();

        /* if gradinet norm at current point is less than epsilon then break. no minimize */
        double gradient_norm = gradientNorm();
        if (gradient_norm < epsilon())
            break;

        k++;

        /* calculating unit vectors */
        for (unsigned int i=0; i<m_g.size(); i++) m_g[i] = m_g[i] / gradient_norm;

        /* R1 minimization in direct of antigradient */
        m_alpha = minimize();

        print();

        for (unsigned int i=0; i<m_x.size(); i++)
        {
            m_x[i] = m_x[i] - m_alpha * m_g[i];
        }

        /* calculating distance previous and new point */
    } while (distance() > epsilon());
}

double SteepestDescentGradient::fx(double alpha)
{
    std::vector<double> x1 = m_x;
    for (unsigned int i=0; i < m_x.size(); i++)
    {
        x1[i] = m_x[i] - alpha * m_g[i];
    }
    return m_fn->fx(x1);
}

void SteepestDescentGradient::print()
{
    if (k == 1)
    {
        printf("No\t|x1      \t|x2      \t|f(x)      \t|grad1      \t|grad2      \t|grad_norm  \t|alpha  \t");
        printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+-------------\n");
    }

    double y = function()->fx(m_x);
    double nr = gradientNorm();

    printf("%d\t", k);
    m_x[0]>=0.0 ? printf("|+%.10f\t", fabs(m_x[0])) : printf("|%.10f\t", m_x[0]);
    m_x[1]>=0.0 ? printf("|+%.10f\t", fabs(m_x[1])) : printf("|%.10f\t", m_x[1]);
    y>=0.0 ? printf("|%+10.6f\t", y) : printf("|%10.6f\t", y);
    m_g[0]>=0.0 ? printf("|%+10.6f\t", m_g[0]) : printf("|%10.6f\t", m_g[0]);
    m_g[1]>=0.0 ? printf("|%+10.6f\t", m_g[1]) : printf("|%10.6f\t", m_g[1]);
    nr>=0.0 ? printf("|%+10.6f\t", nr) : printf("|%10.6f\t", nr);
    m_alpha>=0.0 ? printf("|%+10.6f\t", m_alpha) : printf("|%10.6f\t", m_alpha);
    printf("\n");
}
