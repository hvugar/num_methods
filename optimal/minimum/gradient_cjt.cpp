#include "gradient_cjt.h"

ConjugateGradient::ConjugateGradient() : GradientMethod()
{
}

ConjugateGradient::~ConjugateGradient()
{}

void ConjugateGradient::setX(const DoubleVector &x)
{
    GradientMethod::setX(x);
    s.resize(x.size(), 0.0);
}

void ConjugateGradient::calculate()
{
    unsigned int n = 0;
    double gr0_mod = 0.0;
    double gr1_mod = 0.0;
    double gr2_mod = 0.0;

    do
    {
        // Gradient of objectiv function in current point
        calculateGradient();

        double gradNorm = m_g.L2Norm();
        if (gradNorm < epsilon())
        {
            puts("Optimisation ends, because L2 norm of gradient is less than epsilon...");
            break;
        }

        iterationCount++;

        // Module of gradient
        gr0_mod = 0.0;
        for (unsigned int i=0; i<m_g.size(); i++) gr0_mod = gr0_mod + (m_g[i]*m_g[i]);
        //gr0_mod = sqrt(gr0_mod);

        if (n == 0)
        {
            gr1_mod = gr0_mod;
            // First direction is antigradient
            for (unsigned int i=0; i<m_g.size(); i++) s[i] = -m_g[i];
        }
        else
        {
            gr2_mod = gr0_mod;
            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;
            // Direction in next iteration (n != 0)
            for (unsigned int i=0; i<m_g.size(); i++) s[i] = -m_g[i] + s[i] * w;
        }

        //double sn = 0.0;
        //for (unsigned int i=0; i<s.size(); i++) sn = sn + (s[i]*s[i]);
        //sn = sqrt(sn);
        //for (unsigned int i=0; i<s.size(); i++) s[i] = s[i] / sn;

        //s.L2Normalize();

        m_alpha = minimize();

        print();

        for (unsigned int i=0; i<m_x.size(); i++)
        {
            m_x[i] = m_x[i] + m_alpha * s[i];
        }

        if ( n == (m_x.size()-1) ) { n = 0; } else { n++; }

        /* calculating distance previous and new point */
    } while (distance() > epsilon());
}

double ConjugateGradient::minimize()
{
    double alpha0 = 0.0;
    R1Minimize r1;
    r1.setF(this);
    r1.setX0(alpha0);
    r1.setStep(min_step);
    r1.setEpsilon(min_epsilon);
    r1.straightLineSearch();
    double alpha = r1.goldenSectionSearch();
    if ( fx(alpha) > fx(alpha0) ) alpha = alpha0;
    return alpha;
}

double ConjugateGradient::distance() const
{
    double dist = 0.0;
    for (unsigned int i=0; i<m_x.size(); i++)
    {
        dist = dist + ((m_alpha * s[i]) * (m_alpha * s[i]));
    }
    dist = sqrt(dist);
    return dist;
}

double ConjugateGradient::fx(double alpha)
{
    DoubleVector x(m_x.size());
    for (unsigned int i=0; i < m_x.size(); i++)
    {
        x[i] = m_x[i] + alpha * s[i];
    }
    return m_fn->fx(x);
}

void ConjugateGradient::print()
{
    if (iterationCount == 1)
    {
        printf("No\t|x1      \t|x2      \t|f(x)      \t|s1      \t|s2      \t|grad_norm  \t|alpha  \t");
        printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+-------------\n");
    }

    double y = function()->fx(m_x);
    double nr = gradientNorm();

    printf("%d\t", iterationCount);
    m_x[0]>=0.0 ? printf("|+%.10f\t", fabs(m_x[0])) : printf("|%.10f\t", m_x[0]);
    m_x[1]>=0.0 ? printf("|+%.10f\t", fabs(m_x[1])) : printf("|%.10f\t", m_x[1]);
    y>=0.0 ? printf("|%+.10f\t", y) : printf("|%.10f\t", y);
    s[0]>=0.0 ? printf("|%+.10f\t", s[0]) : printf("|%0.10f\t", s[0]);
    s[1]>=0.0 ? printf("|%+.10f\t", s[1]) : printf("|%0.10f\t", s[1]);
    nr>=0.0 ? printf("|%+.10f\t", nr) : printf("|%.10f\t", nr);
    m_alpha>=0.0 ? printf("|%+.10f\t", m_alpha) : printf("|%.10f\t", m_alpha);
    printf("\n");
}
