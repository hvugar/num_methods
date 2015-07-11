#include "gradient.h"

Gradient::Gradient() : m_fn(NULL)
{
    m_alpha = 0.0;
    m_epsilon = 0.0;
    grad_step = 0.0;
    min_step = 0.0;
    min_epsilon = 0.0;
    k = M = 0;
}

Gradient::~Gradient()
{
}

void Gradient::setFunction(RnFunction *f)
{
    m_fn = f;
}

RnFunction* Gradient::function() const
{
    return m_fn;
}

void Gradient::setX(const std::vector<double> &x)
{
    m_x = x;
    m_g.resize(x.size(), 0.0);
}

const std::vector<double>& Gradient::x() const
{
    return m_x;
}

double Gradient::epsilon() const
{
    return m_epsilon;
}

void Gradient::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void Gradient::calculateGradient()
{
    double h = grad_step;
    for (unsigned i=0; i<m_x.size(); i++)
    {
        m_x[i] = m_x[i] - h;
        double f1 = m_fn->fx(m_x);
        m_x[i] = m_x[i] + 2*h;
        double f2 = m_fn->fx(m_x);
        m_x[i] = m_x[i] - h;

        m_g[i] = (f2 - f1) / (2 * h);
    }
}

void Gradient::setR1MinimizeEpsilon(double step, double epsilon)
{
    min_step = step;
    min_epsilon = epsilon;
}

void Gradient::setGradientStep(double step)
{
    grad_step = step;
}

int Gradient::count() const
{
    return k;
}

double Gradient::gradientNorm() const
{
    double grad_norm = 0.0;
    for (unsigned int i=0; i<m_g.size(); i++)
    {
        grad_norm = grad_norm + m_g[i]*m_g[i];
    }
    grad_norm = sqrt(grad_norm);
    return grad_norm;
}

double Gradient::distance() const
{
    double dist = 0.0;
    for (unsigned int i=0; i<m_x.size(); i++)
    {
        dist = dist + (m_alpha * m_g[i]) * (m_alpha * m_g[i]);
    }
    dist = sqrt(dist);
    return dist;
}
