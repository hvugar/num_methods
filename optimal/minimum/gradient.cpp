#include "gradient.h"

GradientMethod::GradientMethod() : m_fn(NULL)
{
    m_alpha = 0.0;
    m_epsilon = 0.0;
    grad_step = 0.0;
    min_step = 0.0;
    min_epsilon = 0.0;
    iterationCount = M = 0;
}

GradientMethod::~GradientMethod()
{
}

void GradientMethod::setFunction(RnFunction *f)
{
    m_fn = f;
}

RnFunction* GradientMethod::function() const
{
    return m_fn;
}

void GradientMethod::setX(const DoubleVector &x)
{
    m_x = x;
    m_g.resize(x.size(), 0.0);
}

const DoubleVector& GradientMethod::x() const
{
    return m_x;
}

double GradientMethod::epsilon() const
{
    return m_epsilon;
}

void GradientMethod::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void GradientMethod::calculateGradient()
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

void GradientMethod::setR1MinimizeEpsilon(double step, double epsilon)
{
    min_step = step;
    min_epsilon = epsilon;
}

void GradientMethod::setGradientStep(double step)
{
    grad_step = step;
}

int GradientMethod::count() const
{
    return iterationCount;
}

double GradientMethod::gradientNorm() const
{
    double grad_norm = 0.0;
    for (unsigned int i=0; i<m_g.size(); i++)
    {
        grad_norm = grad_norm + m_g[i]*m_g[i];
    }
    grad_norm = sqrt(grad_norm);
    return grad_norm;
}

double GradientMethod::distance() const
{
    double dist = 0.0;
    for (unsigned int i=0; i<m_x.size(); i++)
    {
        dist = dist + (m_alpha * m_g[i]) * (m_alpha * m_g[i]);
    }
    dist = sqrt(dist);
    return dist;
}
