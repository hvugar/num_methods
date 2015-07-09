#include "gradient.h"

Gradient::Gradient() : mfn(NULL)
{
    malpha = 0.0;
    mepsilon = 0.0;
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
    mfn = f;
}

RnFunction* Gradient::function() const
{
    return mfn;
}

void Gradient::setX(const std::vector<double> &x)
{
    mx = x;
    mg = x;
}

const std::vector<double>& Gradient::x() const
{
    return mx;
}

double Gradient::epsilon() const
{
    return mepsilon;
}

void Gradient::setEpsilon(double epsilon)
{
    mepsilon = epsilon;
}

void Gradient::calculateGradient()
{
    double h = grad_step;
    for (unsigned i=0; i<mx.size(); i++)
    {
        mx[i] = mx[i] - h;
        double f1 = mfn->fx(mx);
        mx[i] = mx[i] + 2*h;
        double f2 = mfn->fx(mx);
        mx[i] = mx[i] - h;
        mg[i] = (f2 - f1) / (2 * h);
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
    for (unsigned int i=0; i<mg.size(); i++)
    {
        grad_norm = grad_norm + mg[i]*mg[i];
    }
    grad_norm = sqrt(grad_norm);
    return grad_norm;
}

double Gradient::distance() const
{
    double dist = 0.0;
    for (unsigned int i=0; i<mx.size(); i++)
    {
        dist = dist + (malpha * mg[i]) * (malpha * mg[i]);
    }
    dist = sqrt(dist);
    return dist;
}
