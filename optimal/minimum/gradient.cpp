#include "gradient.h"

Gradient::Gradient()
{
    mf = 0;
    mepsilon = 0.000001;
    mcount = 0;
    malpha = 1.0;
}

Gradient::~Gradient()
{
}

void Gradient::setF(RnFunction *f)
{
    mf = f;
}

RnFunction* Gradient::f() const
{
    return mf;
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

void Gradient::calcGradient()
{
    double h = grad_step;
    for (unsigned i=0; i<mx.size(); i++)
    {
        mx[i] = mx[i] - h;
        double f1 = mf->fx(mx);
        mx[i] = mx[i] + 2*h;
        double f2 = mf->fx(mx);
        mx[i] = mx[i] - h;
        mg[i] = (f2 - f1) / (2 * h);
    }
}

class Argmin : public R1Minimize
{
public:
    std::vector<double>& x;
    std::vector<double>& g;
    Gradient* gm;
    Argmin(std::vector<double>& x, std::vector<double>& g, Gradient* gm) :
        R1Minimize(), x(x), g(g), gm(gm) {}

protected:
    double fx(double alpha) {
        std::vector<double> x1;
        for (unsigned int i=0; i < x.size(); i++) x1.push_back(x[i] - alpha * g[i]);
        return gm->f()->fx(x1);
    }
};

double Gradient::minimize()
{
    double x0 = 0.0;
    Argmin r1(mx, mg, this);
    r1.setX0(x0);
    r1.setStep(min_step);
    r1.setEpsilon(min_epsilon);
    r1.straightLineSearch();
    double alpha = r1.goldenSectionSearch();
    return alpha;
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
    return mcount;
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
