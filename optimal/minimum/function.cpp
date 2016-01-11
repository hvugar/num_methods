#include "function.h"

void RnFunction::Gradient(RnFunction *f, double step, const DoubleVector &x, DoubleVector &g)
{
    double h = step;
    for (unsigned int i=0; i<g.size(); i++)
    {
        double cx = x[i];
        const_cast<DoubleVector&>(x)[i] = cx - h;
        double f1 = f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx + h;
        double f2 = f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx;

        g[i] = (f2 - f1) / (2 * h);
    }
}

void IGradient::Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g)
{
    unsigned int n = g.size();
    for (unsigned int i=0; i<n; i++)
    {
        double cx = x[i];
        const_cast<DoubleVector&>(x)[i] = cx - step;
        double f1 = 0.0;//f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx + step;
        double f2 = 0.0;//f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx;
        g[i] = (f2 - f1) / (2 * step);
    }
}

void RnGradient::gradient(const DoubleVector &x, DoubleVector &g)
{
    unsigned int n = g.size();
    for (unsigned int i=0; i<n; i++)
    {
        double cx = x[i];
        const_cast<DoubleVector&>(x)[i] = cx - step;
        double f1 = 0.0;//f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx + step;
        double f2 = 0.0;//f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx;
        g[i] = (f2 - f1) / (2 * step);
    }
}
