#include "function.h"

double sgn(double x)
{
    if (x < 0.0) return -1.0;
    else if (x > 0.0) return +1.0;
    return 0.0;
}


R1Function::~R1Function() {}

R2Function::~R2Function() {}

R3Function::~R3Function() {}

RnFunction::~RnFunction() {}

VectorFunction::~VectorFunction() {}

VectorRnFunction::~VectorRnFunction() {}

MatrixR1Function::~MatrixR1Function() {}

MatrixRnFunction::~MatrixRnFunction() {}

IGradient::~IGradient() {}

void IGradient::Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g)
{
    double h = step;
    for (unsigned int i=0; i<g.length(); i++)
    {
        double cx = x[i];
        const_cast<DoubleVector&>(x)[i] = cx - h;
        double f1 = f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx + h;
        double f2 = f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx;
        g[i] = (f2 - f1) / (2.0 * h);
    }
}

void IGradient::Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g, const std::vector<size_t> &index)
{
    const size_t size = index.size();
    double h = step;
    for (unsigned int j=0; j<size; j++)
    {
        size_t i = index[j];
        double cx = x[i];
        const_cast<DoubleVector&>(x)[i] = cx - h;
        double f1 = f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx + h;
        double f2 = f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx;
        g[i] = (f2 - f1) / (2.0 * h);
    }
}

void IGradient::Gradient(double step, const RnFunction *f, const DoubleVector &x, DoubleVector &g, const std::vector<size_t> &index)
{
    const size_t size = index.size();
    for (unsigned int j=0; j<size; j++)
    {
        size_t i = index[j];
        double cx = x[i];
        const_cast<DoubleVector&>(x)[i] = cx - step;
        double f1 = f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx + step;
        double f2 = f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx;
        g[i] = (f2 - f1) / (2.0 * step);
    }
}

void IGradient::Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g, size_t start, size_t end)
{
    DoubleVector &var_x = const_cast<DoubleVector&>(x);
    for (size_t i=start; i<=end; i++)
    {
        double cx = var_x[i];
        var_x[i] = cx - step;
        double f1 = f->fx(var_x);
        var_x[i] = cx + step;
        double f2 = f->fx(var_x);
        g[i] = (f2 - f1) / (2.0 * step);
        var_x[i] = cx;
    }
}

void IGradient::GradientL(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g, size_t start, size_t end)
{
    DoubleVector &var_x = const_cast<DoubleVector&>(x);
    for (size_t i=start; i<=end; i++)
    {
        double cx = var_x[i];
        var_x[i] = cx - step;
        double f1 = f->fx(var_x);
        var_x[i] = cx;
        double f2 = f->fx(var_x);
        g[i] = (f2 - f1) / step;
    }
}

void IGradient::GradientR(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g, size_t start, size_t end)
{
    DoubleVector &var_x = const_cast<DoubleVector&>(x);
    for (size_t i=start; i<=end; i++)
    {
        double cx = var_x[i];
        var_x[i] = cx + step;
        double f1 = f->fx(var_x);
        var_x[i] = cx;
        double f2 = f->fx(var_x);
        g[i] = (f1 - f2) / step;
    }
}

void IGradient::GradientF(const RnFunction *f, double factor, const DoubleVector &x, DoubleVector &g, size_t start, size_t end)
{
    DoubleVector &var_x = const_cast<DoubleVector&>(x);
    for (size_t i=start; i<=end; i++)
    {
        double cx = var_x[i];
        var_x[i] = cx * (1.0-factor);
        double f1 = f->fx(var_x);
        var_x[i] = cx + (1.0+factor);
        double f2 = f->fx(var_x);
        g[i] = (f2 - f1) / (2.0 * factor*cx);
        var_x[i] = cx;
    }
}
