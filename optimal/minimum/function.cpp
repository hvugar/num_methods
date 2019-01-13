#include "function.h"

R1Function::~R1Function() {}

R2Function::~R2Function() {}

R3Function::~R3Function() {}

RnFunction::~RnFunction() {}

VectorFunction::~VectorFunction() {}

VectorRnFunction::~VectorRnFunction() {}

MatrixR1Function::~MatrixR1Function() {}

MatrixRnFunction::~MatrixRnFunction() {}

IGradient::~IGradient() {}

double sgn(double x)
{
    if (x < 0.0) return -1.0;
    else if (x > 0.0) return +1.0;
    return 0.0;
}

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

void IGradient::Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g, unsigned int *inx, unsigned int size)
{
    double h = step;
    for (unsigned int j=0; j<size; j++)
    {
        unsigned int i = inx[j];
        double cx = x[i];
        const_cast<DoubleVector&>(x)[i] = cx - h;
        double f1 = f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx + h;
        double f2 = f->fx(x);
        const_cast<DoubleVector&>(x)[i] = cx;
        g[i] = (f2 - f1) / (2.0 * h);
    }
}

void IGradient::Gradient(const RnFunction *f, double step, const DoubleVector &x, DoubleVector &g, unsigned int start, unsigned int end)
{
    DoubleVector &var_x = const_cast<DoubleVector&>(x);
    for (unsigned int i=start; i<=end; i++)
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

