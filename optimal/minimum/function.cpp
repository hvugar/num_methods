#include "function.h"

void IGradient::Gradient(RnFunction *f, double step, const DoubleVector &x, DoubleVector &g)
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
        g[i] = (f2 - f1) / (2.0 * h);
        printf("%14.10f %14.10f %14.10f %14.10f %14.10f\n", f1, f2, f2 - f1, (f2-f1)/(2.0*h), g[i]);
    }
}
