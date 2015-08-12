#include "rungekutta.h"
#include <math.h>

void RungeKutta::calculate(R2Function *f, double x0, double x1, double y0, double &y1, double dx)
{}

unsigned int RungeKutta::calculate(R2Function *f, double x0, double x1, double y0, DoubleVector &y, double dx)
{
    if (dx == 0.0) return 0;

    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    int n = (int)(ceil(fabs(x1-x0)/fabs(dx)))+1;
    y.clear();
    y.resize(n);

    if (dx > 0.0)
    {
        y[0] = y0;
        for (int i=1; i<n; i++)
        {
            k1 = f->fx(x0,        y0);
            k2 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k1);
            k3 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k2);
            k4 = f->fx(x0+dx,     y0+dx*k3);

            y0 = y0 + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + dx;
            y[i] = y0;
        }
    }

    if (dx < 0.0)
    {
        y[n-1] = y0;
        for (int i=n-2; i>=0; i--)
        {
            k1 = f->fx(x0,        y0);
            k2 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k1);
            k3 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k2);
            k4 = f->fx(x0+dx,     y0+dx*k3);

            y0 = y0 + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + dx;
            y[i] = y0;
        }
    }

    return n;
}

void RungeKutta::calculate(R2Function *f, double x0, double y0, DoubleVector &y, double dx)
{
    if (dx == 0.0) return;

    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    unsigned int n = y.size();
    if (dx > 0.0) y[0]   = y0;
    if (dx < 0.0) y[n-1] = y0;

    if (dx > 0.0)
    {
        y[0] = y0;
        for (unsigned int i=1; i<n; i++)
        {
            k1 = f->fx(x0,        y0);
            k2 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k1);
            k3 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k2);
            k4 = f->fx(x0+dx,     y0+dx*k3);
            y0 = y0 + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + dx;
            y[i] = y0;
        }
    }

    if (dx < 0.0)
    {
        y[n-1] = y0;
        for (int i=n-2; i>=0; i--)
        {
            k1 = f->fx(x0,        y0);
            k2 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k1);
            k3 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k2);
            k4 = f->fx(x0+dx,     y0+dx*k3);
            y0 = y0 + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + dx;
            y[i] = y0;
        }
    }
}

void RungeKutta::calculate(R2FunctionX f, double x0, double y0, DoubleVector &y, double dx)
{
    if (dx == 0.0) return;

    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    unsigned int n = y.size();
    if (dx > 0.0) y[0]   = y0;
    if (dx < 0.0) y[n-1] = y0;

    if (dx > 0.0)
    {
        y[0] = y0;
        for (unsigned int i=1; i<n; i++)
        {
            k1 = f(x0,        y0);
            k2 = f(x0+dx/2.0, y0+(dx/2.0)*k1);
            k3 = f(x0+dx/2.0, y0+(dx/2.0)*k2);
            k4 = f(x0+dx,     y0+dx*k3);
            y0 = y0 + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + dx;
            y[i] = y0;
        }
    }

    if (dx < 0.0)
    {
        y[n-1] = y0;
        for (int i=n-2; i>=0; i--)
        {
            k1 = f(x0,        y0);
            k2 = f(x0+dx/2.0, y0+(dx/2.0)*k1);
            k3 = f(x0+dx/2.0, y0+(dx/2.0)*k2);
            k4 = f(x0+dx,     y0+dx*k3);
            y0 = y0 + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x0 = x0 + dx;
            y[i] = y0;
        }
    }
}
