#include "rungekutta.h"
#include <math.h>
#include <float.h>

void RungeKutta::calculate(R2Function *f, double x0, double x1, double y0, double &y1, double dx)
{
    C_UNUSED(f);
    C_UNUSED(x0);
    C_UNUSED(x1);
    C_UNUSED(y0);
    C_UNUSED(y1);
    C_UNUSED(dx);
}

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

void CauchyProblem2::rungeKutta(CauchyProblem2 *cp, double x0 UNUSED_PARAM, double y0 UNUSED_PARAM, double h, unsigned int N, DoubleVector &y)
{
    if (abs(h) < DBL_EPSILON) return;
    y.clear();
    y.resize(N+1);

    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    if (h > 0.0)
    {
        double x = cp->x0;
        y[0] = cp->y0;
        DoubleVector cy(1);
        for (unsigned int i=0; i<=N-1; i++)
        {
            cy[0] = y[i];
            k1 = cp->f(x, cy);
            cy[0] = y[i]+(h/2.0)*k1;
            k2 = cp->f(x+h/2.0, cy);
            cy[0] = y[i]+(h/2.0)*k2;
            k3 = cp->f(x+h/2.0, cy);
            cy[0] = y[i]+h*k2;
            k4 = cp->f(x+h, cy);
            y[i+1] = y[i] + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
            x = x + h;
        }
    }

    //    if (h < 0.0)
    //    {
    //        y[n-1] = y0;
    //        for (int i=n-2; i>=0; i--)
    //        {
    //            k1 = f->fx(x0,        y0);
    //            k2 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k1);
    //            k3 = f->fx(x0+dx/2.0, y0+(dx/2.0)*k2);
    //            k4 = f->fx(x0+dx,     y0+dx*k3);
    //            y0 = y0 + (dx/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
    //            x0 = x0 + dx;
    //            y[i] = y0;
    //        }
    //    }
}

void CauchyProblem2::rungeKutta(std::vector<CauchyProblem2*> cps, double x0, double h, unsigned int N, DoubleMatrix &my)
{
    unsigned int M = cps.size();
    my.clear();
    my.resize(M, N+1);

    DoubleVector k1(M);
    DoubleVector k2(M);
    DoubleVector k3(M);
    DoubleVector k4(M);

    if (h > 0.0)
    {
        for (unsigned int j=0; j<M; j++) my.at(j,0) = cps[j]->y0;

        for (unsigned int i=0; i<N; i++)
        {
            DoubleVector y(M);

            // Calculating k1 vector
            double x = x0;
            for (unsigned int j=0; j<M; j++) y[j] = my[j][i];
            for (unsigned int j=0; j<M; j++) k1[j] = cps[j]->f(x, y);

            // Calculating k2 vector
            x = x0+h/2.0;
            for (unsigned int j=0; j<M; j++) y[j] = my[j][i] + (h/2.0) * k1[j];
            for (unsigned int j=0; j<M; j++) k2[j] = cps[j]->f(x, y);

            // Calculating k3 vector
            x = x0+h/2.0;
            for (unsigned int j=0; j<M; j++) y[j] = my[j][i] + (h/2.0) * k2[j];
            for (unsigned int j=0; j<M; j++) k3[j] = cps[j]->f(x, y);

            // Calculating k4 vector
            x = x0+h;
            for (unsigned int j=0; j<M; j++) y[j] = my[j][i] + h * k3[j];
            for (unsigned int j=0; j<M; j++) k4[j] = cps[j]->f(x, y);

            // Calculating y
            for (unsigned int j=0; j<M; j++) my[j][i+1] = my[j][i] + (h/6.0) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);

            x0 = x0 + h;
        }
    }

    if (h < 0.0)
    {
        for (unsigned int j=0; j<M; j++) my.at(j,N) = cps[j]->y0;

        for (unsigned int i=N; i>=1; i--)
        {
            DoubleVector y(M);

            // Calculating k1 vector
            double x = x0;
            for (unsigned int j=0; j<M; j++) y[j] = my[j][i];
            for (unsigned int j=0; j<M; j++) k1[j] = cps[j]->f(x, y);

            // Calculating k2 vector
            x = x0+h/2.0;
            for (unsigned int j=0; j<M; j++) y[j] = my[j][i] + (h/2.0) * k1[j];
            for (unsigned int j=0; j<M; j++) k2[j] = cps[j]->f(x, y);

            // Calculating k3 vector
            x = x0+h/2.0;
            for (unsigned int j=0; j<M; j++) y[j] = my[j][i] + (h/2.0) * k2[j];
            for (unsigned int j=0; j<M; j++) k3[j] = cps[j]->f(x, y);

            // Calculating k4 vector
            x = x0+h;
            for (unsigned int j=0; j<M; j++) y[j] = my[j][i] + h * k3[j];
            for (unsigned int j=0; j<M; j++) k4[j] = cps[j]->f(x, y);

            // Calculating y
            for (unsigned int j=0; j<M; j++) my[j][i-1] = my[j][i] + (h/6.0) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);

            x0 = x0 + h;
        }
    }
}

void CauchyProblem2::euler1(std::vector<CauchyProblem2 *> cps, double x0, double h, unsigned int N, DoubleMatrix &m)
{
    unsigned int C = cps.size();
    m.clear();
    m.resize(C, N+1);

    for (unsigned int j=0; j<C; j++)
    {
        m[j][0] = cps[j]->y0;
    }

    for (unsigned int i=0; i<N; i++)
    {
        DoubleVector y(C);
        for (unsigned int j=0; j<C; j++) y[j] = m[j][i];
        for (unsigned int j=0; j<C; j++) m[j][i+1] = m[j][i] + h * cps[j]->f(x0, y);
        x0 = x0 + h;
    }
}

void CauchyProblem2::euler2(std::vector<CauchyProblem2 *> cps UNUSED_PARAM, double x0 UNUSED_PARAM, double h UNUSED_PARAM, unsigned int N UNUSED_PARAM, DoubleMatrix &m UNUSED_PARAM)
{
    //    for (unsigned int i=0; i<m.size(); i++) m[i].clear();
    //    m.clear();

    //    unsigned int C = cps.size();
    //    m.resize(C);
    //    DoubleVector k1(C);
    //    DoubleVector k2(C);

    //    for (unsigned int j=0; j<C; j++)
    //    {
    //        m[j].resize(N+1);
    //        m[j][0] = cps[j]->y0;
    //    }

    //    for (unsigned int i=0; i<N; i++)
    //    {
    //        double x = x0;
    //        DoubleVector y(C);
    //        for (unsigned int j=0; j<C; j++) y[j] = m[j][i];
    //        for (unsigned int j=0; j<C; j++) k1[j] = cps[j]->f(x, y);
    //        for (unsigned int j=0; j<C; j++) y[j] = m[j][i]+h*k1[j];
    //        for (unsigned int j=0; j<C; j++) k2[j] = cps[j]->f(x+h, y);
    //        for (unsigned int j=0; j<C; j++) m[j][i+1] = m[j][i] + h/2.0 * (k1[j]+k2[j]);
    //        x0 = x0 + h;
    //    }
}
